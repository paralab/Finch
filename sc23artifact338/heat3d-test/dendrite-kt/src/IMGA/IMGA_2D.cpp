#include <IMGA/IMGA_2D.h>
#include <OctToPhysical.h>
#include <sfcTreeLoop_matvec_io.h>
#include <DendriteUtils.h>

IMGA_2D::IMGA_2D(const DA *_octDA,
                 const std::vector<TREENODE> &treePart,
                 const DomainExtents &domainExtents,
                 const GEOMETRY::Geometry *geometry,
                 const std::vector<GeomRefinement> &geomRefinement)
    : domainExtents_(domainExtents), m_octDA(_octDA), m_treePart(treePart), m_geomRefinement(geomRefinement) {
  geometries_.push_back(geometry);
}

void IMGA_2D::addGeometry(const GEOMETRY::Geometry *geometry) {
  geometries_.push_back(geometry);
}

bool IMGA_2D::ifInside(const DENDRITE_REAL *position) const {

  for (const auto &geom : geometries_) {
    if (geom->ifInside(position)) {
      return true;
    }
  }
  return false;
}

void IMGA_2D::fillWithGaussPoints() {
  static_assert(DIM == 2, "IMGA_2D only supported with 2D");
  DENDRITE_UINT numGP = 0;
  static constexpr DENDRITE_UINT gaussPointsPerLine = 2;
  static const DENDRITE_REAL log2 = log(2);
  for (int geomID = 0; geomID < geometries_.size(); geomID++) {
    const auto &geom = geometries_[geomID];
    const auto &translations = geom->getTranslations();
    const DENDRITE_UINT start = geom->getMSH()->getStart();
    const DENDRITE_UINT finish = geom->getMSH()->getEnd();
    const DENDRITE_UINT sz = finish - start;
    const DENDRITE_REAL octantArea = computeOctantArea(domainExtents_, m_geomRefinement[geomID].octantLevel);
    std::vector<DENDRITE_UINT> numSplit(sz, 0);
    const auto &lines = geom->getMSH()->getLines();
    for (int i = start; i < finish; i++) {
      const DENDRITE_REAL lenLines = computeLineLength(lines[i].lineCoord);
      const DENDRITE_REAL ratio = lenLines / (m_geomRefinement[geomID].ratioArea * octantArea);
      numSplit[i] = (std::min(std::max((int) std::ceil(log(ratio) / log2), 0), (int) m_geomRefinement[geomID].maxSplitIteration));
    }
    DENDRITE_UINT numSplittedLines = 0;
    for (const auto &split:numSplit) {
      numSplittedLines += (1u << (split)); /// 4^x = 2^(2*x)
    }
    numGP += numSplittedLines * gaussPointsPerLine * translations.size();
  }

  /// Step1 : Fill with Gauss points that each processor owns. At this point the GP is not aligned
  /// with octree partitions.
  std::vector<GaussPoint> gaussPoints(numGP);
  GaussPoint GP{};
  DENDRITE_UINT counter = 0;
  std::vector<GEOMETRY::Lines> splittedLines;
  for (DENDRITE_UINT geomID = 0; geomID < geometries_.size(); geomID++) {
    const auto &geom = geometries_[geomID];
    const auto &lines = geom->getMSH()->getLines();
    const auto &translations = geom->getTranslations();
    const DENDRITE_UINT start = geom->getMSH()->getStart();
    const DENDRITE_UINT finish = geom->getMSH()->getEnd();
    GP.geomID = geomID;
    DENDRITE_REAL _pos[DIM];
    const DENDRITE_REAL octantArea = computeOctantArea(domainExtents_, m_geomRefinement[geomID].octantLevel);
    for (DENDRITE_UINT i = start; i < finish; i++) {
      const auto &line = lines[i];
      const DENDRITE_REAL lenLines = computeLineLength(lines[start + i].lineCoord);
      const DENDRITE_REAL ratio = lenLines / (m_geomRefinement[geomID].ratioArea * octantArea);
      const DENDRITE_UINT numSplit = (std::min(std::max((int) std::ceil(log(ratio) / log2), 0), (int) m_geomRefinement[geomID].maxSplitIteration));
      GP.lineID = i;
      std::memcpy(GP.normal, line.normal, sizeof(DENDRITE_REAL) * DIM);

      performLineSplitting(splittedLines, line.lineCoord, numSplit);
      GP.lineLength = lenLines / ((1u << (numSplit))); /// Length of splitted lines
      for (auto &splitLine:splittedLines) {
        for (DENDRITE_UINT k = 0; k < 2; k++) { // 2 coords per line
          for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
            gpPositionLines(splitLine.lineCoord, _pos, k);
          }
          for (const auto &trans: translations) {
            for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
              GP.pos[dim] = _pos[dim] + trans.x(dim);
            }
            std::memcpy(&gaussPoints[counter], &GP, sizeof(GaussPoint));
            counter++;
          }
        }
      }

    }
  }
  // todo, the gp_normal have not been calculated
  //TESTGaussPoints("gp_correct" + std::to_string(TALYFEMLIB::GetMPIRank()) + ".csv", gaussPoints, "w");
  /// Step 2 : Compute owner procs for each of the Gauss Points.
  std::vector<TREENODE> treeNodes;
  std::vector<int> ownerRanks;
  computeTreeNodes(gaussPoints, treeNodes);
  computeOwnerRank(treeNodes, ownerRanks);

  const DENDRITE_UINT numProcs = m_octDA->getNpesAll();
  std::vector<int> sendProcCount(numProcs, 0);
  MPI_Comm comm = MPI_COMM_WORLD;

  /// Computing how much proc i needs to send to poc j
  for (const auto &rank: ownerRanks) {
    sendProcCount[rank]++;
  }

  std::vector<int> receiveProcCount(numProcs, 0);
  /// Step3 : Communicating how many Gauss point each have to receive.
  MPI_Alltoall(sendProcCount.data(), 1, MPI_INT, receiveProcCount.data(), 1, MPI_INT, comm);

  /// Step4 : Sort Gauss point according to the owner rank.
  std::vector<int> permutationIndex(gaussPoints.size(), 0);
  {
    std::fill(permutationIndex.begin(), permutationIndex.end(), 0);
    for (int i = 0; i != permutationIndex.size(); i++) {
      permutationIndex[i] = i;
    }
    std::sort(permutationIndex.begin(), permutationIndex.end(),
              [&](int i, int j) { return (ownerRanks[i] < ownerRanks[j]); });
    std::vector<GaussPoint> sortedGP(gaussPoints.size());
    std::transform(permutationIndex.begin(), permutationIndex.end(), sortedGP.begin(),
                   [&](int i) { return gaussPoints[i]; });
    std::swap(gaussPoints, sortedGP);
  }

  std::vector<int> displacementSendData(numProcs, 0);
  std::vector<int> displacementReceiveData(numProcs, 0);

  std::partial_sum(sendProcCount.begin(), std::prev(sendProcCount.end()), std::next(displacementSendData.begin(), 1)); // Prefix sum
  std::partial_sum(receiveProcCount.begin(), std::prev(receiveProcCount.end()), std::next(displacementReceiveData.begin(), 1));

  int totalGPReceived = std::accumulate(receiveProcCount.begin(), receiveProcCount.end(), 0);
  gaussPoints_.resize(totalGPReceived);
  /// Step 5 : All to all communication
  MPI_Alltoallv(gaussPoints.data(), sendProcCount.data(), displacementSendData.data(), GaussPoint::dataType(),
                gaussPoints_.data(), receiveProcCount.data(), displacementReceiveData.data(), GaussPoint::dataType(), comm);
  TESTGaussPoints("gpInfo_" + std::to_string(TALYFEMLIB::GetMPIRank()) + ".csv", gaussPoints_);
  /// At this point each processor will contain GaussPoint aligned with DA partition.

}

void IMGA_2D::computeTreeNodes(const std::vector<GaussPoint> &gaussPoint, std::vector<TREENODE> &treeNodes) const {
  treeNodes.resize(gaussPoint.size());
  const unsigned int maxOctCoord = (1u << (m_uiMaxDepth));
  OctToPhysical octToPhysical(domainExtents_.fullDADomain);
  const Point<DIM> &scalingFactor = octToPhysical.getScalingFactor();
  unsigned int octCoords[DIM];
  const double physToOct[DIM] = {
      (1u << (m_uiMaxDepth)) / scalingFactor.x(0),
      (1u << (m_uiMaxDepth)) / scalingFactor.x(1),
#ifdef ENABLE_3D
      (1u << (m_uiMaxDepth)) / scalingFactor.x(2),
#endif
  };
  std::array<DENDRITE_UINT, DIM> treeNodeArray;
  for (std::size_t i = 0; i < treeNodes.size(); i++) {
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
      // octree coordinate
      octCoords[dim] = (unsigned int) (gaussPoint[i].pos[dim] * physToOct[dim]);
      /// TODO: Don't know what the hell it is.
      // special fix for positive boundaries... (explanation from Hari, May 7, 2018):
      /*
       * Basically, if the point is on the boundary, then the interpolation needs to happen on the face that overlaps
       * with the boundary. This face can be part of 2 elements, one that is inside the domain and one that is outside.
       * By default we always go for the "right" element, but this is not correct for those on the positive boundaries,
       * hence the error. You can fix it when you compute the TreeNode corresponding to the element. You will need a
       * correction like:
       *
       *   if (xint == xyzFac)
       *     xint = xyzFac - 1;
       *   if (yint == xyzFac)
       *     yint = xyzFac -1;
       *   if (zint == xyzFac)
       *     zint = xyzFac -1;
       */

      if (octCoords[dim] == maxOctCoord) {
        octCoords[dim] = maxOctCoord - 1;
      }
    }
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
      treeNodeArray[dim] = octCoords[dim];
    }
    treeNodes[i] = TREENODE(treeNodeArray, m_uiMaxDepth);
  }
}

void IMGA_2D::computeOwnerRank(const std::vector<TREENODE> &treeNodes, std::vector<int> &ownerRanks) const {
  ownerRanks.resize(treeNodes.size());
  m_octDA->computeTreeNodeOwnerProc(treeNodes.data(), treeNodes.size(), ownerRanks.data());

}

void IMGA_2D::findBackGroundElements() {
//  TESTGaussPoints("gp_" + std::to_string(TALYFEMLIB::GetMPIRank()) + ".csv", gaussPoints_);
  /// All the gauss points here belong to the owned processor only.
  const DENDRITE_UINT numEnteries = gaussPoints_.size();
  std::vector<TREENODE> treeNodes(numEnteries);
  computeTreeNodes(gaussPoints_, treeNodes);
  m_gpInfo.resize(gaussPoints_.size());
  for (DENDRITE_UINT i = 0; i < gaussPoints_.size(); i++) {
    m_gpInfo[i].node = treeNodes[i];
    std::memcpy(m_gpInfo[i].location, gaussPoints_[i].pos, sizeof(DENDRITE_REAL) * DIM);
    std::memcpy(m_gpInfo[i].normal, gaussPoints_[i].normal, sizeof(DENDRITE_REAL) * DIM);
    m_gpInfo[i].geomID = gaussPoints_[i].geomID;
    m_gpInfo[i].elemID = gaussPoints_[i].lineID;
    m_gpInfo[i].elemArea = gaussPoints_[i].lineLength;
  }
  std::sort(m_gpInfo.begin(), m_gpInfo.end());
  /// Delete the unused vectors
  {
    std::vector<TREENODE> _treeNode;
    std::vector<GaussPoint> _gPoint;
    std::swap(_treeNode, treeNodes);
    std::swap(_gPoint, gaussPoints_);
  }
#ifndef NDEBUG
  const size_t sz = m_octDA->getTotalNodalSz();
  auto partFront = m_octDA->getTreePartFront();
  auto partBack = m_octDA->getTreePartBack();
  const auto tnCoords = m_octDA->getTNCoords();
  const DENDRITE_UINT nPe = m_octDA->getNumNodesPerElement();
  const DENDRITE_UINT eleOrder = m_octDA->getElementOrder();
  OctToPhysical octToPhysical(domainExtents_);
  DENDRITE_REAL *coords = new DENDRITE_REAL[nPe * DIM];
  ot::MatvecBaseCoords<DIM> loop(sz, eleOrder, false, 0, tnCoords, m_treePart.data(), m_treePart.size(), *partFront, *partBack);

  DENDRITE_UINT counter = 0;

  std::size_t ptIdx = 0;
  while (!loop.isFinished()) {
    if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
      const double *nodeCoordsFlat = loop.subtreeInfo().getNodeCoords();
      const TREENODE octreeNode = loop.subtreeInfo().getCurrentSubtree();
      std::memcpy(coords, nodeCoordsFlat, sizeof(DENDRITE_REAL) * nPe * DIM);
      octToPhysical.convertCoordsToPhys(coords, nPe);
      for (std::size_t i = ptIdx; i < m_gpInfo.size(); i++) {
        const TREENODE &checkNode = m_gpInfo[i].node;
        if ((octreeNode.isAncestor(checkNode)) or (checkNode == octreeNode)) {

          // Checking if the background Element is correctly identified.
          counter++;
          const double *checkPt = m_gpInfo[i].location;
          bool checkX = (coords[0] <= checkPt[0]) and (coords[(nPe - 1) * DIM + 0] >= checkPt[0]);
          bool checkY = (coords[1] <= checkPt[1]) and (coords[(nPe - 1) * DIM + 1] >= checkPt[1]);
          assert(checkX and checkY);

        } else {
          ptIdx = i;
          break;
        }
      }
      loop.next();
    } else {
      loop.step();
    }
  }
  delete[] coords;

  MPI_Barrier(MPI_COMM_WORLD);
  assert(counter == m_gpInfo.size());
#endif
}

void IMGA_2D::gpPositionLines(const DENDRITE_REAL lineCoord[][2], DENDRITE_REAL *_pos, int gp_no) {
  DENDRITE_REAL x1 = lineCoord[0][0];
  DENDRITE_REAL y1 = lineCoord[0][1];
  DENDRITE_REAL x2 = lineCoord[1][0];
  DENDRITE_REAL y2 = lineCoord[1][1];
  if (gp_no == 0) {
    _pos[0] = x1 + (x2 - x1) / 2 * (1 - sqrt(3) / 3);
    _pos[1] = y1 + (y2 - y1) / 2 * (1 - sqrt(3) / 3);
  } else {
    _pos[0] = x1 + (x2 - x1) / 2 * (1 + sqrt(3) / 3);
    _pos[1] = y1 + (y2 - y1) / 2 * (1 + sqrt(3) / 3);
  }
}

#ifndef NDEBUG
void IMGA_2D::TESTGaussPoints(const std::string &filename, const std::vector<GaussPoint> &gaussPoints_test, const std::string mode) {
  FILE *fp = fopen(filename.c_str(), mode.c_str());
  assert(fp != nullptr);
  int rank = TALYFEMLIB::GetMPIRank();

  fprintf(fp, "gp_x,gp_y,gp_nx,gp_ny,geomID,lineID,lineLen,rank\n");
  for (const auto &gp : gaussPoints_test) {
#if (DIM == 2)
    fprintf(fp, "%.10e,%.10e,"
                "%.10e,%.10e,"
                "%1d,%1d,%1d,%1d\n",
            gp.pos[0], gp.pos[1],
            gp.normal[0], gp.normal[1],
            gp.geomID, gp.lineID, gp.lineLength, rank);
#endif
  }
  fclose(fp);
}

#endif
