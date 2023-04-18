//
// Created by maksbh on 9/12/20.
//

#include <IMGA/IMGA.h>

#include <DendriteUtils.h>

IMGA::IMGA(const DomainExtents &domainExtents, const IBM_METHOD method)
: domainExtents_(domainExtents), method_(method) {
}

IMGA::IMGA(const DomainExtents &domainExtents, const GEOMETRY::Geometry *geometry, const GeomRefinement & geomRefinement, const IBM_METHOD method)
    : domainExtents_(domainExtents),method_(method) {
  geometries_.push_back(geometry);
  m_geomRefinement.push_back(geomRefinement);
}

void IMGA::addGeometry(const GEOMETRY::Geometry *geometry, const GeomRefinement & geomRefinement) {
  geometries_.push_back(geometry);
  m_geomRefinement.push_back(geomRefinement);
}


bool IMGA::ifInside(const DENDRITE_REAL *position) const {

  for (const auto &geom : geometries_) {
    if (geom->ifInside(position)) {
      return true;
    }
  }
  return false;
}

void IMGA::addGaussPoints(std::vector<GaussPoint> &gaussPoints) {
  static constexpr DENDRITE_UINT gaussPointsPerElement = (DIM == 3) ? 3 : 2;
  DENDRITE_UINT numGP = 0;
#if (DIM == 3)
  static const DENDRITE_REAL log4 = log(4);
  for (int geomID = 0; geomID < geometries_.size(); geomID++) {
    const auto &geom = geometries_[geomID];
    const auto &translations = geom->getTranslations();
    const DENDRITE_UINT start = geom->getSTL()->getStart();
    const DENDRITE_UINT finish = geom->getSTL()->getEnd();
    const DENDRITE_UINT sz = finish - start;
    const DENDRITE_REAL octantArea = computeOctantArea(domainExtents_,m_geomRefinement[geomID].octantLevel);
    std::vector<DENDRITE_UINT> numSplit(sz,0);
    const auto &triangles = geom->getSTL()->getTriangles();
    for(int i = start; i < finish; i++){
      const DENDRITE_REAL areaTriangles = computeTriangleArea(triangles[i].triangleCoord);
      const DENDRITE_REAL ratio = areaTriangles/(m_geomRefinement[geomID].ratioArea * octantArea);
      numSplit[i - start] = (std::min(std::max((int)std::ceil(log(ratio)/log4),0),(int)m_geomRefinement[geomID].maxSplitIteration));
    }
    DENDRITE_UINT numSplittedTriangles = 0;
    for(const auto & split:numSplit){
      numSplittedTriangles += (1u << (2*split)); /// 4^x = 2^(2*x)
    }
    numGP += numSplittedTriangles * gaussPointsPerElement * translations.size();
  }


  gaussPoints.resize(numGP);
  GaussPoint GP{};
  DENDRITE_UINT counter = 0;
  std::vector<GEOMETRY::Triangles> splittedTriangles;
  for (DENDRITE_UINT geomID = 0; geomID < geometries_.size(); geomID++) {
    const auto &geom = geometries_[geomID];
    const auto &triangles = geom->getSTL()->getTriangles();
    const auto &translations = geom->getTranslations();
    const DENDRITE_UINT start = geom->getSTL()->getStart();
    const DENDRITE_UINT finish = geom->getSTL()->getEnd();
    GP.geomID = geomID;
    DENDRITE_REAL _pos[DIM];
    const DENDRITE_REAL  octantArea = computeOctantArea(domainExtents_,m_geomRefinement[geomID].octantLevel);
    for (DENDRITE_UINT i = start; i < finish; i++) {
      const auto &tri = triangles[i];
      const DENDRITE_REAL areaTriangles = computeTriangleArea(triangles[i].triangleCoord);
      const DENDRITE_REAL ratio = areaTriangles/(m_geomRefinement[geomID].ratioArea* octantArea);
      const DENDRITE_UINT numSplit = (std::min(std::max((int)std::ceil(log(ratio)/log4),0),(int)m_geomRefinement[geomID].maxSplitIteration));
      GP.elemID = i;
      std::memcpy(GP.normal, tri.normal, sizeof(DENDRITE_REAL) * DIM);

      performTriangleSplitting(splittedTriangles,tri.triangleCoord,numSplit);
      GP.elemArea = areaTriangles / ((1u << (2 * numSplit)) * 1.0); /// Area of splitted triangles
      for(auto &splitTri:splittedTriangles) {
        for (DENDRITE_UINT k = 0; k < 3; k++) { // 3 coords per triangle
          for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
            _pos[dim] = 0.5 * (splitTri.triangleCoord[k][dim] + splitTri.triangleCoord[(k + 1) % 3][dim]);
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
#endif
#if (DIM == 2)
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
      numSplit[i - start] = (std::min(std::max((int) std::ceil(log(ratio) / log2), 0), (int) m_geomRefinement[geomID].maxSplitIteration));
    }
    DENDRITE_UINT numSplittedLines = 0;
    for (const auto &split:numSplit) {
      numSplittedLines += (1u << (split));
    }
    numGP += numSplittedLines * gaussPointsPerElement * translations.size();
  }
  gaussPoints.resize(numGP);
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
      const DENDRITE_REAL lenLines = computeLineLength(lines[i].lineCoord);
      const DENDRITE_REAL ratio = lenLines / (m_geomRefinement[geomID].ratioArea * octantArea);
      const DENDRITE_UINT numSplit = (std::min(std::max((int) std::ceil(log(ratio) / log2), 0), (int) m_geomRefinement[geomID].maxSplitIteration));
      GP.elemID = i;
      std::memcpy(GP.normal, line.normal, sizeof(DENDRITE_REAL) * DIM);
      /// We need to flip the normals if we are solving inside the object.
      if(geom->getRetainSide() == RetainSide::IN){
        for(double & normal : GP.normal){
          normal *=-1;
        }
      }
      performLineSplitting(splittedLines, line.lineCoord, numSplit);
      GP.elemArea = lenLines / ((1u << (numSplit))); /// Length of splitted lines
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
#endif
}

void IMGA::addGaussPoints(std::vector<GaussPoint> & gaussPoints, const std::vector<TALYFEMLIB::ZEROPTV> & positions){
  gaussPoints.resize(positions.size());
  for(int i = 0; i < positions.size(); i++){
    std::memcpy(gaussPoints[i].pos,positions[i].data(),sizeof(double)*DIM);
    gaussPoints[i].elemID = i; // Hacks
    gaussPoints[i].geomID = TALYFEMLIB::GetMPIRank();
  }

}

void IMGA::removeOutOfBoundaryGaussPoints(std::vector<GaussPoint> &gaussPoints){
  std::vector<bool> retainGP(gaussPoints.size(),false);
  for(int i = 0; i < gaussPoints.size(); i++){
    if(gaussPoints[i].isInsideDomain(domainExtents_)){
      retainGP[i] = true;
    }
  }
  bool retainAll = std::all_of(retainGP.begin(),retainGP.end(),[](bool x){return x;}); // all the values are true.
  // Source: https://stackoverflow.com/questions/6506659/is-there-anything-like-stdand-or-stdor

  if(retainAll){
    return;
  }
  /// Step4 : Sort Gauss point according to retain or not retain.
  {
    std::vector<int> permutationIndex(gaussPoints.size(), 0);
    std::fill(permutationIndex.begin(), permutationIndex.end(), 0);
    for (int i = 0; i != permutationIndex.size(); i++) {
      permutationIndex[i] = i;
    }
    std::sort(permutationIndex.begin(), permutationIndex.end(),
              [&](int i, int j) { return (retainGP[i] < retainGP[j]); });
    std::vector<bool> sortedRetainGP(gaussPoints.size());
    std::transform(permutationIndex.begin(), permutationIndex.end(), sortedRetainGP.begin(),
                   [&](int i) { return retainGP[i]; });
    std::swap(retainGP, sortedRetainGP);
    std::vector<GaussPoint> sortedGP(gaussPoints.size());
    std::transform(permutationIndex.begin(), permutationIndex.end(), sortedGP.begin(),
                   [&](int i) { return gaussPoints[i]; });
    std::swap(gaussPoints, sortedGP);
  }


  const auto location = std::find(retainGP.begin(),retainGP.end(),true);
  const int numEntriesToRemove = location - retainGP.begin();
  gaussPoints.erase(gaussPoints.begin(),gaussPoints.begin()+numEntriesToRemove);

}

void  IMGA::addPoints(DA * octDA, const std::vector<TALYFEMLIB::ZEROPTV> & points){
  static_assert((DIM == 3) or (DIM == 2),"IMGA supprorted in 2D or 3D");
  std::vector<GaussPoint> gaussPoints;
  /// Step1 : Fill with Gauss points that each processor owns. At this point the GP is not aligned
  /// with octree partitions.
  addGaussPoints(gaussPoints,points);

  removeOutOfBoundaryGaussPoints(gaussPoints);
  //todo, the gp_normal shown above needs to be oriented consistently.

  /// Step 2 : Compute owner procs for each of the Gauss Points.
  std::vector<TREENODE> treeNodes;
  std::vector<int> ownerRanks;
  computeTreeNodes(gaussPoints, treeNodes);
  computeOwnerRank(octDA,treeNodes, ownerRanks);

  const DENDRITE_UINT numProcs = TALYFEMLIB::GetMPISize();
  std::vector<int> sendProcCount(numProcs, 0);
  MPI_Comm comm = MPI_COMM_WORLD;

  /// Computing how much proc i needs to send to poc j
  for (const auto &rank: ownerRanks) {
    sendProcCount[rank]++;
  }


  std::vector<int> receiveProcCount(numProcs,0);
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
}
void IMGA::fillWithGaussPoints(DA * octDA) {

  static_assert((DIM == 3) or (DIM == 2),"IMGA supprorted in 2D or 3D");
  std::vector<GaussPoint> gaussPoints;
  /// Step1 : Fill with Gauss points that each processor owns. At this point the GP is not aligned
  /// with octree partitions.
  addGaussPoints(gaussPoints);
  removeOutOfBoundaryGaussPoints(gaussPoints);
  //todo, the gp_normal shown above needs to be oriented consistently.

  /// Step 2 : Compute owner procs for each of the Gauss Points.
  std::vector<TREENODE> treeNodes;
  std::vector<int> ownerRanks;
  computeTreeNodes(gaussPoints, treeNodes);
  computeOwnerRank(octDA,treeNodes, ownerRanks);

  const DENDRITE_UINT numProcs = TALYFEMLIB::GetMPISize();
  std::vector<int> sendProcCount(numProcs, 0);
  MPI_Comm comm = MPI_COMM_WORLD;

  /// Computing how much proc i needs to send to poc j
  for (const auto &rank: ownerRanks) {
    sendProcCount[rank]++;
  }


  std::vector<int> receiveProcCount(numProcs,0);
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

  /// At this point each processor will contain GaussPoint aligned with DA partition.

}

void IMGA::computeTreeNodes(const std::vector<GaussPoint> &gaussPoint, std::vector<TREENODE> &treeNodes) const {
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

void IMGA::computeOwnerRank(const DA * octDA, const std::vector<TREENODE> &treeNodes, std::vector<int> &ownerRanks) const {
  ownerRanks.resize(treeNodes.size());
  octDA->computeTreeNodeOwnerProc(treeNodes.data(), treeNodes.size(), ownerRanks.data());

}

void IMGA::findBackGroundElements(const DA * octDA, const std::vector<TREENODE> & treePart) {
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
    m_gpInfo[i].elemID = gaussPoints_[i].elemID;
    m_gpInfo[i].elemArea = gaussPoints_[i].elemArea;
  }
  std::sort(m_gpInfo.begin(), m_gpInfo.end());


  /// Delete the unused vectors
  {
    std::vector<TREENODE> _treeNode;
    std::vector<GaussPoint> _gPoint;
    std::swap(_treeNode, treeNodes);
    std::swap(_gPoint, gaussPoints_);
  }


  const size_t sz = octDA->getTotalNodalSz();
  auto partFront = octDA->getTreePartFront();
  auto partBack = octDA->getTreePartBack();
  const auto tnCoords = octDA->getTNCoords();
  const DENDRITE_UINT nPe = octDA->getNumNodesPerElement();
  const DENDRITE_UINT eleOrder = octDA->getElementOrder();
  OctToPhysical octToPhysical(domainExtents_);
  DENDRITE_REAL *coords = new DENDRITE_REAL[nPe * DIM];
  ot::MatvecBaseCoords<DIM> loop(sz, eleOrder, false, 0, tnCoords, treePart.data(), treePart.size(), *partFront, *partBack);
  DENDRITE_UINT counter = 0;
  DENDRITE_UINT lclElemID = 0;

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
          m_gpInfo[i].localElemID = lclElemID;
#ifndef NDEBUG
          // Checking if the background Element is correctly identified.
          counter++;
          const double *checkPt = m_gpInfo[i].location;
          bool checkX = (coords[0] <= checkPt[0]) and (coords[(nPe - 1) * DIM + 0] >= checkPt[0]);
          bool checkY = (coords[1] <= checkPt[1]) and (coords[(nPe - 1) * DIM + 1] >= checkPt[1]);
          bool checkZ = true;
#if (DIM == 3)
          checkZ = (coords[2] <= checkPt[2]) and (coords[(nPe - 1) * DIM + 2] >= checkPt[2]);
#endif
#endif
          assert(checkX and checkY and checkZ);
        } else {
          ptIdx = i;
          break;
        }
      }
      lclElemID++;
      loop.next();
    } else {
      loop.step();
    }
  }
  delete[] coords;

  MPI_Barrier(MPI_COMM_WORLD);
#ifndef NDEBUG
  assert(counter == m_gpInfo.size());
#endif
  // We can remove duplicate Gauss points.
  // Claim 1: std::unique is valid here because we have already sorted Gauss point
  //          based on the location. So, if tree nodes corresponding to the two nodes
  //          are same. That means they are same Gauss point. But do they encounter FP error?



}
