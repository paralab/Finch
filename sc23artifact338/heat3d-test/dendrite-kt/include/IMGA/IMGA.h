//
// Created by maksbh on 9/12/20.
//

#ifndef DENDRITEKT_IMGA_H
#define DENDRITEKT_IMGA_H

#include <Geometry/Geometry.h>
#include <NodeAndValues.h>
#include <IMGA/IMGAUtils.h>
#include <IMGA/IMGADataTypes.h>
#include <IMGA/SurfaceValues.h>
#include <OctToPhysical.h>
#include <sfcTreeLoop_matvec_io.h>
class IMGA {

  std::vector<const GEOMETRY::Geometry *> geometries_;
  const DomainExtents &domainExtents_;
  std::vector<GeomRefinement> m_geomRefinement;

  std::vector<GaussPoint> gaussPoints_; /// vector with Gauss points
  /**
   * @brief Add Gauss points based on the geometries
   * @param [out] gaussPoints Added Gauss points
   */
  void addGaussPoints(std::vector<GaussPoint> &gaussPoints);
  /**
   * @brief add Gauss Points based on the points
   * @param [out] gaussPoints  Added gauss points
   * @param [in] positions positions corresponding to Gauss points
   */
  void addGaussPoints(std::vector<GaussPoint> &gaussPoints, const std::vector<TALYFEMLIB::ZEROPTV> & positions);

  void removeOutOfBoundaryGaussPoints(std::vector<GaussPoint> &gaussPoints);
  template<int dof>
  void removeOutOfBoundaryGaussPoints(std::vector<GaussPoint> &gaussPoints, std::vector<SurfaceValues<dof>> & surfaceValues){
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

      std::vector<SurfaceValues<dof>> sortedSurfaceValues(gaussPoints.size());
      std::transform(permutationIndex.begin(), permutationIndex.end(), sortedSurfaceValues.begin(),
                     [&](int i) { return surfaceValues[i]; });
      std::swap(surfaceValues,sortedSurfaceValues);
    }

    const auto location = std::find(retainGP.begin(),retainGP.end(),true);
    const int numEntriesToRemove = location - retainGP.begin();
    assert(numEntriesToRemove > 0);
    gaussPoints.erase(gaussPoints.begin(),gaussPoints.begin()+numEntriesToRemove);
    surfaceValues.erase(surfaceValues.begin(),surfaceValues.begin()+numEntriesToRemove);
  }

  const IBM_METHOD method_;
  protected:
  std::vector<NodeAndValues<DENDRITE_REAL> > m_gpInfo; /// Gauss points data

  void findBackGroundElements(const DA * octDA, const std::vector<TREENODE> & treePart);



  void fillWithGaussPoints(DA * octDA);

  void computeTreeNodes(const std::vector<GaussPoint> &gaussPoint, std::vector<TREENODE> &treeNodes) const;

  void computeOwnerRank(const DA * octDA, const std::vector<TREENODE> &treeNodes, std::vector<int> &ownerRanks) const;

 public:
  IMGA(const DomainExtents &domainExtents, const IBM_METHOD method = IBM_METHOD::NITSCHE);
  IMGA(const DomainExtents &domainExtents, const GEOMETRY::Geometry *geometry,  const GeomRefinement & geomRefinement,const IBM_METHOD method = IBM_METHOD::NITSCHE);

  void addGeometry(const GEOMETRY::Geometry *geometry, const GeomRefinement & geomRefinement);
  void addPoints(DA * octDA,const std::vector<TALYFEMLIB::ZEROPTV> & points);
  bool ifInside(const DENDRITE_REAL *position) const;

  void inline initIMGAComputation(DA * octDA, const std::vector<TREENODE> & treePart, const std::vector<TALYFEMLIB::ZEROPTV> & position) {
    addPoints(octDA,position);
    findBackGroundElements(octDA, treePart);
  }

  void inline initIMGAComputation(DA * octDA, const std::vector<TREENODE> & treePart) {
    fillWithGaussPoints(octDA);
    findBackGroundElements(octDA, treePart);
  }


  inline const std::vector<const GEOMETRY::Geometry *> & getGeometries() const{
    return geometries_;
  }
  inline const std::vector<NodeAndValues<DENDRITE_REAL> > & getSurfaceGaussPoints() const{
    return m_gpInfo;
  }
  inline IBM_METHOD getIBMMethod() const{
    return method_;
  }


  template<int dof>
  void findBackGroundElements(const DA * octDA,  const std::vector<TREENODE> & treePart, std::vector<SurfaceValues<dof>> & surfaceValues){
    //  TESTGaussPoints("gp_" + std::to_string(TALYFEMLIB::GetMPIRank()) + ".csv", gaussPoints_);
    /// All the gauss points here belong to the owned processor only.
    const DENDRITE_UINT numEnteries = gaussPoints_.size();
    std::vector<TREENODE> treeNodes(numEnteries);
    computeTreeNodes(gaussPoints_, treeNodes);
    m_gpInfo.clear();
    m_gpInfo.resize(gaussPoints_.size());
    for (DENDRITE_UINT i = 0; i < gaussPoints_.size(); i++) {
      m_gpInfo[i].node = treeNodes[i];
      m_gpInfo[i].localElemID = -1;
      std::memcpy(m_gpInfo[i].location, gaussPoints_[i].pos, sizeof(DENDRITE_REAL) * DIM);
      std::memcpy(m_gpInfo[i].normal, gaussPoints_[i].normal, sizeof(DENDRITE_REAL) * DIM);
      m_gpInfo[i].geomID = gaussPoints_[i].geomID;
      m_gpInfo[i].elemArea = gaussPoints_[i].elemArea;
    }
    /// Delete the unused vectors
    {
      std::vector<TREENODE> _treeNode;
      std::vector<GaussPoint> _gPoint;
      std::swap(_treeNode, treeNodes);
      std::swap(_gPoint, gaussPoints_);
    }

    // Sort
    {
      std::vector<int> permutationIndex(m_gpInfo.size(), 0);
      std::fill(permutationIndex.begin(), permutationIndex.end(), 0);
      for (int i = 0; i != permutationIndex.size(); i++) {
        permutationIndex[i] = i;
      }
      std::sort(permutationIndex.begin(), permutationIndex.end(),
                [&](const int &i, const int &j) { return (m_gpInfo[i] < m_gpInfo[j]); });
      std::vector<NodeAndValues<DENDRITE_REAL> > tempGP(m_gpInfo.size());
      std::transform(permutationIndex.begin(), permutationIndex.end(), tempGP.begin(),
                     [&](int i) { return m_gpInfo[i]; });
      std::swap(m_gpInfo, tempGP);

      std::vector<SurfaceValues<dof>> tempsurfaceValues(surfaceValues.size());
      std::transform(permutationIndex.begin(), permutationIndex.end(), tempsurfaceValues.begin(),
                     [&](int i) { return surfaceValues[i]; });
      std::swap(surfaceValues, tempsurfaceValues);
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
  }



  template<int dof>
  inline void translateGaussPoints(DA * octDA, const std::vector<TREENODE> & treeParition,std::vector<SurfaceValues<dof>> & surfaceValues,
                                   const std::vector<TALYFEMLIB::ZEROPTV> &translations,
                                   const std::vector<TALYFEMLIB::ZEROPTV> &centers_of_mass,
                                   const std::vector<std::vector<TALYFEMLIB::ZEROPTV>> &RotMats) {
    assert(surfaceValues.size() == m_gpInfo.size());
    gaussPoints_.resize(m_gpInfo.size());
    for(int i = 0; i < gaussPoints_.size(); i++){
      std::memcpy(gaussPoints_[i].pos,m_gpInfo[i].location,sizeof(DENDRITE_REAL)*DIM);
      std::memcpy(gaussPoints_[i].normal,m_gpInfo[i].normal,sizeof(DENDRITE_REAL)*DIM);
      gaussPoints_[i].elemID = m_gpInfo[i].elemID;
      gaussPoints_[i].elemArea = m_gpInfo[i].elemArea;
      gaussPoints_[i].geomID = m_gpInfo[i].geomID;
    }
    // rotation
    for (auto &gp: gaussPoints_) {
      TALYFEMLIB::ZEROPTV newLocation;
      for (int dim = 0; dim < DIM; dim++) {
        for (int k = 0; k < DIM; k++) {
          newLocation(dim) += (gp.pos[k] - centers_of_mass[gp.geomID](k))*RotMats[gp.geomID][dim](k);
        }
      }
      for (int dim = 0; dim < DIM; dim++) {
        gp.pos[dim] = newLocation[dim] + centers_of_mass[gp.geomID][dim];
      }
    }
    // translation
    for (auto &gp: gaussPoints_) {
      for (int d = 0; d < DIM; d++) {
        gp.pos[d] += translations[gp.geomID].data()[d];
      }
    }
    removeOutOfBoundaryGaussPoints(gaussPoints_);
    std::vector<TREENODE> treeNodes;
    std::vector<int> ownerRanks;
    computeTreeNodes(gaussPoints_, treeNodes);
    computeOwnerRank(octDA, treeNodes, ownerRanks);


    int numProcs = octDA->getNpesAll();
    std::vector<int> sendProcCount(numProcs, 0);
    MPI_Comm comm = MPI_COMM_WORLD;

    /// Computing how much proc i needs to send to poc j
    for (const auto &rank: ownerRanks) {
      sendProcCount[rank]++;
    }
    std::vector<int> receiveProcCount(numProcs, 0);
    /// Step3 : Communicating how many Gauss point each have to receive.
    MPI_Alltoall(sendProcCount.data(), 1, MPI_INT, receiveProcCount.data(), 1, MPI_INT, comm);

    {
      std::vector<int> permutationIndex(gaussPoints_.size(), 0);
      std::fill(permutationIndex.begin(), permutationIndex.end(), 0);
      for (int i = 0; i != permutationIndex.size(); i++) {
        permutationIndex[i] = i;
      }
      std::sort(permutationIndex.begin(), permutationIndex.end(),
                [&](int i, int j) { return (ownerRanks[i] < ownerRanks[j]); });
      std::vector<GaussPoint> sortedGP(gaussPoints_.size());
      std::transform(permutationIndex.begin(), permutationIndex.end(), sortedGP.begin(),
                     [&](int i) { return gaussPoints_[i]; });

      std::vector<SurfaceValues<dof>> sortedSurfaceValues(surfaceValues.size());

      std::transform(permutationIndex.begin(), permutationIndex.end(), sortedSurfaceValues.begin(),
                     [&](int i) { return surfaceValues[i]; });

      std::swap(gaussPoints_, sortedGP);
      std::swap(surfaceValues, sortedSurfaceValues);

      // Now the Gauss Points and surfaceValues are sorted based on owner rank
    }


    std::vector<int> sendCounts(octDA->getNpesAll(),0);
    // Now lets do P2P communication.


    for(const auto & ownRank:ownerRanks){
      sendCounts[ownRank]++;
    }

    std::vector<int> receiveCount(octDA->getNpesAll());
    MPI_Alltoall(sendCounts.data(), 1, MPI_INT, receiveCount.data(), 1, MPI_INT, octDA->getCommActive());

    std::vector<int> recvDisplacements(receiveCount.size(), 0);
    std::vector<int> sendDisplacements(sendCounts.size(), 0);

    for (int i = 1; i < sendDisplacements.size(); i++) {
      sendDisplacements[i] = sendDisplacements[i - 1] + sendCounts[i - 1];
      recvDisplacements[i] = recvDisplacements[i - 1] + receiveCount[i - 1];
    }

    std::vector<MPI_Request *> recieveRequests;
    int totalGaussPointSize = std::accumulate(receiveCount.begin(), receiveCount.end(), 0);
    std::vector<GaussPoint> movedGaussPoint(totalGaussPointSize);
    std::vector<SurfaceValues<dof>> movedSurfaceValues(totalGaussPointSize);
    // Start receiving communication
    for (int i = 0; i < receiveCount.size(); i++) {
      if ((receiveCount[i] > 0) and (i != octDA->getRankAll())) {
        MPI_Request *req1 = new MPI_Request();
        recieveRequests.push_back(req1);
        MPI_Irecv(&movedGaussPoint[recvDisplacements[i]], receiveCount[i], GaussPoint::dataType(), i, MPI_ANY_TAG,
                  octDA->getCommActive(), req1);

        MPI_Request *req2 = new MPI_Request();
        recieveRequests.push_back(req2);
        MPI_Irecv(&movedSurfaceValues[recvDisplacements[i]], receiveCount[i], SurfaceValues<dof>::dataType(), i, MPI_ANY_TAG,
                  octDA->getCommActive(), req2);

      }
    }
    // Start sending communication
    for (int i = 0; i < receiveCount.size(); i++) {
      if ((sendCounts[i] > 0) and (i != octDA->getRankAll())) {
        MPI_Send(&gaussPoints_[sendDisplacements[i]], sendCounts[i], GaussPoint::dataType(), i,
                 octDA->getRankAll(), octDA->getCommActive());
        MPI_Send(&surfaceValues[sendDisplacements[i]], sendCounts[i], SurfaceValues<dof>::dataType(), i,
                 octDA->getRankAll(), octDA->getCommActive());

      }
    }


    // Copy own data that
    int rank = octDA->getRankAll();
    assert(sendCounts[rank] == receiveCount[rank]);

    std::memcpy(&movedGaussPoint[recvDisplacements[rank]], &gaussPoints_[sendDisplacements[rank]],
                sizeof(GaussPoint) * sendCounts[rank]);
    std::memcpy(&movedSurfaceValues[recvDisplacements[rank]], &surfaceValues[sendDisplacements[rank]],
                sizeof(SurfaceValues<dof>) * sendCounts[rank]);

    for(int i = 0; i < recieveRequests.size(); i++){
      MPI_Wait(recieveRequests[i],MPI_STATUS_IGNORE);
    }


    for(int i = 0; i < recieveRequests.size(); i++){
      delete recieveRequests[i];
    }

    std::swap(gaussPoints_,movedGaussPoint);
    std::swap(movedSurfaceValues,surfaceValues);
    assert(surfaceValues.size()==gaussPoints_.size());
    // Now at this moment, Gauss point and surface values are at moved according to the new position.

    this->findBackGroundElements(octDA,treeParition,surfaceValues);



  }
};

#endif //DENDRITEKT_IMGA_H
