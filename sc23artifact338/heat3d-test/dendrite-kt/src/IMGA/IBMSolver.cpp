////
//// Created by maksbh on 7/12/20.
////
//
//#include <IMGA/IBMSolver.h>
//#include <sfcTreeLoop_matvec_io.h>
//
//IBMSolver::IBMSolver(const DA* octDA, const STL::Geometry * geometry, const DomainInfo & domainInfo)
//:m_geometry(geometry),m_octDA(octDA),m_domain(domainInfo){
//  fillAndPartitionGaussPoints();
//  getBackGroundElement();
//}
//
//void IBMSolver::fillGaussPointsWithoutCommunication(std::vector<Point<DIM>> & points){
//  /// TODO:
//}
//void IBMSolver::fillGaussPointsWithCommunication(std::vector<Point<DIM>> & points){
//  /// TODO:
//}
//
//void IBMSolver::fillAndPartitionGaussPoints() {
//  // Stage 1: Filling with Gauss Points
//  // There are two options:
//  // 1. Every processor fills in the Gauss points, and keeps only the local ones.
//  // Advantage: No communication
//  // Disadvantage: Can be time consuming
//  // Function to be implemented: fillGaussPointWithoutCommunication
//  // 2. Every processor fills in only the set of triangles and compute the owner proc.
//  // Advantage: Very few operation
//  // Disadvantage: Communication
//  // Function to be implemented: fillGaussPointWithCommunication
//
//  ///@songhze: TODO : Implement two separate functions for both of them.
//
//
//  // Once you have the gauss Point, below  needs to be changed accordingly
//
//  // Allocating points : Need to update with Gauss points.
//  static constexpr DENDRITE_UINT numEnteries = 1000;
//  std::vector<Point<DIM>> pointsData(numEnteries);
//  std::array<DENDRITE_REAL,DIM> value;
//  for(int i = 0; i < numEnteries; i++){
//    for(int d = 0; d < DIM; d++){
//      value[d] = static_cast<DENDRITE_REAL>(rand()/(RAND_MAX*1.0));
//    }
//    pointsData[i] = Point<DIM>(value);
//  }
//
//  // Computing owner ranks. This part will remain same.
//  std::vector<TREENODE> treeNodes(numEnteries);
//  std::vector<int> ownerRank(numEnteries);
//  computeTreeNodesAndCordinates(m_domain,pointsData,treeNodes);
//  m_octDA->computeTreeNodeOwnerProc(treeNodes.data(),treeNodes.size(),ownerRank.data());
//
//  DENDRITE_UINT  rank = m_octDA->getRankActive();
//  std::size_t  myNumGaussPoint = std::count(ownerRank.begin(),ownerRank.end(),rank);
//  gaussPoints_.resize(myNumGaussPoint);
//
//  std::size_t counter = 0;
//  for(std::size_t i = 0; i < numEnteries; i++){
//    if(ownerRank[i] == rank){
//      gaussPoints_[counter].node = treeNodes[i];
//      gaussPoints_[counter].geomID = i;
//      for(DENDRITE_UINT dim = 0; dim < DIM; dim++){
//        gaussPoints_[counter].values[dim] = pointsData[i].x(dim);
//      }
//      counter++;
//    }
//  }
//  assert(counter == myNumGaussPoint);
//}
//
//
//void IBMSolver::getBackGroundElement() {
//  std::sort(gaussPoints_.begin(),gaussPoints_.end());
//
//  const size_t sz = m_octDA->getTotalNodalSz();
//  auto partFront = m_octDA->getTreePartFront();
//  auto partBack = m_octDA->getTreePartBack();
//  const auto tnCoords = m_octDA->getTNCoords();
//  const DENDRITE_UINT nPe = m_octDA->getNumNodesPerElement();
//  const DENDRITE_UINT eleOrder = m_octDA->getElementOrder();
//  OctToPhysical octToPhysical(m_domain);
//  DENDRITE_REAL *coords = new DENDRITE_REAL[nPe*DIM];
//  ot::MatvecBaseCoords <DIM> loop(sz,eleOrder, false,0,tnCoords,*partFront,*partBack);
//#ifndef NDEBUG
//  DENDRITE_UINT counter = 0;
//#endif
//  std::size_t ptIdx = 0;
//  DENDRITE_UINT lclEleID = 0;
//  while(!loop.isFinished()){
//    if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
//      const double *nodeCoordsFlat = loop.subtreeInfo().getNodeCoords();
//      const TREENODE  octreeNode = loop.subtreeInfo().getCurrentSubtree();
//      std::memcpy(coords,nodeCoordsFlat, sizeof(DENDRITE_REAL)*nPe*DIM);
//      octToPhysical.convertCoordsToPhys(coords,nPe);
//      for(std::size_t i = ptIdx; i < gaussPoints_.size(); i++){
//        const TREENODE &checkNode = gaussPoints_[i].node;
//        if ((octreeNode.isAncestor(checkNode)) or (checkNode == octreeNode)){
//          gaussPoints_[i].backGroundElementID = lclEleID;
//#ifndef NDEBUG
//          // Checking if the background Element is correctly identified.
//          counter++;
//          const double *checkPt = gaussPoints_[i].values;
//          bool checkX = (coords[0] <= checkPt[0]) and (coords[(nPe-1)*DIM + 0] >= checkPt[0]);
//          bool checkY = (coords[1] <= checkPt[1]) and (coords[(nPe-1)*DIM + 1] >= checkPt[1]);
//          bool checkZ = (coords[2] <= checkPt[2]) and (coords[(nPe-1)*DIM + 2] >= checkPt[2]);
//          assert(checkX and checkY and checkZ);
//#endif
//        }
//        else{
//          lclEleID++;
//          ptIdx = i;
//          break;
//        }
//      }
//      loop.next();
//    }
//    else{
//      loop.step();
//    }
//  }
//  delete[] coords;
//
//#ifndef NDEBUG
//  assert(counter == gaussPoints_.size());
//#endif
//}
//
//const std::vector<NodeAndValues<DENDRITE_REAL,DIM>> & IBMSolver::getGaussPoints() const{
//  return gaussPoints_;
//}
//
//void IBMSolver::computeTreeNodesAndCordinates(const DomainInfo & domainInfo, const std::vector<Point<DIM>> & points,std::vector<TREENODE>& treeNodes) {
//    assert(points.size() == treeNodes.size());
//    const unsigned int maxOctCoord = (1u << (m_uiMaxDepth));
//    OctToPhysical octToPhysical(domainInfo);
//    const Point<DIM> &scalingFactor = octToPhysical.getScalingFactor();
//    unsigned int octCoords[DIM];
//    const double physToOct[DIM] = {
//        (1u << (m_uiMaxDepth)) / scalingFactor.x(0),
//        (1u << (m_uiMaxDepth)) / scalingFactor.x(1),
//#ifdef ENABLE_3D
//        (1u << (m_uiMaxDepth)) / scalingFactor.x(2),
//#endif
//    };
//    std::array<DENDRITE_UINT,DIM> treeNodeArray;
//    for (std::size_t i = 0; i < points.size(); i++) {
//      for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
//        // octree coordinate
//        octCoords[dim] = (unsigned int) (points[i].x(dim) * physToOct[dim]);
//
//        /// TODO: Don't know what the hell it is.
//        // special fix for positive boundaries... (explanation from Hari, May 7, 2018):
//        /*
//         * Basically, if the point is on the boundary, then the interpolation needs to happen on the face that overlaps
//         * with the boundary. This face can be part of 2 elements, one that is inside the domain and one that is outside.
//         * By default we always go for the "right" element, but this is not correct for those on the positive boundaries,
//         * hence the error. You can fix it when you compute the TreeNode corresponding to the element. You will need a
//         * correction like:
//         *
//         *   if (xint == xyzFac)
//         *     xint = xyzFac - 1;
//         *   if (yint == xyzFac)
//         *     yint = xyzFac -1;
//         *   if (zint == xyzFac)
//         *     zint = xyzFac -1;
//         */
//
//        if (octCoords[dim] == maxOctCoord) {
//          octCoords[dim] = maxOctCoord - 1;
//        }
//      }
//      for(DENDRITE_UINT dim = 0; dim <DIM; dim++){
//        treeNodeArray[dim] = octCoords[dim];
//      }
//      treeNodes[i] = TREENODE(treeNodeArray,m_uiMaxDepth);
//    }
//}