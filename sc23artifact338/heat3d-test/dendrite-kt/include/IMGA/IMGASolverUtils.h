//
// Created by maksbh on 7/29/21.
//

#ifndef DENDRITEKT_IMGASOLVERUTILS_H
#define DENDRITEKT_IMGASOLVERUTILS_H

#include <IMGA/Marker.h>
static void getIBMDirichletNodes(const DA * octDA, const std::vector<TREENODE> &treeNodes,std::vector<PetscInt> & dirichletNodes, const Marker * marker){
  if(octDA->isActive()) {
    dirichletNodes.clear();
    const auto & markers = marker->getMarkers();

    std::vector<PetscInt> nonDirichletNodes;
    using CoordT = typename ot::DA<DIM>::C;
    using ot::RankI;

    const size_t ghostedNodalSz = octDA->getTotalNodalSz();
    const TREENODE *odaCoords = octDA->getTNCoords();
    const std::vector<RankI> &ghostedGlobalNodeId = octDA->getNodeLocalToGlobalMap();
    const bool visitEmpty = false;
    const unsigned int padLevel = 0;

    ot::MatvecBaseIn<DIM, RankI, false> treeLoopIn(ghostedNodalSz,
                                                   1,                // node id is scalar
                                                   octDA->getElementOrder(),
                                                   visitEmpty,
                                                   padLevel,
                                                   odaCoords,
                                                   &(*ghostedGlobalNodeId.cbegin()),
                                                   &(*treeNodes.cbegin()),
                                                   treeNodes.size(),
                                                   *octDA->getTreePartFront(),
                                                   *octDA->getTreePartBack());

    int eleCounter = 0;
    while (!treeLoopIn.isFinished()) {
      const ot::TreeNode<CoordT, DIM> subtree = treeLoopIn.getCurrentSubtree();
      const auto subtreeInfo = treeLoopIn.subtreeInfo();
      if (treeLoopIn.isPre() && subtreeInfo.isLeaf()) {
        const RankI *nodeIdsFlat = subtreeInfo.readNodeValsIn();
        const std::vector<bool>& nodeNonhangingIn = subtreeInfo.readNodeNonhangingIn();
        assert((markers[eleCounter] > ElementMarker::INTERCEPTED_ELEMENT) and (markers[eleCounter] <=ElementMarker::INTERCEPTED_GP));
        if(markers[eleCounter] != IN_GP){
          const DENDRITE_UINT nPe = octDA->getNumNodesPerElement();
          for(int i = 0; i < nPe ; i++){
            assert(nodeNonhangingIn.size() == nPe);
              nonDirichletNodes.push_back(nodeIdsFlat[i]);
          }
        }
        eleCounter++;
      }
      treeLoopIn.step();
    }
    assert(eleCounter == markers.size());
    std::sort(nonDirichletNodes.begin(),nonDirichletNodes.end());
    auto last = std::unique(nonDirichletNodes.begin(), nonDirichletNodes.end());
    nonDirichletNodes.erase(last,nonDirichletNodes.end());
    Vec nonDirichletVec;
    octDA->petscCreateVector(nonDirichletVec,false,false,1);
    VecSet(nonDirichletVec,0.0);
    std::vector<PetscScalar> values(nonDirichletNodes.size(),1.0);
    VecSetValues(nonDirichletVec,nonDirichletNodes.size(),nonDirichletNodes.data(),values.data(),ADD_VALUES);
    VecAssemblyBegin(nonDirichletVec);
    VecAssemblyEnd(nonDirichletVec);

    PetscInt start, end;
    VecGetOwnershipRange(nonDirichletVec,&start,&end);
    const PetscScalar  * array;
    VecGetArrayRead(nonDirichletVec,&array);
    for(int nodeID = start; nodeID < end; nodeID++){
      if(FEQUALS(array[nodeID - start],0.0)){
        dirichletNodes.push_back(nodeID);
      }
    }
    VecRestoreArrayRead(nonDirichletVec,&array);
    VecDestroy(&nonDirichletVec);
  }

}
#endif //DENDRITEKT_IMGASOLVERUTILS_H
