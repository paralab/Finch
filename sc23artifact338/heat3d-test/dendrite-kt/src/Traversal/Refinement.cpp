//
// Created by maksbh on 6/13/20.
//

#include <Traversal/Refinement.h>
#include <DendriteUtils.h>
#include <intergridTransfer.h>
#include <PETSc/PetscUtils.h>
Refinement::Refinement(DA *octDA, const std::vector<TREENODE> &treePart, const DomainExtents &domainInfo)
  : Traversal(octDA, treePart, domainInfo) {
  assert(octDA->getLocalElementSz() == treePart.size());
  refineFlags_.resize(treePart.size(), ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
  m_refinementType = RefinementType::POSITION_BASED;
  coords_.resize(octDA->getNumNodesPerElement());

}
Refinement::Refinement(DA * octDA, const std::vector<TREENODE> & treePart, const VecInfo & v, const DomainExtents & domainInfo)
  :Traversal(octDA,treePart,v,domainInfo){
  assert(octDA->getLocalElementSz() == treePart.size());
  refineFlags_.resize(treePart.size(), ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
  m_refinementType = RefinementType::VALUE_BASED;
  coords_.resize(octDA->getNumNodesPerElement());
}

void Refinement::traverseOperation(TALYFEMLIB::FEMElm &fe) {
  assert(counter_ < refineFlags_.size());
  coordsToZeroptv(m_coords, coords_, this->m_octDA->getElementOrder(), true);
  refineFlags_[counter_] = getRefineFlags(fe, coords_);
  counter_++;
}

void Refinement::traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar * values) {
  assert(counter_ < refineFlags_.size());
  coordsToZeroptv(m_coords, coords_, this->m_octDA->getElementOrder(), true);
  refineFlags_[counter_] = getRefineFlags(fe, coords_,values);
  counter_++;
}

DA *Refinement::getRefineDA(std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> &oldDAtreeNode,const double sfcTol) {
  throw  std::logic_error("Use getRefineSubDA");
  return nullptr;
//  bool doRefine = true;
//  if (std::all_of(refineFlags_.begin(), refineFlags_.end(),
//                  [](ot::OCT_FLAGS::Refine v) { return v == ot::OCT_FLAGS::Refine::OCT_NO_CHANGE; })) {
//    doRefine = false;
//  }
//
//  bool gRefine = false;
//  MPI_Allreduce(&doRefine, &gRefine, 1, MPI_CXX_BOOL, MPI_LOR, this->m_octDA->getCommActive());
//
//  if (not(gRefine)) {
//    return NULL;
//  }
//  std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> newDATreeNode;
//  ot::SFC_Tree<DENDRITE_UINT, DIM>::distRemeshWholeDomain(oldDAtreeNode, refineFlags_, newDATreeNode,
//                                                          sfcTol, MPI_COMM_WORLD);
//  newDATreeNode_.resize(newDATreeNode.size());
//  std::copy(newDATreeNode.begin(), newDATreeNode.end(), newDATreeNode_.begin());
//  DA *newDA = new DA(newDATreeNode, MPI_COMM_WORLD, this->m_octDA->getElementOrder(), 100, sfcTol);
//  std::swap(newDATreeNode_, oldDAtreeNode);
//  return newDA;
}

DA *Refinement::getRefineSubDA(DistTREE &oldDistTree,const double sfcTol) {
  bool doRefine = true;
  if (std::all_of(refineFlags_.begin(), refineFlags_.end(),
                  [](ot::OCT_FLAGS::Refine v) { return v == ot::OCT_FLAGS::Refine::OCT_NO_CHANGE; })) {
    doRefine = false;
  }

  bool gRefine = false;
  MPI_Allreduce(&doRefine, &gRefine, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);

  if (not(gRefine)) {
    TALYFEMLIB::PrintStatus("Active Processors [No Remesh] = ", m_octDA->getNpesActive());
    return NULL;
  }
  DistTREE::distRemeshSubdomain(oldDistTree, refineFlags_, newDistTree_,surrDistTree_,ot::RemeshPartition::SurrogateOutByIn,sfcTol);
  DENDRITE_UINT oldDASz = m_octDA->getLocalElementSz();
  DENDRITE_UINT newDASz = newDistTree_.getFilteredTreePartSz();
  doRefine = (oldDASz!=newDASz);
  MPI_Allreduce(&doRefine, &gRefine, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
  if (not(gRefine)) {
    TALYFEMLIB::PrintStatus("Active Processors [No Remesh] = ", m_octDA->getNpesActive());
    return NULL;
  }
  DA *newDA = new DA(newDistTree_, MPI_COMM_WORLD, this->m_octDA->getElementOrder(), 100, sfcTol);
  std::swap(newDistTree_, oldDistTree);
  TALYFEMLIB::PrintStatus("Active Processors [Remesh] = ", newDA->getNpesActive());
  return newDA;


}

DA *Refinement::getForceRefineSubDA(DistTREE &oldDistTree, const double sfcTol) {
  std::fill(refineFlags_.begin(), refineFlags_.end(), ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
  DistTREE::distRemeshSubdomain(oldDistTree, refineFlags_, newDistTree_, surrDistTree_, ot::RemeshPartition::SurrogateOutByIn,sfcTol);

  DA *newDA = new DA(newDistTree_, MPI_COMM_WORLD, this->m_octDA->getElementOrder(), 100, sfcTol);
  TALYFEMLIB::PrintStatus("Active Processors [Forced] = ", newDA->getNpesActive());
  std::swap(newDistTree_,oldDistTree);
  return newDA;


}

void Refinement::initRefinement() {
  this->traverse();
}

void Refinement::petscIntergridTransfer(DA *newDA, const DistTREE &newDistTree, std::vector<VecInfo> & vec){
  if(vec.size() == 1){
    petscIntergridTransfer(newDA,newDistTree,vec[0].v,vec[0].ndof);
    return;
  }
  int ndof = 0;
  for(const auto & v: vec){
    ndof+= v.ndof;
  }
  Vec outVec;
  VecInfo vecInfoOutVec(outVec,ndof,0);

  if(m_octDA->isActive()) {
    PETSc::recombineVec(m_octDA, vec, vecInfoOutVec, false);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  petscIntergridTransfer(newDA,newDistTree,vecInfoOutVec.v,vecInfoOutVec.ndof);

  if(m_octDA->isActive()) {
    for(auto & v: vec){
      VecDestroy(&v.v);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for(auto & v: vec){
    newDA->petscCreateVector(v.v, false, false,v.ndof);
  }
  if(newDA->isActive()){
    PETSc::scatterVec(newDA,vec,vecInfoOutVec,true);
    VecDestroy(&vecInfoOutVec.v);
  }


}

void Refinement::petscIntergridTransfer(DA *newDA, const DistTREE &newDistTree, std::vector<VecUtils> & vec){
  if(vec.size() == 1){
    petscIntergridTransfer(newDA,newDistTree,vec[0].petscVec,vec[0].ndof);
    return;
  }
  int ndof = 0;
  for(const auto & v: vec){
    ndof+= v.ndof;
  }
  Vec outVec;
  VecUtils vecInfoOutVec(outVec,ndof);

  if(m_octDA->isActive()) {
    PETSc::recombineVec(m_octDA, vec, vecInfoOutVec, false);
  }

  petscIntergridTransfer(newDA,newDistTree,vecInfoOutVec.petscVec,vecInfoOutVec.ndof);

  if(m_octDA->isActive()) {
    for(auto & v: vec){
      VecDestroy(&v.petscVec);
    }
  }
  for(auto & v: vec){
    newDA->petscCreateVector(v.petscVec, false, false,v.ndof);
  }
  if(newDA->isActive()){
    PETSc::scatterVec(newDA,vec,vecInfoOutVec, true);
    VecDestroy(&vecInfoOutVec.petscVec);
  }


}

void Refinement::petscIntergridTransfer(DA *newDA, const DistTREE &newDistTree, Vec &inVec, const DENDRITE_UINT ndof) {
  static std::vector<VECType> fineGhosted, surrGhosted;
  newDA->template createVector<VECType>(fineGhosted, false, true, ndof);
  surrDA_->template createVector<VECType>(surrGhosted, false, true, ndof);
  std::fill(fineGhosted.begin(), fineGhosted.end(), 0);

  VECType *fineGhostedPtr = fineGhosted.data();
  VECType *surrGhostedPtr = surrGhosted.data();
  PetscScalar *array = nullptr;
  if(m_octDA->isActive()) {
    VecGetArray(inVec, &array);
  }

  ot::distShiftNodes(*this->m_octDA, array,
                     *surrDA_, surrGhostedPtr + ndof * surrDA_->getLocalNodeBegin(),
                     ndof);

  surrDA_->template readFromGhostBegin<VECType>(surrGhostedPtr, ndof);
  surrDA_->template readFromGhostEnd<VECType>(surrGhostedPtr, ndof);

  fem::MeshFreeInputContext<VECType, TREENODE>
    inctx{surrGhostedPtr,
          surrDA_->getTNCoords(),
          (unsigned) surrDA_->getTotalNodalSz(),
          &(*surrDistTree_.getTreePartFiltered().cbegin()),
          surrDistTree_.getTreePartFiltered().size(),
          *surrDA_->getTreePartFront(),
          *surrDA_->getTreePartBack()};

  fem::MeshFreeOutputContext<VECType, TREENODE>
    outctx{fineGhostedPtr,
           newDA->getTNCoords(),
           (unsigned) newDA->getTotalNodalSz(),
           &(*newDistTree.getTreePartFiltered().cbegin()),
           newDistTree.getTreePartFiltered().size(),
           *newDA->getTreePartFront(),
           *newDA->getTreePartBack()};

  // Hack for the current version
  std::vector<char> outDirty;
  newDA->template createVector<char>(outDirty, false, true, 1);
  newDA->setVectorByScalar(&(*outDirty.begin()), (std::array<char, 1>{0}).data(), false, true, 1);

  const RefElement *refel = newDA->getReferenceElement();
  fem::locIntergridTransfer(inctx, outctx, ndof, refel, &(*outDirty.begin()));

  newDA->template writeToGhostsBegin<VECType>(fineGhostedPtr, ndof, &(*outDirty.cbegin()));
  newDA->template writeToGhostsEnd<VECType>(fineGhostedPtr, ndof, false, &(*outDirty.cbegin()));
  if(m_octDA->isActive()) {
    VecRestoreArray(inVec, &array);
    VecDestroy(&inVec);
  }

  newDA->petscCreateVector(inVec, false, false, ndof);
  double *newDAVecArray;
  if(newDA->isActive()) {
    VecGetArray(inVec, &newDAVecArray);
  }
  newDA->template ghostedNodalToNodalVec<VECType>(fineGhostedPtr, newDAVecArray, true, ndof);
  if(newDA->isActive()) {
    VecRestoreArray(inVec, &newDAVecArray);
  }
}

void Refinement::initPetscIntergridTransfer() {
  surrDA_ = new DA(surrDistTree_, MPI_COMM_WORLD, this->m_octDA->getElementOrder());
}

void Refinement::finializeIntergridTransfer(){
  delete surrDA_;
}