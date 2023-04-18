//
// Created by maksbh on 5/9/20.
//

#include <Traversal/Traversal.h>
#include <DendriteUtils.h>

Traversal::Traversal(ot::DA<DIM> * octDA,const std::vector<TREENODE> & treePart, const DomainExtents & domain)
:m_octToPhysical(domain),m_treePart(treePart){
  m_octDA = octDA;
  m_traversalType = TRAVERSAL_TYPE::COORDS;
  init();

}

Traversal::Traversal(ot::DA<DIM> * octDA, const std::vector<TREENODE> & treePart, const VecInfo & v,const DomainExtents & domain)
:m_octToPhysical(domain),m_treePart(treePart){
  m_octDA= octDA;
  v_ = v;
  m_traversalType = TRAVERSAL_TYPE::VALUES;
  init();
}

void Traversal::constructTalyGrid(const double *coords){
  if(eleOrder_ == 1) {
    sync_.syncCoords<1>(coords, &m_grid);
    return;
  }
  if(eleOrder_ == 2) {
    sync_.syncCoords<2>(coords, &m_grid);
    return;
  }
  throw  TALYFEMLIB::TALYException() << "NotSupported in func" << " " << __func__ << "for eleOrder = " << eleOrder_ << "\n" ;

}

bool Traversal::traverseByCoords() {
  const size_t sz = m_octDA->getTotalNodalSz();
  auto partFront = m_octDA->getTreePartFront();
  auto partBack = m_octDA->getTreePartBack();
  const auto tnCoords = m_octDA->getTNCoords();
  int elemID=0;

  ot::MatvecBaseCoords <DIM> loop(sz,eleOrder_, false,0,tnCoords,&(*m_treePart.cbegin()),m_treePart.size(),*partFront,*partBack);
  while(!loop.isFinished()){
    if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
      const double *nodeCoordsFlat = loop.subtreeInfo().getNodeCoords();
      m_BoundaryOctant = loop.subtreeInfo().isElementBoundary();
      m_level = loop.getCurrentSubtree().getLevel();
      std::memcpy(m_coords,nodeCoordsFlat, sizeof(DENDRITE_REAL)*nPe_*DIM);
      m_octToPhysical.convertCoordsToPhys(m_coords,nPe_);
      constructTalyGrid(m_coords);
      TALYFEMLIB::FEMElm fe(&m_grid, TALYFEMLIB::BASIS_ALL);
      fe.refill(0,this->getRelativeOrder(elemID));
      traverseOperation(fe);
      elemID++;
      loop.next();

    }
    else{
      loop.step();

    }
  }
  return true;
}

bool Traversal::traverseByValues(){
  if(m_octDA->isActive()) {
    VecGetArrayRead(v_.v, &v_.val);
    int lsize;
    VecGetSize(v_.v, &lsize);
    PetscScalar *ghostedArray;
    m_octDA->nodalVecToGhostedNodal(v_.val, ghostedArray, false, v_.ndof);
    m_octDA->readFromGhostBegin(ghostedArray, v_.ndof);
    m_octDA->readFromGhostEnd(ghostedArray, v_.ndof);
    const size_t sz = m_octDA->getTotalNodalSz();
    auto partFront = m_octDA->getTreePartFront();
    auto partBack = m_octDA->getTreePartBack();
    const auto tnCoords = m_octDA->getTNCoords();
    const int &ndof = v_.ndof;
    DENDRITE_UINT elemID = 0;
    ot::MatvecBase<DIM, PetscScalar> treeloop(sz, ndof, eleOrder_, tnCoords, ghostedArray, &(*m_treePart.cbegin()),
                                              m_treePart.size(), *partFront, *partBack);

    while (!treeloop.isFinished()) {
      if (treeloop.isPre() && treeloop.subtreeInfo().isLeaf()) {
        const double *nodeCoordsFlat = treeloop.subtreeInfo().getNodeCoords();

        const PetscScalar *nodeValsFlat = treeloop.subtreeInfo().readNodeValsIn();
        m_level = treeloop.getCurrentSubtree().getLevel();
        m_BoundaryOctant = treeloop.subtreeInfo().isElementBoundary();
        std::memcpy(m_coords, nodeCoordsFlat, sizeof(DENDRITE_REAL) * nPe_ * DIM);
        m_octToPhysical.convertCoordsToPhys(m_coords, nPe_);
        constructTalyGrid(m_coords);
        TALYFEMLIB::FEMElm fe(&m_grid, TALYFEMLIB::BASIS_ALL);
        fe.refill(0, this->getRelativeOrder(elemID));
        traverseOperation(fe, nodeValsFlat);
        elemID++;
        treeloop.next();
      } else
        treeloop.step();
    }


    size_t writtenSz = treeloop.finalize(ghostedArray);

    if (sz > 0 && writtenSz == 0)
      std::cerr << "Warning: matvec() did not write any data! Loop misconfigured?\n";

    VecRestoreArrayRead(v_.v, &v_.val);
    delete [] ghostedArray;
  }

  return true;
}

const DENDRITE_UINT Traversal::getNdof() const{
  assert(m_traversalType == TRAVERSAL_TYPE::VALUES);
  return v_.ndof;
}


bool Traversal::traverse() {

  using namespace TALYFEMLIB;

  if(m_traversalType==TRAVERSAL_TYPE::COORDS){
    return (traverseByCoords());
  }
  else if(m_traversalType == TRAVERSAL_TYPE::VALUES){

    return (traverseByValues());
  }
  return false;
}

void Traversal::init(){
  using namespace TALYFEMLIB;

  nPe_ = m_octDA->getNumNodesPerElement();
  eleOrder_ = m_octDA->getElementOrder();
  m_grid.redimArrays(nPe_,1);
  for(DENDRITE_UINT i = 0; i < nPe_; i++){
      m_grid.node_array_[i] = new TALYFEMLIB::NODE();
  }
#if (DIM == 3)
  m_elem = new  ELEM3dHexahedral();
#elif(DIM == 4)
  m_elem = new  ELEM4dTesseract();
#elif(DIM == 2)
  m_elem = new ELEM2dBox();
#else
  throw TALYFEMLIB::TALYException() <<"Element for nsd = " << DIM  << "Not supported \n";
#endif

  m_grid.elm_array_[0] = m_elem;
  node_id_array_ = new int[nPe_];
  for(int i = 0; i < nPe_; i++){
    node_id_array_[i] = i;
  }
  m_elem->redim(nPe_,node_id_array_);
  m_coords = new DENDRITE_REAL[nPe_*DIM];
}
Traversal::~Traversal() {
  delete [] m_coords;
}