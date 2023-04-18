//
// Created by maksbh on 9/18/19.
//

#ifndef DENDRITEKT_TALYMAT_H
#define DENDRITEKT_TALYMAT_H
#include <DataTypes.h>
#include <feMatrix.h>
#include <TalyDendroSync.h>
#include <Basis/Basis.h>
#include <Basis/MatVec.h>
#include <Basis/Mat.h>
#include <PETSc/VecInfo.h>
#include <MatVecCommon.h>
#include "TimeInfo.h"
#include "OctToPhysical.h"
#ifdef PROFILING
#include <Profiling/GlobalProfiler.h>
#endif
/**
 *
 * @tparam Equation Equation class
 * @tparam NodeData Node Data class
 * [Convention] If you are looping over the dimensions: please use d for the loop.
 */
template<class Equation, class NodeData>
class TalyMat : public feMatrix<TalyMat<Equation, NodeData>, DIM> {

  TalyDendroSync sync_;
  TALYFEMLIB::GRID *talyGrid_;
  TALYFEMLIB::GridField<NodeData> *talyGf_;
  TALYFEMLIB::ELEM *talyElem_;

  Equation *talyEq_;
  TALYFEMLIB::ZeroMatrix<DENDRITE_REAL> Ae_;
  DENDRITE_UINT eleOrder_;
  DENDRITE_UINT nPe_;
  bool surfAssembly_;

  double timerTalyMat = 0;
  double timerTalyMatVec = 0;
  double timerMatVec= 0;
  double timerMat= 0;
  std::vector<VecInfo> syncValues_;
  TalyMatVecCommon<NodeData> & matVecCommon_;
  const TimeInfo * ti_ = nullptr;
  DENDRITE_REAL * physCoords_ = nullptr;
  TALYFEMLIB::FEMElm fe_;
#ifdef TENSOR
  const RefElement *refElement;
  TensorMatVec::MatVec *tensorMatVec;
  TensorMat::Mat *tensorMat;
  DENDRITE_REAL *mat;
  DENDRITE_REAL * outTemp;
#endif
#ifdef IBM
  std::vector<DENDRITE_REAL> surfaceValues_;
#endif
protected:

  typedef feMatrix<TalyMat<Equation, NodeData>, DIM> Parent;
  using Parent::m_uiOctDA;
  using Parent::m_uiDof;
  using Parent::m_uiPtMin;
  using Parent::m_uiPtMax;

 public:
  using Parent::setProblemDimensions;

/**
 *
 * @tparam Equation Equation class
 * @tparam NodeData NodeData class
 * @param da        octant
 * @param eq        pointer to equation class
 * @param grid      pointer to grid object
 * @param gf        pointer to gridfield
 * @param ndof      ndof corresponding to equation class
 * @param surfAssembly  true for Integrands4Side Robins/Neumann boundary condition.
 */
  TalyMat(ot::DA<DIM> *da, const std::vector<TREENODE> & treePart, Equation *eq, TALYFEMLIB::GRID *grid,
          TALYFEMLIB::GridField<NodeData> *gf, unsigned int dof, bool surfAssembly, TalyMatVecCommon<NodeData> & matVecCommon);

  /**@biref elemental compute vec for rhs*/
  /**
   *
   * @param in Input vector
   * @param out Output assembled Matrix \times in :local matVec
   * @param coords coordinates of size [nPe * DIM]
   * @param scale not needed. Basically is the Jacobian.
   */
  virtual void elementalMatVec(const VECType *in, VECType *out,const DENDRITE_UINT ndof, const double *coords, double scale, bool isElementalBoundary) override ;

  /**
   * @brief Called before MatVec, used to initialize stuff for matVec computation
   * @param in  vec in (pass null ptr)
   * @param out  vec out (pass null ptr)
   * @param scale
   * @return true
   */
  bool preMatVec(const VECType *in, VECType *out, double scale = 1.0);

  /**
   * @brief Called after MatVec, used to finalize stuff after matVec computation
   * @param in  vec in (pass null ptr)
   * @param out  vec out (pass null ptr)
   * @param scale
   * @return true
   */
  bool postMatVec(const VECType *in, VECType *out, double scale = 1.0);

  /**
   * @brief Called before Matrix computation, used to initialize stuff for matrix computation
   * @return
   */
  bool preMat();

  /**
   * @brief Called after matrix computation, used to finalize stuff after matrrix computation
   * @return
   */
  bool postMat();

 /**
  * @brief compute local elemental matrix
  * @param records records for matrix computation
  * @param coords  coordinates
  */
  void getElementalMatrix(std::vector<ot::MatRecord> &records, const double *coords, bool isElementalBoundary) ;

  ~TalyMat();

  /**
   * @brief The vectors for synchronizing the node Data like non-linear guess and previous solution
   * @param vecs vectors for synchronizing
   */
  void setVector(const std::vector<VecInfo> &vecs);

  /**
   * @brief vector to set placeholder for non-linear class
   * @param v
   */
  void setPlaceholder(const Vec &v) override;

  /**
   * @brief set time construct
   * @param ti time construct
   */
  void setTime(const TimeInfo * ti);

  void performSurfaceMatAssembly();
};

template<class Equation, class NodeData>
TalyMat<Equation, NodeData>::TalyMat(ot::DA<DIM> *da, const std::vector<TREENODE> & treePart,Equation *eq, TALYFEMLIB::GRID *grid,
                                     TALYFEMLIB::GridField<NodeData> *gf, unsigned int ndof, bool surfAssembly, TalyMatVecCommon<NodeData> & matVecCommon)
    : feMatrix<TalyMat<Equation, NodeData>, DIM>(da, &treePart,ndof),
      talyEq_(eq),
      talyGf_(gf),
      talyGrid_(grid),
      fe_(talyGrid_,TALYFEMLIB::BASIS_ALL),
      talyElem_(grid->elm_array_[0]),
      surfAssembly_(surfAssembly),
      eleOrder_(da->getElementOrder()),
      nPe_(da->getNumNodesPerElement()),
      matVecCommon_(matVecCommon){
#ifdef TENSOR
        refElement = da->getReferenceElement();
        tensorMatVec = new TensorMatVec::MatVec(eleOrder_,nPe_,refElement);
        tensorMat = new TensorMat::Mat(eleOrder_,nPe_,refElement);
        mat = new double[nPe_*nPe_];
        outTemp = new DENDRITE_REAL[nPe_*nPe_];
#endif
        physCoords_ = new DENDRITE_REAL[nPe_*DIM];

}
template<class Equation, class NodeData>
void TalyMat<Equation, NodeData>::setVector(const std::vector<VecInfo> &vecs) {
    syncValues_ = vecs;

}

template <class Equation, class NodeData>
void TalyMat<Equation, NodeData>::setTime(const TimeInfo *ti) {
  ti_ = ti;
}
template<class Equation, class NodeData>
TalyMat<Equation, NodeData>::~TalyMat(){
  delete [] physCoords_;

}

template<class Equation, class NodeData>
bool TalyMat<Equation,NodeData>::preMatVec(const VECType* in,VECType* out,double scale)
{

  if(ti_ != nullptr){
    talyEq_->set_t(ti_->getCurrentTime());
    talyEq_->set_dt(ti_->getCurrentStep());
  }
  matVecCommon_.initMatVecOperation(m_uiOctDA,syncValues_,TalyMatVecCommon<NodeData>::CALL::MATVEC_ASSEMBLY);
  Ae_.redim(nPe_*m_uiDof,nPe_*m_uiDof);
#ifdef IBM
  if((surfaceValues_.empty()) and (matVecCommon_.getDofSurface() > 0)) {
    surfaceValues_.resize(matVecCommon_.getDofSurface());
  }
#endif
  return true;


}

template<class Equation, class NodeData>
void TalyMat<Equation, NodeData>::elementalMatVec(const VECType *in, VECType *out, const DENDRITE_UINT ndof, const double *coords,
                                                  double scale, bool isElementalBoundary) {
#ifdef IBM
  throw std::runtime_error ("Matrix - free not supported with IBM");
#endif
  memset(out, 0, sizeof(DENDRITE_REAL) * nPe_ * ndof);
  matVecCommon_.convertToPhys(coords,physCoords_);
#if TENSOR
  memset(outTemp, 0, sizeof(DENDRITE_REAL) * nPe_ * ndof);
  memset(mat, 0, sizeof(DENDRITE_REAL) * nPe_ * nPe_ * ndof);
  timerMat += talyEq_->Integrands_Ae(tensorMat, mat, coords);
#else

  Ae_.fill(0.0);
  if(eleOrder_ == 1){
    sync_.syncCoords<1>(physCoords_,talyGrid_);
  }
  else if(eleOrder_ == 2){
    sync_.syncCoords<2>(physCoords_,talyGrid_);
  }
  matVecCommon_.mapToGridField(talyGf_);
  fe_.refill(0,matVecCommon_.getRelativeOrder());
  memset(out,0, sizeof(DENDRITE_REAL)*nPe_*m_uiDof);

  while (fe_.next_itg_pt()) {
    talyEq_->Integrands_Ae(fe_,Ae_);
  }

  if(isElementalBoundary and surfAssembly_){
    performSurfaceMatAssembly();
  }
  matVecCommon_.finishElementOperation();
  for (int i = 0; i < nPe_ * m_uiDof; i++) {
    for (int j = 0; j < nPe_ * m_uiDof; j++) {
      out[i] += in[j] * Ae_(i, j);
    }
  }

#endif
}
template<class Equation, class NodeData>
bool TalyMat<Equation, NodeData>::postMatVec(const VECType *in, VECType *out, double scale){
  matVecCommon_.finalizeMatVecoperation();
  return true;
}

template<class Equation, class NodeData>
bool TalyMat<Equation, NodeData>::preMat(){
  if(ti_ != nullptr){
    talyEq_->set_t(ti_->getCurrentTime());
    talyEq_->set_dt(ti_->getCurrentStep());
  }
  Ae_.redim(nPe_*m_uiDof,nPe_*m_uiDof);
  matVecCommon_.initMatVecOperation(m_uiOctDA,syncValues_,TalyMatVecCommon<NodeData>::CALL::MAT_ASSEMBLY);
#ifdef IBM
  if((surfaceValues_.empty()) and (matVecCommon_.getDofSurface() > 0)) {
    surfaceValues_.resize(matVecCommon_.getDofSurface());
  }
#endif
  return true;
}

template<class Equation, class NodeData>
bool TalyMat<Equation, NodeData>::postMat(){
  matVecCommon_.finalizeMatVecoperation();
  return true;
}

template<class Equation, class NodeData>
void TalyMat<Equation, NodeData>::getElementalMatrix(std::vector<ot::MatRecord> &records, const double *coords, bool isElementalBoundary){
  Ae_.fill(0.0);

  matVecCommon_.convertToPhys(coords,physCoords_);
  if(eleOrder_ == 1){
    sync_.syncCoords<1>(physCoords_,talyGrid_);
  }
  else if(eleOrder_ == 2){
    sync_.syncCoords<2>(physCoords_,talyGrid_);
  }
  matVecCommon_.mapToGridField(talyGf_);
#ifdef IBM
  /// TODO: Adaptive quadrature
  /// Volume Assembly
  fe_.refill(0,0);

  if(matVecCommon_.getIBMMethod() == IBM_METHOD::NITSCHE) {
    /// Out element
    if (matVecCommon_.getElementMarkerType() == ElementMarker::OUT_GP) {
      while (fe_.next_itg_pt()) {
        talyEq_->Integrands_Ae(fe_, Ae_);
      }
    }

    /// Intercepted elements
    if (matVecCommon_.getElementMarkerType() == ElementMarker::INTERCEPTED_GP) {
      while (fe_.next_itg_pt()) {
        if (not(matVecCommon_.isInside(fe_.position()))) {
          talyEq_->Integrands_Ae(fe_, Ae_);
        }
      }
    }


    /// Surface STL/MSH assembly
    NodeAndValues<DENDRITE_REAL> gaussPoint;
    while (true) {
      if (matVecCommon_.getLocalSurfaceGaussPoint(fe_, gaussPoint,surfaceValues_.data())) {
        TALYFEMLIB::ZEROPTV h, position;
        for (int d = 0; d < DIM; d++) {
          h(d) = physCoords_[DIM * (nPe_ - 1) + d] - physCoords_[d];
          position(d) = physCoords_[d];
        }
        talyEq_->ibm_Integrands4side_Ae(fe_, Ae_, gaussPoint, position, h,surfaceValues_);
      } else {
        break;
      }
    }
  }
  else {
//    if (matVecCommon_.getElementMarkerType() != ElementMarker::IN_GP) {
      while (fe_.next_itg_pt()) {
        talyEq_->Integrands_Ae(fe_, Ae_);
      }
//    }
  }
  if(surfAssembly_ and isElementalBoundary){
    performSurfaceMatAssembly();
  }


#else
  fe_.refill(0,matVecCommon_.getRelativeOrder());
#ifdef PROFILING
  PetscLogEventBegin(Profiling::matGaussPointAssembly,0,0,0,0);
#endif

  while (fe_.next_itg_pt()) {
   talyEq_->Integrands_Ae(fe_,Ae_);
  }
#ifdef PROFILING
  PetscLogEventEnd(Profiling::matGaussPointAssembly,0,0,0,0);
#endif
  if(surfAssembly_ and isElementalBoundary){
    performSurfaceMatAssembly();
  }
#endif
  ot::MatRecord matRecord;
  for(int dofi = 0;dofi < m_uiDof; dofi++) {
    for(int dofj = 0;dofj < m_uiDof; dofj++) {
      for (int i = 0; i < nPe_; i++) {
        for (int j = 0; j < nPe_; j++) {
          matRecord.setColID(j);
          matRecord.setRowID(i);
          matRecord.setRowDim(dofi);
          matRecord.setColDim(dofj);
          matRecord.setMatValue(Ae_(i*m_uiDof + dofi, j*m_uiDof + dofj));
          records.push_back(matRecord);
        }
      }
    }
  }
  matVecCommon_.finishElementOperation();
}

template<class Equation, class NodeData>
void TalyMat<Equation, NodeData>::setPlaceholder(const Vec &v) {
  syncValues_[0].v = v;
}

template<class Eqaution, class NodeData>
void TalyMat<Eqaution, NodeData>::performSurfaceMatAssembly() {


  std::array<BoundarySurface,2*DIM> boundarySurfaces;
  matVecCommon_.getSurfaceID(physCoords_,boundarySurfaces);
  const int *surf_arr = talyElem_->GetSurfaceCheckArray();
  const int surf_row_len = talyElem_->GetSurfaceCheckArrayRowLength();

  for (int i = 0; i < 2*DIM; i++) {
    const DENDRITE_UINT & boundaryType = boundarySurfaces[i].boundaryType;
    int surf_id = surf_arr[i * surf_row_len];  // id determines which nodes on the surface
    if (boundaryType != -1) {
      TALYFEMLIB::SurfaceIndicator surf(surf_id);
      TALYFEMLIB::ZEROPTV  normal = talyElem_->CalculateNormal(talyGrid_, surf_id);
      surf.set_normal(normal);

      fe_.refill_surface(talyElem_, &surf, 0);
      // loop over surface gauss points
      while (fe_.next_itg_pt()) {
        talyEq_->Integrands4side_Ae(fe_, boundaryType,boundarySurfaces[i].id, Ae_);
      }
    }
  }
}
#endif //DENDRITEKT_TALYMAT_H
