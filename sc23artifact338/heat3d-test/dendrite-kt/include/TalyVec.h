//
// Created by maksbh on 9/18/19.
//

//
#ifndef DENDRITEKT_TALYVEC_H
#define DENDRITEKT_TALYVEC_H

/**
 * @brief computes the RHS vector.
 */

#include <DataTypes.h>
#include <feVector.h>
#include <TalyDendroSync.h>
#include <Basis/Vec.h>
#include <TimeInfo.h>
#include <MatVecCommon.h>
/**
 *
 * @tparam Equation Equation class
 * @tparam NodeData Node Data class
 * [Convention] If you are looping over the dimensions: please use d for the loop.
 */
template<class Equation, class NodeData>
class TalyVec : public feVector<TalyVec<Equation, NodeData>, DIM> {
  TalyDendroSync sync_;
  TALYFEMLIB::GRID *talyGrid_;
  TALYFEMLIB::GridField<NodeData> *talyGf_;
  TALYFEMLIB::ELEM *talyElem_;

  Equation *talyEq_;
  TALYFEMLIB::ZEROARRAY<DENDRITE_REAL> be_;
  DENDRITE_UINT eleOrder_;
  DENDRITE_UINT nPe_;
  bool surfAssembly_;

  TALYFEMLIB::FEMElm fe_;

  double timer_ = 0;
  std::vector<VecInfo> syncValues_;

  DENDRITE_REAL  * physCoords_;
  TalyMatVecCommon<NodeData> & matVecCommon_;
  const TimeInfo * ti_ = nullptr;
#ifdef IBM
  std::vector<DENDRITE_REAL> surfaceValues_;
#endif
  void setPlaceholder(const Vec &v);
#ifdef TENSOR
  const RefElement *refElement;
  TensorVec::Vec *tensorVec;
#endif
 protected:

  typedef feVector<TalyVec<Equation, NodeData>, DIM> Parent;
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
  TalyVec(ot::DA<DIM> *da, const std::vector<TREENODE> & treePart, Equation *eq, TALYFEMLIB::GRID *grid,
          TALYFEMLIB::GridField<NodeData> *gf, unsigned int dof, bool surfAssembly, TalyMatVecCommon<NodeData> & talyMatVecCommon);

  /**@biref elemental compute vec for rhs*/
  /**
   *
   * @param in Input vector
   * @param out Output assembled local element vector for rhs
   * @param coords coordinates of size [nPe * DIM]
   * @param scale not needed. Basically is the Jaccobian.
   */
  virtual void elementalComputeVec(const VECType *in,
                                   VECType *out,const DENDRITE_UINT ndof,
                                   const double *coords,
                                   double scale , bool isElementalBoundary) override;



  /**
  * @brief Called before computeVec, used to initialize stuff for computeVec computation
  * @param in  vec in (pass null ptr)
  * @param out  vec out (pass null ptr)
  * @param scale
  * @return true
  */
  bool preComputeVec(const VECType *in, VECType *out, double scale = 1.0) ;

  /**
   * @brief Called after computeVec, used to initialize stuff for computeVec computation
   * @param in
   * @param out
   * @param scale
   * @return
   */

  bool postComputeVec(const VECType *in, VECType *out, double scale = 1.0);

  ~TalyVec();
  /**
   * @brief The vectors for synchronizing the node Data like non-linear guess and previous solution
   * @param vecs vectors for synchronizing
   */
  void setVector(const std::vector<VecInfo> & v);

  /**
   * @brief set time construct
   * @param ti time construct
   */
  void setTime(const TimeInfo * ti);

  void performSurfaceVecAssembly();


};

template<class Equation, class NodeData>
TalyVec<Equation, NodeData>::TalyVec(ot::DA<DIM> *da, const std::vector<TREENODE> & treeNode, Equation *eq, TALYFEMLIB::GRID *grid,
                                     TALYFEMLIB::GridField<NodeData> *gf, unsigned int ndof, bool surfAssembly,TalyMatVecCommon<NodeData> & matVecCommon)
    : feVector<TalyVec<Equation, NodeData>, DIM>(da, &treeNode,ndof),
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
  tensorVec = new TensorVec::Vec(eleOrder_,nPe_,refElement);
#endif
  physCoords_ = new DENDRITE_REAL[nPe_*DIM];
}

template <class Equation, class NodeData>
void TalyVec<Equation, NodeData>::setTime(const TimeInfo *ti) {
  ti_ = ti;
}

template<class Equation, class NodeData>
TalyVec<Equation, NodeData>::~TalyVec(){

delete [] physCoords_;

}

template<class Equation, class NodeData>
void TalyVec<Equation, NodeData>::setVector(const std::vector<VecInfo> &v) {
  syncValues_ = v;
}

template <class Equation, class NodeData>
bool TalyVec<Equation,NodeData>::preComputeVec(const VECType *in, VECType *out, double scale) {
  if(ti_ != nullptr){
    talyEq_->set_t(ti_->getCurrentTime());
    talyEq_->set_dt(ti_->getCurrentStep());
  }


  be_.redim(nPe_*m_uiDof);
  matVecCommon_.initMatVecOperation(m_uiOctDA,syncValues_,TalyMatVecCommon<NodeData>::CALL::VEC_ASSEMBLY);
#ifdef IBM
  if((surfaceValues_.empty()) and (matVecCommon_.getDofSurface() > 0)) {
    surfaceValues_.resize(matVecCommon_.getDofSurface());
  }
#endif
  return true;

}

template <class Equation, class NodeData>
void TalyVec<Equation,NodeData>::elementalComputeVec(const VECType *in,
                                                     VECType *out,
                                                     const DENDRITE_UINT ndof,
                                                     const double *coords,
                                                     double scale, bool isElementalBoundary) {

  matVecCommon_.convertToPhys(coords,physCoords_);
#if TENSOR

  timer_ += talyEq_->Integrands_Be(tensorVec, in ,out, coords);
#else
  be_.fill(0.0);

  if(eleOrder_ == 1){
    sync_.syncCoords<1>(physCoords_,talyGrid_);
  }
  else if(eleOrder_ == 2){
    sync_.syncCoords<2>(physCoords_,talyGrid_);
  }
  matVecCommon_.mapToGridField(talyGf_);

#ifdef IBM
    bool doGlobalAssemble = false;
    fe_.refill(0,0);

    if(matVecCommon_.getIBMMethod() == IBM_METHOD::NITSCHE) {
        /// TODO: Adaptive quadrature
        if(matVecCommon_.getElementMarkerType() == ElementMarker::OUT_GP){
            while (fe_.next_itg_pt()) {
                talyEq_->Integrands_be(fe_, be_);
            }
        }
        if(matVecCommon_.getElementMarkerType() == ElementMarker::INTERCEPTED_GP) {
            while (fe_.next_itg_pt()) {
                if (not(matVecCommon_.isInside(fe_.position()))) {
                    talyEq_->Integrands_be(fe_, be_);
                    doGlobalAssemble = true;
                }
            }
        }
        NodeAndValues<DENDRITE_REAL> gaussPoint;
        /// Surface STL/MSH assembly (if any)
        while (true){
            if(matVecCommon_.getLocalSurfaceGaussPoint(fe_,gaussPoint,surfaceValues_.data())){
                doGlobalAssemble = true;
                TALYFEMLIB::ZEROPTV h,position;
                for(int d = 0; d < DIM ; d++){
                    h(d) =  physCoords_[DIM*(nPe_ - 1) + d] - physCoords_[d];
                    position(d) = physCoords_[d];
                }
                talyEq_->ibm_Integrands4side_be(fe_, be_ ,gaussPoint,position,h);
            }
            else{
                break;
            }
        }
    } else{
        while (fe_.next_itg_pt()) {
            talyEq_->Integrands_be(fe_, be_);
        }
    }
    if(surfAssembly_ and isElementalBoundary){
        performSurfaceVecAssembly();
    }

#else

  fe_.refill(0,matVecCommon_.getRelativeOrder());
#ifdef PROFILING
  PetscLogEventBegin(Profiling::vecGaussPointAssembly,0,0,0,0);
#endif

  while (fe_.next_itg_pt()) {
    talyEq_->Integrands_be(fe_, be_);
  }
#ifdef PROFILING
  PetscLogEventEnd(Profiling::vecGaussPointAssembly,0,0,0,0);
#endif
  if(isElementalBoundary and surfAssembly_){
    performSurfaceVecAssembly();
  }
#endif
  for(int i = 0; i < nPe_*m_uiDof; i++){
    out[i] = be_(i);
  }
  matVecCommon_.finishElementOperation();
#endif

}

template <class Equation, class NodeData>
bool TalyVec<Equation,NodeData>::postComputeVec(const VECType *in, VECType *out, double scale) {
  matVecCommon_.finalizeMatVecoperation();
  return true;
}

template<class Equation, class NodeData>
void TalyVec<Equation, NodeData>::setPlaceholder(const Vec &v) {
  syncValues_[0].v = v;
//  matVecCommon.setPlaceHolder(v);
}

template<class Eqaution, class NodeData>
void TalyVec<Eqaution, NodeData>::performSurfaceVecAssembly() {


  std::array<BoundarySurface, 2 * DIM> boundarySurfaces;
  matVecCommon_.getSurfaceID(physCoords_, boundarySurfaces);
  const int *surf_arr = talyElem_->GetSurfaceCheckArray();
  const int surf_row_len = talyElem_->GetSurfaceCheckArrayRowLength();


  for (int i = 0; i < 2*DIM; i++) {
    const DENDRITE_UINT & boundaryType = boundarySurfaces[i].boundaryType;
    int surf_id = surf_arr[i * surf_row_len];  // id determines which nodes on the surface
    if (boundaryType != -1) {
      TALYFEMLIB::SurfaceIndicator surf(surf_id);
      TALYFEMLIB::ZEROPTV  normal = talyElem_->CalculateNormal(talyGrid_, surf_id);

      /*
      if(boundaryType >= BoundaryTypes::WALL::MAX_WALL_TYPE_BOUNDARY) {
        surf.set_normal(normal*-1);
      }
      else
      {
        surf.set_normal(normal);
      }
      */

      surf.set_normal(normal);
      
      fe_.refill_surface(talyElem_, &surf, 0);
      // loop over surface gauss points
      while (fe_.next_itg_pt()) {
        talyEq_->Integrands4side_be(fe_, boundaryType,boundarySurfaces[i].id, be_);
      }
    }
  }
}
#endif //DENDRITEKT_TALYVEC_H
