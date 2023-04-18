//
// Created by maksbh on 4/1/21.
//

#ifndef DENDRITEKT_IMGATRAVERSAL_H
#define DENDRITEKT_IMGATRAVERSAL_H

#include <Traversal/Traversal.h>
#include <IMGA/IMGA.h>
enum IMGATraversalType:bool{
  POSTPROCESSING = false,
  INTERPOLATE = true
};

class IMGATraversal: public Traversal{
  bool hasLocalGaussPoint = true;
  std::vector<NodeAndValues<DENDRITE_REAL>>::const_iterator gpInfoIt_;
  std::vector<NodeAndValues<DENDRITE_REAL>>::const_iterator gpInfoEnd_;

  DENDRITE_UINT lclElemID = 0;
  const IMGATraversalType  traversalType_;

public:
  IMGATraversal(DA * octDA, const IMGA * imga, const std::vector<TREENODE> &treePart, const VecInfo &v, const DomainExtents &domain, const IMGATraversalType _traversalType = IMGATraversalType::POSTPROCESSING);

  void imgaTraverse();

  virtual void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;

  virtual void imgaTraversalOperation(const TALYFEMLIB::FEMElm & fe,const NodeAndValues<DENDRITE_REAL> & gaussPoint,const TALYFEMLIB::ZEROPTV & h,const PetscScalar * values) = 0;
};

IMGATraversal::IMGATraversal(DA * octDA, const IMGA * imga, const std::vector<TREENODE> &treePart, const VecInfo &v, const DomainExtents &domain,const IMGATraversalType _traversalType)
:Traversal(octDA,treePart,v,domain),traversalType_(_traversalType){
  const auto & gpInfo_ = imga->getSurfaceGaussPoints();
  if(gpInfo_.empty()){
    hasLocalGaussPoint = false;
  }
  else {
    gpInfoIt_ = gpInfo_.cbegin();
    gpInfoEnd_ = gpInfo_.cend();
  }
  MPI_Barrier(octDA->getCommActive());

}

void IMGATraversal::imgaTraverse() {
    this->traverse();
}

void IMGATraversal::traverseOperation(TALYFEMLIB::FEMElm & fe,const PetscScalar * values){
  if(hasLocalGaussPoint) {
    while ((gpInfoIt_ < gpInfoEnd_) and (gpInfoIt_->localElemID==lclElemID)) {
      const auto gaussPoint = *gpInfoIt_;

      /// TODO Fix it. Seems stupid;
      TALYFEMLIB::ZEROPTV ptvg, ptvl;
      std::memcpy(ptvg.data(), gaussPoint.location, sizeof(double)*DIM);
      GetLocalPtv(fe, ptvg, ptvl);
      fe.calc_at(ptvl);
      const auto &coords = this->m_coords;
      TALYFEMLIB::ZEROPTV h;
      const DENDRITE_UINT nPe = this->m_octDA->getNumNodesPerElement();
      for (int d = 0; d < DIM; d++) {
        h.data()[d] = coords[(nPe - 1)*DIM + d] - coords[d];
      }
      if (traversalType_==IMGATraversalType::POSTPROCESSING) {
        static constexpr DENDRITE_REAL weight = (DIM==2) ? 2.0 : 3.0;
        fe.set_jacc_x_w(gaussPoint.elemArea/weight);
      }
      imgaTraversalOperation(fe, *gpInfoIt_, h, values);
      gpInfoIt_ = std::next(gpInfoIt_);

    }
  }
  lclElemID++;
}


#endif //DENDRITEKT_IMGATRAVERSAL_H
