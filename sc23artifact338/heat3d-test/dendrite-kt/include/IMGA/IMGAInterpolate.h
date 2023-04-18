//
// Created by maksbh on 7/30/21.
//

#ifndef DENDRITEKT_INTERPOLATE_H
#define DENDRITEKT_INTERPOLATE_H

#include <IMGA/IMGATraversal.h>
#include <IMGA/SurfaceValues.h>
#include <DendriteUtils.h>
template<int dof>
class IMGAInterpolate : public IMGATraversal {
  std::vector<SurfaceValues<dof>> &surfaceValues_;
  int counter = 0;
public:

  /**
   *
   * @param octDA octDA
   * @param treePart treePartition nodes
   * @param imga imga construct
   * @param _surfaceValues [In/Out] interpolated surface values
   * @param v vector
   * @param domain domainExtents
   */
  IMGAInterpolate(DA *octDA, const std::vector<TREENODE> &treePart, const IMGA *imga,
                  std::vector<SurfaceValues<dof>> &_surfaceValues, const VecInfo &v, const DomainExtents &domain);

  void imgaTraversalOperation(const TALYFEMLIB::FEMElm & fe,const NodeAndValues<DENDRITE_REAL> & gaussPoint,const TALYFEMLIB::ZEROPTV & h,const PetscScalar * values) override;

};

template<int dof>
IMGAInterpolate<dof>::IMGAInterpolate(DA *octDA, const std::vector<TREENODE> &treePart, const IMGA *imga,
                                      std::vector<SurfaceValues<dof>> &_surfaceValues, const VecInfo &v,
                                      const DomainExtents &domain)
  :IMGATraversal(octDA, imga, treePart, v, domain, IMGATraversalType::INTERPOLATE), surfaceValues_(_surfaceValues) {
    assert(v.ndof == dof);
    const std::size_t  numLocalGaussPoints = imga->getSurfaceGaussPoints().size();
    surfaceValues_.resize(numLocalGaussPoints);
    this->imgaTraverse();
}

template<int dof>
void IMGAInterpolate<dof>::imgaTraversalOperation(const TALYFEMLIB::FEMElm & fe,const NodeAndValues<DENDRITE_REAL> & gaussPoint,const TALYFEMLIB::ZEROPTV & h,const PetscScalar * values){
  calcValueFEM(fe,dof,values,surfaceValues_[counter].values);
  counter++;
}



#endif //DENDRITEKT_INTERPOLATE_H
