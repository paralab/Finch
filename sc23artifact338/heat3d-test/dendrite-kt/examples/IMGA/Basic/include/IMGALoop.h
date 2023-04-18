//
// Created by maksbh on 4/1/21.
//

#ifndef DENDRITEKT_IMGALOOP_H
#define DENDRITEKT_IMGALOOP_H

#include <IMGA/IMGATraversal.h>
#include <DendriteUtils.h>

struct GPValues{
  DENDRITE_REAL position[DIM];
  DENDRITE_REAL value;
};

class IMGALoop: public IMGATraversal{
  std::vector<GPValues> gpValues_;
  int counter = 0;
public:
  IMGALoop(DA * octDA, const IMGA * imga, const std::vector<TREENODE> &treePart, const VecInfo &v, const DomainExtents &domain);

  void imgaTraversalOperation(const TALYFEMLIB::FEMElm & fe,const NodeAndValues<DENDRITE_REAL> & gaussPoint, const TALYFEMLIB::ZEROPTV & h, const PetscScalar * values) override;

  void computeBoundaryError(DENDRITE_REAL * globalError);
};

IMGALoop::IMGALoop(DA *octDA, const IMGA *imga, const std::vector<TREENODE> &treePart, const VecInfo &v,
                   const DomainExtents &domain)
                   :IMGATraversal(octDA,imga,treePart,v,domain){
  gpValues_.resize(imga->getSurfaceGaussPoints().size());
  this->imgaTraverse();
}

void IMGALoop::imgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint, const TALYFEMLIB::ZEROPTV & h,
                                      const PetscScalar *values) {
  DENDRITE_REAL val_c;
  calcValueFEM(fe,1,values,&val_c);
  for(int dim = 0; dim < DIM; dim++){
    gpValues_[counter].position[dim] = gaussPoint.location[dim];
  }
  gpValues_[counter].value = val_c;
  counter++;
}

void IMGALoop::computeBoundaryError(DENDRITE_REAL * globalError) {
  DENDRITE_REAL error[2]{0.0,-INFINITY};
  std::ofstream fout("boundaryValues_"+ TALYFEMLIB::ToString(TALYFEMLIB::GetMPIRank()));
  for(const auto & gp:gpValues_){
    error[0] += fabs(gp.value);
    if(error[1] < fabs(gp.value)){
      error[1] = fabs(gp.value);
    }

    fout << gp.position[0] << " " << gp.position[1] << " " << gp.value << "\n";
  }
  fout.close();
  MPI_Reduce(&error[0],&globalError[0],1,MPI_DOUBLE,MPI_SUM,0,m_octDA->getCommActive());
  MPI_Reduce(&error[1],&globalError[1],1,MPI_DOUBLE,MPI_MAX,0,m_octDA->getCommActive());

}
#endif //DENDRITEKT_IMGALOOP_H
