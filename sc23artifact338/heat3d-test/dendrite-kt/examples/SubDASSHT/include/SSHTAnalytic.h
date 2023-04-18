//
// Created by maksbh on 9/19/19.
//

#ifndef DENDRITEKT_ANALYTICSOLUTION_H
#define DENDRITEKT_ANALYTICSOLUTION_H

#include <PostProcessing/postProcessing.h>
#include <DendriteUtils.h>
#include <Traversal/Traversal.h>

class TSHTAnalytic : public Traversal {

  DENDRITE_REAL err = 0;
 public:
  TSHTAnalytic(ot::DA<DIM> *da, const VecInfo & inVec);

  DENDRITE_REAL getAnalyticValues(const TALYFEMLIB::ZEROPTV &pos);

  void traverseOperation(TALYFEMLIB::FEMElm & fe,const PetscScalar * values) override;

  DENDRITE_REAL getL2Error();
};

TSHTAnalytic::TSHTAnalytic(ot::DA<DIM> *da, const VecInfo & inVec)
    : Traversal(da, inVec) {

  this->traverse();
}

void TSHTAnalytic::traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values) {
  double val_c = 0;
  while (fe.next_itg_pt()) {
    calcValueFEM(fe, 1, values, &val_c);
    DENDRITE_REAL val_a = getAnalyticValues(fe.position());
    err += (val_c - val_a) * (val_c - val_a) *fe.detJxW();
  }
}

DENDRITE_REAL TSHTAnalytic::getL2Error() {
  DENDRITE_REAL l2_error;
  MPI_Reduce(&err, &l2_error, 1, MPI_DOUBLE, MPI_SUM, 0, this->m_octDA->getCommActive());
  return (sqrt(l2_error));
}

DENDRITE_REAL TSHTAnalytic::getAnalyticValues(const TALYFEMLIB::ZEROPTV &pos) {
  DENDRITE_REAL x(pos.x()), y(pos.y()), z(pos.z());
  return (sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z));

}
#endif //DENDRITEKT_ANALYTICSOLUTION_H
