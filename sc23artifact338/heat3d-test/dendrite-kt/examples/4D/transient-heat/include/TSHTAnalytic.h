//
// Created by maksbh on 9/19/19.
//

#ifndef DENDRITEKT_ANALYTICSOLUTION_H
#define DENDRITEKT_ANALYTICSOLUTION_H

#include <PostProcessing/postProcessing.h>
#include <DendriteUtils.h>

class TSHTAnalytic : public postProcessor {

  DENDRITE_REAL err = 0;
 public:
  TSHTAnalytic(ot::DA<DIM> *da, const Vec & inVec, int ndof);

  DENDRITE_REAL getAnalyticValues(const TALYFEMLIB::ZEROPTV &pos);

  void calcFunction(TALYFEMLIB::FEMElm &fe, DENDRITE_REAL *val) override;

  DENDRITE_REAL getL2Error();
};

TSHTAnalytic::TSHTAnalytic(ot::DA<DIM> *da, const Vec & inVec, int ndof)
    : postProcessor(da, inVec, ndof) {

  this->loop();
}

void TSHTAnalytic::calcFunction(TALYFEMLIB::FEMElm &fe, DENDRITE_REAL *val) {
  double val_c;
  while (fe.next_itg_pt()) {
    calcValueFEM(fe, 1, val, &val_c);
    DENDRITE_REAL val_a = getAnalyticValues(fe.position());
   // std::cout << val_c << " " << val_a << " " << fe.position() << "\n";
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
  if(DIM == 3) {
    return (sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z));
  }
  else if(DIM == 4) {
    return (sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) * exp(- pos.t()));
  }

}
#endif //DENDRITEKT_ANALYTICSOLUTION_H
