//
// Created by maksbh on 5/21/20.
//

#ifndef DENDRITEKT_HTANALYTIC_H
#define DENDRITEKT_HTANALYTIC_H
#include <PostProcessing/postProcessing.h>
#include <DendriteUtils.h>
#include <TimeInfo.h>
class HTAnalytic : public postProcessor {
  DENDRITE_REAL err = 0;
  const TimeInfo * ti_;
 public:
  HTAnalytic(ot::DA<DIM> *da, const Vec & inVec, int ndof,const  TimeInfo* ti);

  DENDRITE_REAL getAnalyticValues(const TALYFEMLIB::ZEROPTV &pos);

  void calcFunction(TALYFEMLIB::FEMElm &fe, DENDRITE_REAL *val) override;

  DENDRITE_REAL getL2Error();
};

HTAnalytic::HTAnalytic(ot::DA<DIM> *da, const Vec & inVec, int ndof, const  TimeInfo* ti)
    : postProcessor(da, inVec, ndof) {
  ti_ = ti;
  this->loop();
}

void HTAnalytic::calcFunction(TALYFEMLIB::FEMElm &fe, DENDRITE_REAL *val) {
  double val_c;
  while (fe.next_itg_pt()) {
    calcValueFEM(fe, 1, val, &val_c);
    DENDRITE_REAL val_a = getAnalyticValues(fe.position());
    err += (val_c - val_a) * (val_c - val_a) *fe.detJxW();
  }
}

DENDRITE_REAL HTAnalytic::getL2Error() {
  DENDRITE_REAL l2_error;
  MPI_Reduce(&err, &l2_error, 1, MPI_DOUBLE, MPI_SUM, 0, this->m_octDA->getCommActive());
  return (sqrt(l2_error));
}

DENDRITE_REAL HTAnalytic::getAnalyticValues(const TALYFEMLIB::ZEROPTV &pos) {
  return (sin(M_PI * pos.x()) * sin(M_PI * pos.y()) * sin(M_PI * pos.z()) * exp(-ti_->getCurrentTime()));
}

#endif //DENDRITEKT_HTANALYTIC_H
