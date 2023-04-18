#pragma once

#include <talyfem/fem/cequation.h>
#include "BTNodeData.h"

class BTEquation : public TALYFEMLIB::CEquation<BTNodeData> {
 public:
  void Solve(double delta_t, double current_time) override {
    assert(false);
  }
  void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                  TALYFEMLIB::ZEROARRAY<double> &be) override {
    assert(false);
  }

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae);
  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be);
};

// note: BC on x- and x+ walls at xl = -0.4 and xr = 0.4
// dirichlet 0 (because nonlinear), IC -0.668371 for left and 3.017089 for right

void BTEquation::Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae) {
  using namespace TALYFEMLIB;

  const int spatial_dims = fe.nsd();
  const int nbf = fe.nbf();
  const double detJxW = fe.detJxW();
  const double lambda = -M_PI * M_PI;

  const double u = this->p_data_->valueFEM(fe, 0);

  const double lambda_exp_u = lambda * exp(u);

  double N;
  for (int a = 0; a < nbf; a++) {
    for (int b = 0; b < nbf; b++) {
      N = 0;
      for (int k = 0; k < spatial_dims; k++) {
        N += fe.dN(a, k) * fe.dN(b, k);
      }
      Ae(a, b) += (N - lambda_exp_u * fe.N(a) * fe.N(b)) * detJxW;
    }
  }
}

void BTEquation::Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be) {
  using namespace TALYFEMLIB;

  const int spatial_dims = fe.nsd();
  const int nbf = fe.nbf();
  const double detJxW = fe.detJxW();
  const double lambda = -M_PI * M_PI;

  double u = this->p_data_->valueFEM(fe, 0);

  const double lambda_exp_u = lambda * exp(u);

  ZEROPTV du;
  for (int i = 0; i < spatial_dims; i++) {
    du(i) = this->p_data_->valueDerivativeFEM(fe, 0, i);
  }

  for (int a = 0; a < nbf; a++) {
    double du_dw = 0;
    for (int i = 0; i < spatial_dims; i++) {
      du_dw += du(i) * fe.dN(a, i);
    }
    // note: no forcing term because it cancels out for this particular manufactured solution
    be(a) += -(-du_dw + lambda_exp_u * fe.N(a)) * detJxW;
  }
}
