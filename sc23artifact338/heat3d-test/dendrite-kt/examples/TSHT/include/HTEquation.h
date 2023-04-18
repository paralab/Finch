//
// Created by maksbh on 5/20/20.
//

#ifndef DENDRITEKT_TSHTEQUATION_H
#define DENDRITEKT_TSHTEQUATION_H
#include <talyfem/fem/cequation.h>
#include "HTNodeData.h"
#include "HTInputData.h"
using namespace TALYFEMLIB;
class HTEquation : public TALYFEMLIB::CEquation<HTNodeData> {
 public:

  explicit HTEquation(const HTInputData * idata)
      : TALYFEMLIB::CEquation<HTNodeData>(false, TALYFEMLIB::kAssembleGaussPoints) {
    K_ = 1.0 / (3 * M_PI * M_PI);
    idata_ = idata;
    if (idata_->forcing) {
      K_ = 1.0;
    }
  }
  // never called, but required for HTEquation to be instantiable (pure virtual in CEquation)
  virtual void Solve(double dt, double t) { assert(false); }

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae) {
    using namespace TALYFEMLIB;
    // # of dimensions: 1D, 2D, or 3D
    const int n_dimensions = fe.nsd();
    // # of basis functions
    const int n_basis_functions = fe.nbf();
    // (determinant of J) cross W
    const double detJxW = fe.detJxW();
    // thermal diffusivity in heat equation
    double k_val = K_;

    // in order to assemble the gauss point, we loop over each pair of basis
    // functions and calculate the individual contributions to the 'Ae' matrix
    // and 'be' vector.
    for (int a = 0; a < n_basis_functions; a++) {
      for (int b = 0; b < n_basis_functions; b++) {
        double M = 0.0;
        if (idata_ && idata_->BDF2 and t_ > 0.99 * dt_) {
          M = 3.0 / 2.0 * fe.N(a) * fe.N(b) * detJxW;
        } else {
          M = fe.N(a) * fe.N(b) * detJxW;
        }
        double N = 0;
        for (int k = 0; k < n_dimensions; k++) {
          N += k_val * fe.dN(a, k) * fe.dN(b, k) * detJxW;
        }
        // Add term to the A element matrix
        Ae(a, b) += M / dt_ + N;
      }
    }
  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be) {
    using namespace TALYFEMLIB;
    // # of basis functions
    const int n_basis_functions = fe.nbf();
    // (determinant of J) cross W
    const double detJxW = fe.detJxW();

    const double u_pre_curr = p_data_->valueFEM(fe, U_PRE);
    const double u_pre_pre_curr = p_data_->valueFEM(fe, U_PRE_PRE);
    double forcing = 0.0;
    if (idata_ && idata_->forcing) {
      forcing = this->forcing(fe, t_ + dt_);
    }
    for (int a = 0; a < n_basis_functions; a++) {
      if (idata_ && idata_->BDF2 and t_ > 0.99 * dt_) {
        be(a) += fe.N(a) / 2.0 / dt_ * (4 * u_pre_curr - u_pre_pre_curr) * detJxW + fe.N(a) * forcing * detJxW;
      } else {
        be(a) += fe.N(a) / dt_ * u_pre_curr * detJxW + fe.N(a) * forcing * detJxW;
      }
    }
  }

 private:
  double K_;
  const HTInputData *idata_;

  static double forcing(const FEMElm &fe, double t) {
    ZEROPTV p = fe.position();
    double forcing = 0.0;
    if (DIM == 2) {
      forcing = (-2 * M_PI * sin(2 * M_PI * t) * sin(M_PI * p.x()) * sin(M_PI * p.y())) +
        (DIM * M_PI * M_PI * sin(M_PI * p.x()) * sin(M_PI * p.y()) * cos(2 * M_PI * t));
    } else if (DIM == 3) {
     forcing = (-2 * M_PI * sin(2 * M_PI * t) * sin(M_PI * p.x()) * sin(M_PI * p.y()) * sin(M_PI * p.z())) +
        (DIM * M_PI * M_PI * sin(M_PI * p.x()) * sin(M_PI * p.y()) * sin(M_PI * p.z()) * cos(2 * M_PI * t));
    }
    return forcing;
  }
};
#endif //DENDRITEKT_TSHTEQUATION_H
