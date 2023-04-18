#pragma once

#include <talyfem/fem/cequation.h>
#include "SSHTNodeData.h"

class SSHTEquation : public TALYFEMLIB::CEquation<SSHTNodeData> {
 public:
  void Solve(double dt, double t) override {
    assert(false);
  }

  void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                  TALYFEMLIB::ZEROARRAY<double> &be) override {
    assert(false);
  }

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae) {
    using namespace TALYFEMLIB;
    // # of dimensions: 1D, 2D, or 3D
    const int n_dimensions = fe.nsd();
    // # of basis functions
    const int n_basis_functions = fe.nbf();
    // (determinant of J) cross W
    const double detJxW = fe.detJxW();

    for (int a = 0; a < n_basis_functions; a++) {
      for (int b = 0; b < n_basis_functions; b++) {
        double N = 0;
        for (int k = 0; k < n_dimensions; k++) {

          N += fe.dN(a, k) * fe.dN(b, k) * detJxW;
        }
        Ae(a, b) -= N;
      }
    }

  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be) {
    using namespace TALYFEMLIB;
    // # of basis functions
    const int n_basis_functions = fe.nbf();
    // (determinant of J) cross W
    const double detJxW = fe.detJxW();

    const ZEROPTV p = fe.position();
    double force = calc_d2u_at(p);

    for (int a = 0; a < n_basis_functions; a++) {
      be(a) += fe.N(a) * force * detJxW;
    }

  }

 protected:
  double calc_d2u_at(const TALYFEMLIB::ZEROPTV &pt) const {
    return -3 * M_PI * M_PI * sin(M_PI * pt.x()) * sin(M_PI * pt.y()) * sin(M_PI * pt.z());
  }
};
