//
// Created by maksbh on 7/28/20.
//

#ifndef DENDRITEKT_SSHTEQUATION_H
#define DENDRITEKT_SSHTEQUATION_H
#include <talyfem/fem/cequation.h>
#include "SSHTNodeData.h"
#include <Boundary/DomainBoundary.h>
class SSHTEquation : public TALYFEMLIB::CEquation<SSHTNodeData> {
 public:
  void Solve(double dt, double t) override final{
    assert(false);
  }

  void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                  TALYFEMLIB::ZEROARRAY<double> &be) override final {
    assert(false);
  }

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae) {
    using namespace TALYFEMLIB;
    // # of dimensions: 1D, 2D, or 3D
    const int n_dimensions = DIM;
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

  void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT side_idx, const DENDRITE_UINT id, TALYFEMLIB::ZeroMatrix<double> &Ae) override {

  }

  void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT side_idx, const DENDRITE_UINT id, TALYFEMLIB::ZEROARRAY<double> &be) override {
    // neumann conditions on X-, X+, Z-, and Z+ walls
    using namespace TALYFEMLIB;

    if (side_idx == BoundaryTypes::WALL::X_MINUS || side_idx == BoundaryTypes::WALL::X_PLUS ||
    side_idx == BoundaryTypes::WALL::Y_MINUS || side_idx == BoundaryTypes::WALL::Y_PLUS ) {
      const ZEROPTV &p = fe.position();
      const ZEROPTV &normal = fe.surface()->normal();

      double flux = calc_grad_u_at(p).innerProduct(normal);

      for (int i = 0; i < fe.nbf(); i++) {
        be(i) -= fe.N(i) * fe.detJxW() * flux;
      }
    }
  }

 protected:
  double calc_d2u_at(const TALYFEMLIB::ZEROPTV &pt) const {

#if(DIM == 3)
    return -3 * M_PI * M_PI * sin(M_PI * pt.x()) * sin(M_PI * pt.y()) * sin(M_PI * pt.z());
#elif(DIM == 2)
    return -2 * M_PI * M_PI * sin(M_PI * pt.x()) * sin(M_PI * pt.y());
#elif(DIM == 4)
    return -4 * M_PI * M_PI * sin(M_PI * pt.x()) * sin(M_PI * pt.y()) * sin(M_PI * pt.z()) * sin(M_PI * pt.t());
#endif
  }

  inline TALYFEMLIB::ZEROPTV calc_grad_u_at(const TALYFEMLIB::ZEROPTV &pt) const {
#if (DIM == 2)
    return TALYFEMLIB::ZEROPTV(M_PI * cos(M_PI * pt.x()) * sin(M_PI * pt.y()),
                               M_PI * sin(M_PI * pt.x()) * cos(M_PI * pt.y()),
                               0.0);
#elif (DIM == 3)
    return TALYFEMLIB::ZEROPTV(M_PI * cos(M_PI * pt.x()) * sin(M_PI * pt.y()) * sin(M_PI * pt.z()),
                               M_PI * sin(M_PI * pt.x()) * cos(M_PI * pt.y()) * sin(M_PI * pt.z()),
                               M_PI * sin(M_PI * pt.x()) * sin(M_PI * pt.y()) * cos(M_PI * pt.z()));
#endif
  }

};
#endif //DENDRITEKT_SSHTEQUATION_H
