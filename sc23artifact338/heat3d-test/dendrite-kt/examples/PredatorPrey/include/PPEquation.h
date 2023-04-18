/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //



#pragma once
#include <talyfem/fem/cequation.h>
#include "PPNodeData.h"

class PPEquation : public TALYFEMLIB::CEquation<PPNodeData> {
 public:

  explicit PPEquation()
      : TALYFEMLIB::CEquation<PPNodeData>(false, TALYFEMLIB::kAssembleGaussPoints) {}

  // never called, but required for PPEquation to be instantiable (pure virtual in CEquation)
  virtual void Solve(double dt, double t) { assert(false); }

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae) {
    // # of basis functions
    const int n_basis_functions = fe.nbf();
    // (determinant of J) cross W
    const double detJxW = fe.detJxW();
    // degree of freedom
    const int ndof = 2;

    // in order to assemble the gauss point, we loop over each pair of basis
    // functions and calculate the individual contributions to the 'Ae' matrix
    // and 'be' vector.

    for (int a = 0; a < n_basis_functions; a++) {
      for (int b = 0; b < n_basis_functions; b++) {
        double M = fe.N(a) * fe.N(b) * detJxW;

        Ae(a * ndof + 0, b * ndof + 0) += M * (1.0 / dt_ - 1.0);
        Ae(a * ndof + 0, b * ndof + 1) += -2 * M ;
        Ae(a * ndof + 1, b * ndof + 0) += 1 * M ;
        Ae(a * ndof + 1, b * ndof + 1) += M * (1.0 / dt_ - 1.0);
      }
    }

  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be) {
    // # of basis functions
    const int n_basis_functions = fe.nbf();
    // (determinant of J) cross W
    const double detJxW = fe.detJxW();
    // degree of freedom
    const int ndof = 2;

    const double u_pre_curr = p_data_->valueFEM(fe, U_PRE);
    const double v_pre_curr = p_data_->valueFEM(fe, V_PRE);

    for (int a = 0; a < n_basis_functions; a++) {
      be(a * ndof + 0) +=
          fe.N(a) / dt_ * u_pre_curr * detJxW ;
      be(a * ndof + 1) +=
          fe.N(a) / dt_ * v_pre_curr * detJxW ;
    }
  }

};
