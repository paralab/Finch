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

#include <talyfem/grid/zeroptv.h>

namespace TALYFEMLIB {

/**
 * Implementation for quadratic basis function for tetrahedral elements.
 */
class TetQuadraticBasisImpl {
 public:
  static constexpr int nsd = 3;  ///< number of spatial dimensions
  static constexpr int nbf = 10;  ///< number of shape functions
  static constexpr int nbf_per_node = 1;  ///< shape funcs per node
  static constexpr double jacobian_scale = 1.0 / 6.0;  ///< scale jacobian by

  /**
   * Calculate shape functions.
   * @param localPt integration point
   * @param[out] n_out shape functions evaluated at localPt (output)
   */
  static void calc_N(const ZEROPTV& localPt, double (&n_out)[nbf]) {
    n_out[1] = localPt.x()*(localPt.x()*2.0-1.0);
    n_out[2] = localPt.y()*(localPt.y()*2.0-1.0);
    n_out[3] = localPt.z()*(localPt.z()*2.0-1.0);
    n_out[0] = (localPt.x()+localPt.y()+localPt.z()-1.0)*(localPt.x()*2.0+localPt.y()*2.0+localPt.z()*2.0-1.0);
    n_out[5] = localPt.x()*localPt.y()*4.0;
    n_out[8] = localPt.y()*localPt.z()*4.0;
    n_out[9] = localPt.x()*localPt.z()*4.0;
    n_out[4] = localPt.x()*(localPt.x()+localPt.y()+localPt.z()-1.0)*-4.0;
    n_out[7] = localPt.z()*(localPt.x()+localPt.y()+localPt.z()-1.0)*-4.0;
    n_out[6] = -localPt.y()*(localPt.x()*4.0+localPt.y()*4.0+localPt.z()*4.0-4.0);
  }

  /**
   * Calculate derivative of N in isoparametric space.
   * @param localPt integration point
   * @param[out] dnde_out derivative of N in isoparametric space (output)
   */
  static void calc_dNde(const ZEROPTV& localPt, double (&dnde_out)[nbf][nsd]) {
    dnde_out[1][0] = localPt.x()*4.0-1.0;
    dnde_out[1][1] = 0.0;
    dnde_out[1][2] = 0.0;

    dnde_out[2][0] = 0.0;
    dnde_out[2][1] = localPt.y()*4.0-1.0;
    dnde_out[2][2] = 0.0;

    dnde_out[3][0] = 0.0;
    dnde_out[3][1] = 0.0;
    dnde_out[3][2] = localPt.z()*4.0-1.0;

    dnde_out[0][0] = localPt.x()*4.0+localPt.y()*4.0+localPt.z()*4.0-3.0;
    dnde_out[0][1] = localPt.x()*4.0+localPt.y()*4.0+localPt.z()*4.0-3.0;
    dnde_out[0][2] = localPt.x()*4.0+localPt.y()*4.0+localPt.z()*4.0-3.0;

    dnde_out[5][0] = localPt.y()*4.0;
    dnde_out[5][1] = localPt.x()*4.0;
    dnde_out[5][2] = 0.0;

    dnde_out[8][0] = 0.0;
    dnde_out[8][1] = localPt.z()*4.0;
    dnde_out[8][2] = localPt.y()*4.0;

    dnde_out[9][0] = localPt.z()*4.0;
    dnde_out[9][1] = 0.0;
    dnde_out[9][2] = localPt.x()*4.0;

    dnde_out[4][0] = localPt.x()*-8.0-localPt.y()*4.0-localPt.z()*4.0+4.0;
    dnde_out[4][1] = localPt.x()*-4.0;
    dnde_out[4][2] = localPt.x()*-4.0;

    dnde_out[7][0] = localPt.z()*-4.0;
    dnde_out[7][1] = localPt.z()*-4.0;
    dnde_out[7][2] = localPt.x()*-4.0-localPt.y()*4.0-localPt.z()*8.0+4.0;

    dnde_out[6][0] = localPt.y()*-4.0;
    dnde_out[6][1] = localPt.x()*-4.0-localPt.y()*8.0-localPt.z()*4.0+4.0;
    dnde_out[6][2] = localPt.y()*-4.0;
  }
};

}  // namespace TALYFEMLIB
