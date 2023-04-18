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

#include <array>

#include <talyfem/grid/zeroptv.h>
#include <talyfem/basis/nonlinear.h>

namespace TALYFEMLIB {

/**
 * Maps tensor product-ordered N values calculated in NonlinearBasis to
 * the TalyFEM-ordered values.
 * Shape function order is like this (diagram taken from Gmsh manual):
 *
 *    3-----6-----2 
 *    |           | 
 *    |           | 
 *    7     8     5 
 *    |           | 
 *    |           | 
 *    0-----4-----1 
 */
static constexpr int quadratic_n_map[4][81] = {
  { 0, 1, 2 },  // 1D
  { 0, 1, 2, 3, 4, 5, 6, 7, 8 },  // 2D
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,  // 3D
    14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26 },

    {0, 1,  2, 3,  4,  5,  6, 7,  8, 9, 10, 11, 12, 13, //4D
          14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
          26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
          38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
          50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61,
          62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
          74, 75, 76, 77, 78, 79, 80}
};

/**
 * Implementation for the quadratic "box" basis functions (1D line, 2D box, 3D
 * hexahedron).
 */
template <int _nsd>
class BoxQuadraticBasisImpl {
 public:
  static constexpr int nsd = _nsd;  ///< number of spatial dimensions
  static constexpr int nbf = constexpr_pow(3, _nsd);  ///< number of shape funcs
  static constexpr int nbf_per_node = 1;  ///< shape funcs per node
  static constexpr double jacobian_scale = 1.0;  ///< multiplier for jacobian

  /**
   * Calculate shape functions.
   * @param localPt integration point
   * @param[out] n_out shape functions evaluated at localPt (output)
   */
  static void calc_N(const ZEROPTV& localPt, double (&n_out)[nbf]) {
    NonlinearBasisImpl<nsd, 3>::template calc_N<N_1D>(quadratic_n_map[nsd - 1], localPt, n_out);
  }

  /**
   * Calculate derivative of N in isoparametric space.
   * @param localPt integration point
   * @param[out] dnde_out derivative of N in isoparametric space (output)
   */
  static void calc_dNde(const ZEROPTV& localPt, double (&dnde_out)[nbf][nsd]) {
    NonlinearBasisImpl<nsd, 3>::template calc_dNde<N_1D, dN_1D>(quadratic_n_map[nsd - 1], localPt, dnde_out);
  }

  /**
   * Calculate second derivative of N in isoparametric space.
   * @param localPt integration point
   * @param[out] d2nde_out second derivative of N in isoparametric space
   */
  static void calc_d2Nde(const ZEROPTV& localPt, double (&d2nde_out)[nbf][nsd * (nsd + 1) / 2]) {
    NonlinearBasisImpl<nsd, 3>::template calc_d2Nde<N_1D, dN_1D, d2N_1D>(quadratic_n_map[nsd - 1], localPt, d2nde_out);
  }

 private:
  static std::array<double, 3> N_1D(double x) {
    return std::array<double, 3> {{
      0.5 * x * (x - 1),
      (1 - x * x),
      0.5 * x * (x + 1)
    }};
  }

  static std::array<double, 3> dN_1D(double x) {
    return std::array<double, 3> {{
      0.5 * (2 * x - 1),
      -2.0 * x,
      0.5 * (2 * x + 1)
    }};
  }

  static std::array<double, 3> d2N_1D(double x) {
    return std::array<double, 3> {{ 1, -2, 1 }};
  }
};

}  // namespace TALYFEMLIB
