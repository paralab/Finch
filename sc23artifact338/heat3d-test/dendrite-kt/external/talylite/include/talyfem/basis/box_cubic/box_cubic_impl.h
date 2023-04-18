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
 */
static constexpr int cubic_n_map[4][256] = {
    { 0, 2, 3, 1 },  // 1D
    { 0, 4, 5, 1, 11, 12, 13, 6, 10, 15, 14, 7, 3, 9, 8, 2 }, // 2D
    { 0, 16, 17, 1, 23, 24, 25, 18, 22, 27, 26, 19, 3, 21, 20, 2, 8, 40, 41, 9,
      47, 48, 49, 42, 46, 51, 50, 43, 11, 45, 44, 10, 12, 52, 53, 13, 59, 60,
      61, 54, 58, 63, 62, 55, 15, 57, 56, 14, 4, 28, 29, 5, 35, 36, 37, 30, 34,
      39, 38, 31, 7, 33, 32, 6 },  // 3D
    {0, 1,  2, 3,  4,  5,  6, 7,  8, 9, 10, 11, 12, 13, //4D
     14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
     26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
     38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
     50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61,
     62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
     74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85,
     86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97,
     98, 99, 100, 101, 102, 103, 104, 105, 106, 107,
     108, 109, 110, 111, 112, 113, 114, 115, 116, 117,
     118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
     128, 129, 130, 131, 132, 133, 134, 135, 136, 137,
     138, 139, 140, 141, 142, 143, 144, 145, 146, 147,
     148, 149, 150, 151, 152, 153, 154, 155, 156, 157,
     158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
     168, 169, 170, 171, 172, 173, 174, 175, 176, 177,
     178, 179, 180, 181, 182, 183, 184, 185, 186, 187,
     188, 189, 190, 191, 192, 193, 194, 195, 196, 197,
     198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
     208, 209, 210, 211, 212, 213, 214, 215, 216, 217,
     218, 219, 220, 221, 222, 223, 224, 225, 226, 227,
     228, 229, 230, 231, 232, 233, 234, 235, 236, 237,
     238, 239, 240, 241, 242, 243, 244, 245, 246, 247,
     248, 249, 250, 251, 252, 253, 254, 255}

};

/**
 * Implementation for the cubic "box" basis functions (1D line, 2D box, 3D
 * hexahedron).
 */
template <int _nsd>
class BoxCubicBasisImpl {
 public:
  static constexpr int nsd = _nsd;  ///< number of spatial dimensions
  static constexpr int nbf = constexpr_pow(4, _nsd);  ///< number of shape funcs
  static constexpr int nbf_per_node = 1;  ///< shape funcs per node
  static constexpr double jacobian_scale = 1.0;  ///< multiplier for jacobian

  /**
   * Calculate shape functions.
   * @param localPt integration point
   * @param[out] n_out shape functions evaluated at localPt (output)
   */
  static void calc_N(const ZEROPTV& localPt, double (&n_out)[nbf]) {
    NonlinearBasisImpl<nsd, 4>::template calc_N<N_1D>(cubic_n_map[nsd - 1], localPt, n_out);
  }

  /**
   * Calculate derivative of N in isoparametric space.
   * @param localPt integration point
   * @param[out] dnde_out derivative of N in isoparametric space (output)
   */
  static void calc_dNde(const ZEROPTV& localPt, double (&dnde_out)[nbf][nsd]) {
    NonlinearBasisImpl<nsd, 4>::template calc_dNde<N_1D, dN_1D>(cubic_n_map[nsd - 1], localPt, dnde_out);
  }

  /**
   * Calculate second derivative of N in isoparametric space.
   * @param localPt integration point
   * @param[out] d2nde_out second derivative of N in isoparametric space
   */
  static void calc_d2Nde(const ZEROPTV& localPt, double (&d2nde_out)[nbf][nsd * (nsd + 1) / 2]) {
    NonlinearBasisImpl<nsd, 4>::template calc_d2Nde<N_1D, dN_1D, d2N_1D>(cubic_n_map[nsd - 1], localPt, d2nde_out);
  }

 private:
  static std::array<double, 4> N_1D(double x) {
    return std::array<double, 4> {{
                                      (-9. / 16.) * (x * x * x - x * x - (1. / 9.) * x + (1. / 9.)),
                                      (27. / 16.) * (x * x * x - (1. / 3.) * x * x - x + (1. / 3.)),
                                      (-27. / 16.) * (x * x * x + (1. / 3.) * x * x - x - (1. / 3.)),
                                      (9. / 16.) * (x * x * x + x * x - (1. / 9.) * x - (1. / 9.))
                                  }};
  }

  static std::array<double, 4> dN_1D(double x) {
    return std::array<double, 4> {{
                                      (-9. / 16.) * (3 * x * x - 2 * x - (1. / 9.)),
                                      (27. / 16.) * (3 * x * x - (2. / 3.) * x - 1),
                                      (-27. / 16.) * (3 * x * x + (2. / 3.) * x - 1),
                                      (9. / 16.) * (3 * x * x + 2 * x - (1. / 9.))
                                  }};
  }

  // not verified - blindly copied from old library
  static std::array<double, 4> d2N_1D(double x) {
    return std::array<double, 4> {{
                                      (-9. / 16.) * (6 * x - 2),
                                      (27. / 16.) * (6 * x - (2. / 3.)),
                                      (-27. / 16.) * (6 * x + (2. / 3.)),
                                      (9. / 16.) * (6 * x + 2)
                                  }};
  }
};

}  // namespace TALYFEMLIB
