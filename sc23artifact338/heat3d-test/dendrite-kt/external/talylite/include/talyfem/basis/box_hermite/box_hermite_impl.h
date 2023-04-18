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
#include <talyfem/basis/box_linear/box_linear_impl.h>
#include <talyfem/basis/constants.h>  // for constexpr_pow()

namespace TALYFEMLIB {

/**
 * Maps tensor product-ordered N values calculated in NonlinearBasis to
 * the TalyFEM-ordered values.
 * This is outside BoxHermiteBasisImpl to give it external linkage
 * (which works because this file is eventually included by a .cpp file).
 * If this is moved inside BoxHermiteBasisImpl, you will probably get
 * undefined symbol errors when you try to compile user code.
 * - could be moved inside now that this is a templated class
 */
constexpr int hermite_n_map[4][256] = {
    { 0, 1, 2, 3 },
    { 0, 1, 4, 5, 2, 3, 6, 7, 12, 13, 8, 9, 14, 15, 10, 11 },
    { 0,1,8,9,2,4, 10, 12, 16, 17, 24, 25, 18, 20, 26, 28,3,5, 11, 13,6,7, 14,
      15, 19, 21, 27, 29, 22, 23, 30, 31, 32, 33, 40, 41, 34, 36, 42, 44, 48,
      49, 56, 57, 50, 52, 58, 60, 35, 37, 43, 45, 38, 39, 46, 47, 51, 53, 59,
      61, 54, 55, 62, 63 },
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
 * Implementation for hermite "box" basis functions (1D line, 2D box, 3D hex).
 */
template <int _nsd>
class BoxHermiteBasisImpl {
 public:
  static constexpr int nsd = _nsd;  ///< number of spatial dimensions
  static constexpr int nbf = constexpr_pow(4, _nsd);  ///< number of shape funcs
  /// Shape functions per node
  static constexpr int nbf_per_node = constexpr_pow(2, _nsd);
  static constexpr double jacobian_scale = 1.0;  ///< multiplier for jacobian

  /**
   * Calculate the global position of localPtv.
   * Currently uses linear basis function.
   * @param localPtv integration point
   * @param[in] n previously computed N values (unused)
   * @param elem node accessor
   * @param[out] pos_out global position of localPtv
   */
  static void calc_position(const ZEROPTV& localPtv, const double (&n)[nbf], const ElemNodes& elem, ZEROPTV* pos_out) {
    *pos_out = ZEROPTV(0, 0, 0);

/*
    // TODO hardcoded element size
    double hx = 1.0 / 10.0;
    for (int dir = 0; dir < nsd; dir++) {
      for (ElemNodeID a = 0; a < nbf; a++) {
        double h = (a % 2 == 0) ? hx : 0;
        (*pos_out)(dir) += n[a] * elem.get(a / nbf_per_node)->location()(dir) * h;
      }
    }*/

    // linear
    double n_linear[nbf / nbf_per_node];
    LinearBasis::calc_N(localPtv, n_linear);

    for (int dir = 0; dir < nsd; dir++) {
      for (int bf = 0; bf < (nbf / nbf_per_node); bf++) {
        (*pos_out)(dir) += n_linear[bf] * elem.node_pos(bf)(dir);
      }
    }

    // linear 1D
    // TODO HACK HACK HACK
    //const double n0_linear = (1 + -localPtv.x()) / 2.0;
    //const double n1_linear = (1 + localPtv.x()) / 2.0;
    //pos_out->x() = n0_linear * elem.get(0)->location().x() + n1_linear * elem.get(1)->location().x();
  }

  /**
   * Calculate shape functions.
   * @param localPt integration point
   * @param[out] n_out shape functions evaluated at localPt (output)
   */
  static void calc_N(const ZEROPTV& localPt, double (&n_out)[nbf]) {
    NonlinearBasisImpl<nsd, 4>::template calc_N<N_1D>(hermite_n_map[nsd - 1], localPt, n_out);
  }

  /**
   * Calculate derivative of N in isoparametric space.
   * @param localPt integration point
   * @param[out] dnde_out derivative of N in isoparametric space (output)
   */
  static void calc_dNde(const ZEROPTV& localPt, double (&dnde_out)[nbf][nsd]) {
    NonlinearBasisImpl<nsd, 4>::template calc_dNde<N_1D, dN_1D>(hermite_n_map[nsd - 1], localPt, dnde_out);
  }

  /**
   * Calculate the jacobian matrix.
   * Currently uses linear calculation.
   * @param localPtv integration point
   * @param[in] dnde derivative of N in isoparametric space at localPtv (unused)
   * @param elem node accessor
   * @param[out] dxde_out jacobian matrix
   */
  static void calc_dXde(const ZEROPTV& localPtv, const double(&dnde)[nbf][nsd], const ElemNodes& elem, double (&dxde_out)[nsd][nsd]) {
    //PrintInfo("Hermite calc_dXde, dNde[0][0]: ", dnde[0][0], ", dNde[1][0]: ", dnde[1][0], ", dNde[2][0]: ", dnde[2][0], ", dNde[3][0]: ", dnde[3][0]);

    // TODO hardcoded element size
    //double hx = 1.0 / 10.0;

    // linear dNde
    double dnde_linear[nbf / nbf_per_node][nsd];
    LinearBasis::calc_dNde(localPtv, dnde_linear);

    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        dxde_out[i][j] = 0;

        // 1D hermite
        /*for (ElemNodeID a = 0; a < nbf; a++) {
          double h = (a % 2 == 0) ? 1 : hx;
          dxde_out[i][j] += dnde[a][j] * elem.get(a / nbf_per_node)->location()(i) * h;
        }*/

        // linear generic
        for (int bf = 0; bf < (nbf / nbf_per_node); bf++) {
          dxde_out[i][j] += dnde_linear[bf][j] * elem.node_pos(bf)(i);
        }

        // 1D linear
        //dxde_out[i][j] += -0.5 * elem.get(0)->location()(i);
        //dxde_out[i][j] += +0.5 * elem.get(1)->location()(i);
      }
    }
  }

 private:
  typedef BoxLinearBasisImpl<_nsd> LinearBasis;

  static std::array<double, 4> N_1D(double x) {
    return std::array<double, 4> {{
                                      1.0 - (3.0/4.0)*(x+1)*(x+1) + (1.0/4.0)*(x+1)*(x+1)*(x+1),
                                      (1.0/2.0)*(x+1) - (1.0/2.0)*(x+1)*(x+1) + (1.0/8.0)*(x+1)*(x+1)*(x+1),
                                      1 - (1.0 - (3.0/4.0) * (x+1)*(x+1) + (1.0/4.0)*(x+1)*(x+1)*(x+1)),  // 1 - N(0)
                                      (1.0/8.0)*(x+1)*(x+1)*(x+1) - (1.0/4.0)*(x+1)*(x+1)
                                  }};
  }

  static std::array<double, 4> dN_1D(double x) {
    return std::array<double, 4> {{
                                      (-3.0/2.0)*(x+1) + (3.0/4.0) * (x+1)*(x+1),
                                      (1.0/2.0) - (x+1) + (3.0/8.0) * (x+1)*(x+1),
                                      -((-3.0/2.0)*(x+1) + (3.0/4.0) * (x+1)*(x+1)),  // -dN(0)
                                      (3.0/8.0)*(x+1)*(x+1) - (1.0/2.0)*(x+1)
                                  }};
  }
};

}  // namespace TALYFEMLIB
