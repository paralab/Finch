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
#include <talyfem/basis/constants.h>  // for constexpr_pow()

namespace TALYFEMLIB {

/** Used in calc_N and calc_dNde.
 * This is outside BoxLinearBasisImpl to give it external linkage
 * (which works because this file is eventually included by a .cpp file).
 * If this is moved inside BoxLinearBasisImpl, you will probably get
 * undefined symbol errors when you try to compile user code.
 */
    constexpr double epsilon_a[3][8] = {
            { -1, +1, -1, +1, -1, +1, -1, +1 },  // x
            { -1, -1, +1, +1, -1, -1, +1, +1 },  // y
            { -1, -1, -1, -1, +1, +1, +1, +1 }   // z
    };

#ifdef ENABLE_4D
constexpr double ep4d[4][16] = {
    {-1, +1, -1, +1, -1, +1, -1, +1, -1, +1, -1, +1, -1, +1, -1, +1},
    {-1, -1, +1, +1, -1, -1, +1, +1, -1, -1, +1, +1, -1, -1, +1, +1},
    {-1, -1, -1, -1, +1, +1, +1, +1, -1, -1, -1, -1, +1, +1, +1, +1},
    {-1, -1, -1, -1, -1, -1, -1, -1, +1, +1, +1, +1, +1, +1, +1, +1}
};
#endif

/**
 * Implementation of the linear "box" basis functions (1D line, 2D box, 3D hex).
 */
    template <int _nsd>
    class BoxLinearBasisImpl {
    public:
        static constexpr int nsd = _nsd;  ///< number of dimensions
        static constexpr int nbf = constexpr_pow(2, _nsd);  ///< number of shape funcs
        static constexpr int nbf_per_node = 1;  ///< shape funcs per node
        static constexpr double jacobian_scale = 1.0;  ///< multiplier for jacobian

        /**
         * Calculate shape functions.
         * @param localPt integration point
         * @param[out] n_out shape functions evaluated at localPt (output)
         */
#ifndef ENABLE_4D
        static void calc_N(const ZEROPTV& localPt, double (&n_out)[nbf]) {
          for (int bf = 0; bf < nbf; bf++) {
            n_out[bf] = 1.0;

            if (nsd >= 1)  // 1D component
              n_out[bf] *= (1 + epsilon_a[0][bf] * localPt.x());

            if (nsd >= 2)  // 2D component
              n_out[bf] *= (1 + epsilon_a[1][bf] * localPt.y());

            if (nsd >= 3)  // 3D component
              n_out[bf] *= (1 + epsilon_a[2][bf] * localPt.z());

            n_out[bf] /= nbf;
          }
        }
#else
  static void calc_N(const ZEROPTV& localPt, double (&n_out)[nbf]) {
          for (int bf = 0; bf < nbf; bf++) {
            n_out[bf] = 1.0;
            n_out[bf] *= (1 + ep4d[0][bf] * localPt.x());
            n_out[bf] *= (1 + ep4d[1][bf] * localPt.y());
            n_out[bf] *= (1 + ep4d[2][bf] * localPt.z());
            n_out[bf] *= (1 + ep4d[3][bf] * localPt.t());
            n_out[bf] /= nbf;
          }
        }
#endif

        /**
         * Calculate dNde values at integration point.
         * This could be simplified to use the nsd == 3 case all the time, but
         * that would rely on multiply-by-zero, which may not get optimized out
         * (due to NaN propagation).
         * @param localPt integration point
         * @param[out] dnde_out derivative of N in isoparametric space (output)
         */
#ifndef ENABLE_4D
        static void calc_dNde(const ZEROPTV& localPt, double (&dnde_out)[nbf][nsd]) {
          for (int bf = 0; bf < nbf; bf++) {
            for (int d = 0; d < nsd; d++) {
              switch (nsd) {
                case 1:
                  dnde_out[bf][d] = epsilon_a[0][bf] / nbf;
                      break;
                case 2:
                  switch (d) {  // 2D
                    case 0:
                      dnde_out[bf][d] = epsilon_a[0][bf] * (1 + epsilon_a[1][bf] * localPt.y()) / nbf;
                          break;
                    case 1:
                      dnde_out[bf][d] = (1 + epsilon_a[0][bf] * localPt.x()) * epsilon_a[1][bf] / nbf;
                          break;
                  }
                      break;
                case 3:
                  switch (d) {  // 3D
                    case 0:
                      dnde_out[bf][d] = epsilon_a[0][bf] * (1 + epsilon_a[1][bf] * localPt.y()) * (1 + epsilon_a[2][bf] * localPt.z()) / nbf;
                          break;
                    case 1:
                      dnde_out[bf][d] = (1 + epsilon_a[0][bf] * localPt.x()) * epsilon_a[1][bf] * (1 + epsilon_a[2][bf] * localPt.z()) / nbf;
                          break;
                    case 2:
                      dnde_out[bf][d] = (1 + epsilon_a[0][bf] * localPt.x()) * (1 + epsilon_a[1][bf] * localPt.y()) * epsilon_a[2][bf] / nbf;
                          break;
                  }
                      break;
              }
            }
          }
        }
#else
  static void calc_dNde(const ZEROPTV& localPt, double (&dnde_out)[nbf][nsd]) {
    for (int bf = 0; bf < nbf; bf++) {
      for (int d = 0; d < nsd; d++) {
        switch (d) {
          case 0:dnde_out[bf][d] = ep4d[0][bf] * (1 + ep4d[1][bf] * localPt.y()) * (1 + ep4d[2][bf] * localPt.z())
                * (1 + ep4d[3][bf] * localPt.t()) / nbf;
            break;
          case 1:dnde_out[bf][d] = (1 + ep4d[0][bf] * localPt.x()) * ep4d[1][bf] * (1 + ep4d[2][bf] * localPt.z())
                * (1 + ep4d[3][bf] * localPt.t()) / nbf;
            break;
          case 2:dnde_out[bf][d] = (1 + ep4d[0][bf] * localPt.x()) * (1 + ep4d[1][bf] * localPt.y()) * ep4d[2][bf]
                * (1 + ep4d[3][bf] * localPt.t()) / nbf;
            break;
          case 3:dnde_out[bf][d] =
                     (1 + ep4d[0][bf] * localPt.x()) * (1 + ep4d[1][bf] * localPt.y()) * (1 + ep4d[2][bf] * localPt.z())
                         * ep4d[3][bf] / nbf;
            break;
        }
      }
    }
  }
#endif
  /**
        * Calculate d2Nde values at integration point.
        * This could be simplified to use the nsd == 3 case all the time, but
        * that would rely on multiply-by-zero, which may not get optimized out
        * (due to NaN propagation).
        * @param localPt integration point
        * @param[out] d2nde_out derivative of N in isoparametric space (output)
        */
        static void calc_d2Nde(const ZEROPTV &localPt, double (&d2nde_out)[nbf][nsd * (nsd + 1) / 2]) {

          switch (nsd) {
            case 3: {

              /**
          { 0, 3, 4 },
          { 3, 1, 5 },
          { 4, 5, 2 },
           **/

              for (int bf = 0; bf < nbf; bf++) {

                /** d2N/(d2_xi) **/
                d2nde_out[bf][0] = 0;

                /** d2N/(d2_eta) **/
                d2nde_out[bf][1] = 0;

                /** d2N/(d2_zeta) **/
                d2nde_out[bf][2] = 0;

                /** d2N/(d_xi d_eta) **/
                d2nde_out[bf][3] = epsilon_a[0][bf]  * epsilon_a[1][bf] * (1 + epsilon_a[2][bf] * localPt.z()) / nbf;

                /** d2N/(d_xi d_zeta) **/
                d2nde_out[bf][4] = epsilon_a[0][bf] * (1 + epsilon_a[1][bf] * localPt.y()) * (epsilon_a[2][bf]) / nbf;

                /** d2N/(d_eta d_zeta) **/
                d2nde_out[bf][5] = (1 + epsilon_a[0][bf] * localPt.x()) * epsilon_a[1][bf] *  epsilon_a[2][bf] / nbf;

              }
              break;
            }
            case 2: {
              /**
              { 0, 1 }
              { 1, 2 }
              **/
              for (int bf = 0; bf < nbf; bf++) {


                /** d2N/(d2_xi) **/
                d2nde_out[bf][0] = 0;

                /** d2N/(d_xi d_eta) **/
                d2nde_out[bf][1] = epsilon_a[0][bf]*epsilon_a[1][bf]/nbf;

                /** d2N/(d2_eta ) **/
                d2nde_out[bf][2] = 0;

              }
              break;
            }
            case 1: {
              return;
            }
            default: {
//              std::cout << "d2Nde not implemented\n";
//              assert(false);

            }

          }
        }


    };

}  // namespace TALYFEMLIB
