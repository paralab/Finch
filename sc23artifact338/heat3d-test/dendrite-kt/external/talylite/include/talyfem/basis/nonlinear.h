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

#include <assert.h>

#include <array>

namespace TALYFEMLIB {

/**
 * Helper functions for calculating N/dNde/d2Nde values for non-linear basis
 * functions (quadratic, cubic, ...) for boxes.
 *
 * It's easy to generate higher dimension N values using a tensor product of
 * the 1D gauss points:
 *   N[0] = n_1d_x[0] * n_1d_y[0];
 *   N[1] = n_1d_x[0] * n_1d_y[1];
 *   N[2] = n_1d_x[1] * n_1d_y[0];
 *   N[3] = n_1d_x[1] * n_1d_y[1];
 *
 * This is simple to calculate with two nested loops (or 3 in 3D):
 *   for (int j = 0; j < n_1d_shape_funcs; j++) {
 *     for (int i = 0; i < n_1d_shape_funcs; i++) {
 *       double val = n_1d_x[i] * n_1d_y[j];
 *       // store val in the N array somewhere...
 *     }
 *   }
 *
 * The hard part is ordering these values so that they match the order of
 * the element connectivity ordering. For this, we use the concept of an
 * "N map."
 *
 * This is a map that defines which shape function an iteration of the nested
 * loops (i, j) corresponds to. We flatten the (i, j) indices into 1D space
 * with the classic trick:
 *   int idx = i * n_1d_shape_funcs + j;
 *
 * The "N map" then maps this index to a shape function index.
 *   int N_idx = n_map[i * n_1d_shape_funcs + j];
 * 
 * The complete 2D example looks something like this:
 *   for (int j = 0; j < n_1d_shape_funcs; j++) {
 *     for (int i = 0; i < n_1d_shape_funcs; i++) {
 *       double val = n_1d_x[i] * n_1d_y[j];
 *
 *       int n_map_idx = i * n_1d_shape_funcs + j;
 *       int out_idx = n_map[n_map_idx];
 *
 *       N[out_idx] = val;
 *     }
 *   }
 *
 * Again - the N map is defined by the connectivity order, which depends on
 * how the mesh is generated/loaded. The mesh generators in Grid/grid_types/
 * and the reorder function of the Gmsh loader in FileIO/Gmsh_ascii.cpp are
 * probably the best place to look for diagrams.
 */
template <int nsd, int order>
struct NonlinearBasisImpl {
  static const int nbf = constexpr_pow(order, nsd);  ///< number of shape funcs

  //! type of a function that calculates N in 1 dimension
  typedef std::array<double, order> (* n_1d_func_t)(double x);

  //! type of function that calculates dN in 1 dimension
  typedef std::array<double, order> (* dn_1d_func_t)(double x);

  //! type of function that calculates d2N in 1 dimension
  typedef std::array<double, order> (* d2n_1d_func_t)(double x);

  //! type of the n map array
  typedef const int (&n_map_t)[order * order * order * order];

  /**
   * Calculate the shape function values at an integration point.
   * @param[in] n_map map from the tensor product order to the TalyFEM
   *                  connectivity order
   * @param localPt integration point
   * @param[out] n_out shape function evaluated at localPt (output)
   */
  template <n_1d_func_t n_1d_func>
  inline static void calc_N(n_map_t n_map, const ZEROPTV &localPt,
                            double (&n_out)[nbf]) {
    std::array<std::array<double, order>, nsd> arr;
    for (int dim = 0; dim < nsd; dim++) {
      arr[dim] = n_1d_func(localPt(dim));
    }

    // the outer loop will only do anything if nsd >= 3,
    // the secondmost outer loop will only do anything if nsd >= 2, etc.
    for (int l = 0; l < (nsd >= 4 ? order : 1); l++) {
      for (int k = 0; k < (nsd >= 3 ? order : 1); k++) {
        for (int j = 0; j < (nsd >= 2 ? order : 1); j++) {
          for (int i = 0; i < (nsd >= 1 ? order : 1); i++) {
            const int map_idx = order * order * order * l + order * order * k + order * j + i;
            const int id = n_map[map_idx];
            // 1D
            n_out[id] = arr[0][i];

            if (nsd >= 2)  // 2D
              n_out[id] *= arr[1][j];

            if (nsd >= 3)  // 3D
              n_out[id] *= arr[2][k];

            if (nsd >= 4)  // 4D
              n_out[id] *= arr[3][l];
          }
        }
      }
    }
  }


  // constant
  /**
   * Calculate the derivative of the shape functions in isoparametric space.
   * @param[in] n_map map from tensor product order to TalyFEM order
   * @param localPt integration point
   * @param[out] dnde_out dN/de output
   */
  template <n_1d_func_t n_1d_func, dn_1d_func_t dn_1d_func>
  inline static void calc_dNde(n_map_t n_map, const ZEROPTV &localPt,
                               double (&dnde_out)[nbf][nsd]) {
    // N_1D values
    std::array<std::array<double, order>, nsd> n_1d;
    for (int dim = 0; dim < nsd; dim++) {
      n_1d[dim] = n_1d_func(localPt(dim));
    }

    // dN_1D values
    std::array<std::array<double, order>, nsd> dn_1d;
    for (int dim = 0; dim < nsd; dim++) {
      dn_1d[dim] = dn_1d_func(localPt(dim));
    }

    // the outer loop will only do anything if nsd >= 3,
    // the secondmost outer loop will only do anything if nsd >= 2, etc.
    for (int l = 0; l < (nsd >= 4 ? order : 1); l++) {
      for (int k = 0; k < (nsd >= 3 ? order : 1); k++) {
        for (int j = 0; j < (nsd >= 2 ? order : 1); j++) {
          for (int i = 0; i < (nsd >= 1 ? order : 1); i++) {
            const int map_idx = order * order * order * l + order * order * k + order * j + i;
            const int id = n_map[map_idx];

            // This is what we calculate in the 3D case.
            // dnde_out[id][0] = dn_1d[0][i] * n_1d[1][j] * n_1d[2][k];
            // dnde_out[id][1] = n_1d[0][i] * dn_1d[1][j] * n_1d[2][k];
            // dnde_out[id][2] = n_1d[0][i] * n_1d[1][j] * dn_1d[2][k];
            // The code below will handle 1D, 2D, and 3D.


            dnde_out[id][0] = dn_1d[0][i];

            if (nsd >= 2) {
              dnde_out[id][0] *= n_1d[1][j];
              dnde_out[id][1] = n_1d[0][i] * dn_1d[1][j];
            }
            if (nsd >= 3) {
              dnde_out[id][0] *= n_1d[2][k];
              dnde_out[id][1] *= n_1d[2][k];
              dnde_out[id][2] = n_1d[0][i] * n_1d[1][j] * dn_1d[2][k];
            }

            if (nsd >= 4) {
              dnde_out[id][0] *= n_1d[3][l];
              dnde_out[id][1] *= n_1d[3][l];
              dnde_out[id][2] *= n_1d[3][l];
              dnde_out[id][3] = n_1d[0][i] * n_1d[1][j] * n_1d[2][k] * dn_1d[3][l];
            }
          }
        }
      }
    }
  }

  //! Dimension of dXde2 (second derivative's jacobian matrix)/cof2.
  static const int d2_dim = nsd * (nsd + 1) / 2;

  /**
   * Calculate the second derivative of the shape functions in isoparametric
   * space.
   * @param localPt integration point
   * @param n_map map of which tensor product order shape function goes where
   * @param[out] d2nde_out second derivative of N in insoparametric space
   */
  template <n_1d_func_t n_1d_func, dn_1d_func_t dn_1d_func, d2n_1d_func_t d2n_1d_func>
  // NOLINT(whitespace/line_length)
  inline static void calc_d2Nde(n_map_t n_map, const ZEROPTV &localPt,
                                double (&d2nde_out)[nbf][d2_dim]) {

    // n_values[deriv][axis][bf_1d] = deriv-th derivative of the 1D basis function bf_1d along axis
    // so n_values[0][0][0] = N evaluated at localPt.x() for basis function 0
    //    n_values[1][2][3] = dN evaluated at localPt.z() for basis function 3
    //    n_values[2][1][0] = d2N evaluated at localPt.y() for basis function 0
    std::array<std::array<std::array<double, order>, nsd>, 3> n_values;

    // calculate N_1D values
    for (int dim = 0; dim < nsd; dim++) {
      n_values[0][dim] = n_1d_func(localPt(dim));
    }

    // calculate dN_1D values
    for (int dim = 0; dim < nsd; dim++) {
      n_values[1][dim] = dn_1d_func(localPt(dim));
    }

    // calculate d2N_1D values
    for (int dim = 0; dim < nsd; dim++) {
      n_values[2][dim] = d2n_1d_func(localPt(dim));
    }


    // fill in d2nde by iterating through all the nodes in this order:
    /*
      ^
      | y

      6 - 7 - 8
      |       |
      3 - 4 - 5
      |       |
      0 - 1 - 2    -> x
    */
    // (0, 0), (1, 0), (2, 0), then the next row up (0, 1), (1, 1), (2, 1), etc
    // note that this is NOT the normal order of nodes in an element!
    // i DOES NOT directly correspond to basis function i - we need to use a map
    for (int i = 0; i < nbf; i++) {
      int bf = n_map[i];
      int bf_1d[4] = {
          i % order,
          (i / order) % order,
          ((i / order) / order) % order,
          (((i / order) / order) / order)
      };

      // loop over the possible derivative combinations (x*x, x*y, y*y, ...)
      // what derivative j corresponds to is defined by the deriv_level function
      for (int j = 0; j < d2_dim; j++) {

        // calculate d2nde[bf][j] by multiplying together values in n_values
        // (the number of values to multiply depends on nsd)
        double val = 1.0;
        for (int axis = 0; axis < nsd; axis++) {
          int deriv = deriv_level(j, axis);
          val *= n_values[deriv][axis][bf_1d[axis]];
        }

        d2nde_out[bf][j] = val;
      }
    }
  }

 private:
  /**
   * This controls the order of the derivatives in d2Nde.
   * @param d derivative combination
   * @param axis which axis (x/y/z)
   * @returns which derivative level corresponds to (d, axis)
   *          0 corresponds to N, 1 corresponds to dNde,
   *          and 2 corresponds to d2Nde (in 1D)
   */
  inline static int deriv_level(int d, int axis) {
    if (nsd == 1) {
      return 2;  // x^2
    } else if (nsd == 2) {
      // order: x^2, xy, y^2
      static const int map[3][2]{
          {2, 0},
          {1, 1},
          {0, 2},
      };
      return map[d][axis];
    } else if (nsd == 3) {
      // order: x^2, y^2, z^2,
      // xy, xz, yz
      static const int map[6][3]{
          {2, 0, 0},
          {0, 2, 0},
          {0, 0, 2},
          {1, 1, 0},
          {1, 0, 1},
          {0, 1, 1},
      };
      return map[d][axis];
    } else if(nsd == 4){
      // order: x^2, y^2, z^2,
      // xy, xz, yz
      static const int map[10][4]{
          {2, 0, 0, 0},
          {0, 2, 0, 0},
          {0, 0, 2, 0},
          {0, 0, 0, 2},
          {1, 1, 0, 0},
          {1, 0, 1, 0},
          {1, 0, 0, 1},
          {0, 1, 1, 0},
          {0, 1, 0, 1},
          {0, 0, 1, 1},
      };
      return map[d][axis];
    }
  }
};

}  // namespace TALYFEMLIB
