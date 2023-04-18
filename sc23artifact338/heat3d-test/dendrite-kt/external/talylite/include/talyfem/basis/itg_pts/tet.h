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
#include <talyfem/data_structures/constexpr_array.h>
#include <talyfem/data_structures/matrix4.h>
#include <talyfem/basis/constants.h>

namespace TALYFEMLIB {

/**
 * Volume gauss points and weights for a 3D tetrahedron for N = 2 (4 points).
 * These points are taken from N=2 of this paper:
 * http://www.sciencedirect.com/science/article/pii/0045782584900720
 * We don't have the analytic values for these, and it looks very
 * painful to calculate.
 */
template<>
struct TetItgPts<2, 0> {
  static constexpr int n_itg_pts = 4;  ///< number of integration points
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, 4> itg_pts = {{
    ZEROPTV(0.58541020, 0.13819660, 0.13819660),
    ZEROPTV(0.13819660, 0.58541020, 0.13819660),
    ZEROPTV(0.13819660, 0.13819660, 0.58541020),
    ZEROPTV(0.13819660, 0.13819660, 0.13819660)
  }};
  ///! integration point weights
  static constexpr constexpr_array<double, 4> weights = {{
    1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0
  }};
};

/**
 * Volume integration points and weights for a 3D tetrahedron for N = 3 (quadratic).
 */
template<>
struct TetItgPts<3, 0> {
  static constexpr int n_itg_pts = 8;  ///< number of integration points
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, 8> itg_pts = {{
     ZEROPTV(0.01583591,  0.328054697, 0.328054697),
     ZEROPTV(0.328054697, 0.01583591, 0.328054697),
     ZEROPTV(0.328054697, 0.328054697, 0.01583591),
     ZEROPTV(0.328054697, 0.328054697, 0.328054697),
     ZEROPTV(0.679143178, 0.106952274, 0.106952274),
     ZEROPTV(0.106952274, 0.679143178, 0.106952274),
     ZEROPTV(0.106952274, 0.106952274, 0.679143178),
     ZEROPTV(0.106952274, 0.106952274, 0.106952274)
 }};

  ///! integration point weights
  static constexpr constexpr_array<double, 8> weights = {{
      0.023087995*6, 0.023087995*6, 0.023087995*6, 0.023087995*6,
      0.018578672*6, 0.018578672*6, 0.018578672*6, 0.018578672*6
  }};
};

/**
 * Volume integration points and weights for a 3D tetrahedron for N = 4 (cubic).
 */
template<>
struct TetItgPts<4, 0> {
  static constexpr int n_itg_pts = 10;  ///< number of integration points

  ///! integration points
  static constexpr constexpr_array<ZEROPTV, 10> itg_pts = {{
    ZEROPTV(0.7784952948213300,     0.0738349017262234,     0.0738349017262234),
    ZEROPTV(0.0738349017262234,     0.7784952948213300,     0.0738349017262234),
    ZEROPTV(0.0738349017262234,     0.0738349017262234,     0.7784952948213300),
    ZEROPTV(0.0738349017262234,     0.0738349017262234,     0.0738349017262234),
    ZEROPTV(0.4062443438840510,     0.4062443438840510,     0.0937556561159491),
    ZEROPTV(0.4062443438840510,     0.0937556561159491,     0.4062443438840510),
    ZEROPTV(0.4062443438840510,     0.0937556561159491,     0.0937556561159491),
    ZEROPTV(0.0937556561159491,     0.4062443438840510,     0.4062443438840510),
    ZEROPTV(0.0937556561159491,     0.4062443438840510,     0.0937556561159491),
    ZEROPTV(0.0937556561159491,     0.0937556561159491,     0.4062443438840510)
  }};

  ///! integration point weights
  static constexpr constexpr_array<double, 10> weights = {{
    0.0476331348432089,0.0476331348432089,0.0476331348432089,0.0476331348432089,
    0.1349112434378610,0.1349112434378610,0.1349112434378610,0.1349112434378610,
    0.1349112434378610,0.1349112434378610
  }};
};

/**
 * Hypotenuse surface gauss points for 3D tetrahedron (surface ID 1).
 * Defined as a 4x4 matrix transformation of the 2D triangle gauss points.
 */
template <int order>
struct TetItgPts<order, 1> {  // whatever the other side is called (hypotenuse?)
  //! number of integration points
  static constexpr int n_itg_pts = TriItgPts<order, 0>::n_itg_pts;

  /**
   * Integration points.
   * Matrix found by transforming the following set of points:
   *  (0, 0, 0), (1, 0, 0), (0, 1, 0)
   * to
   *  (0, 0, 1), (1, 0, 0), (0, 1, 0)
   * Using the Python script found here: http://stackoverflow.com/a/27547597
   */
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4 {{
      1, 0, 1, 0,
      0, 1, 1, 0,
      -1, -1, 1, 1,
      0, 0, 0, 1 }} * TriItgPts<order, 0>::itg_pts;

  // barycentric quadratic
  /*static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = {{
    ZEROPTV(0.333333333333332, 0.333333333333334, 0.333333333333334),
    ZEROPTV(0.059715871789770, 0.470142064105115, 0.470142064105115),
    ZEROPTV(0.470142064105115, 0.470142064105115, 0.059715871789770),
    ZEROPTV(0.470142064105115, 0.059715871789770, 0.470142064105115),
    ZEROPTV(0.797426985353088, 0.101286507323456, 0.101286507323456),
    ZEROPTV(0.101286507323457, 0.101286507323456, 0.797426985353087),
    ZEROPTV(0.101286507323457, 0.797426985353087, 0.101286507323456)
  }};*/

  // pengfei transformation quadratic
  /*static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = {{
    ZEROPTV(0.333333333333334, 0.333333333333334, 0.333333333333332),
    ZEROPTV(0.470142064105115, 0.470142064105115, 0.059715871789770),
    ZEROPTV(0.470142064105115, 0.059715871789770, 0.470142064105115),
    ZEROPTV(0.059715871789770, 0.470142064105115, 0.470142064105115),
    ZEROPTV(0.101286507323456, 0.101286507323456, 0.797426985353088),
    ZEROPTV(0.101286507323456, 0.797426985353087, 0.101286507323457),
    ZEROPTV(0.797426985353087, 0.101286507323456, 0.101286507323457)
  }};*/


  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = TriItgPts<order, 0>::weights;

  /*static constexpr constexpr_array<double, n_itg_pts> weights = {{  // ignored
     0.450000000000000/2,0.264788305577012/2,0.264788305577012/2,0.264788305577012/2,
     0.251878361089654/2,0.251878361089654/2,0.251878361089654/2
 }};*/
};

/**
 * YZ-plane surface gauss points for 3D tetrahedron (surface ID 2).
 * Defined as a 4x4 matrix transformation of the 2D triangle gauss points.
 */
template <int order>
struct TetItgPts<order, 2> {  // YZ-plane surface
  //! number of integration points
  static constexpr int n_itg_pts = TriItgPts<order, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::rotateY(-90) * TriItgPts<order, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = TriItgPts<order, 0>::weights;
};

/**
 * XZ-plane surface gauss points for 3D tetrahedron (surface ID 3).
 * Defined as a 4x4 matrix transformation of the 2D triangle gauss points.
 */
template <int order>
struct TetItgPts<order, 3> {  // XZ-plane surface
  //! number of integration points
  static constexpr int n_itg_pts = TriItgPts<order, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::rotateX(90) * TriItgPts<order, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = TriItgPts<order, 0>::weights;
};


//! Used to calculate tet surface 4 GPs (mirror along Y=X)
namespace detail2 {
//! Used by swap_helper
template <std::size_t... Is>
struct indices {};

//! Used by swap_helper
template <std::size_t N, std::size_t... Is>
struct build_indices: build_indices<N-1, N-1, Is...> {
};

//! Used by swap_helper
template <std::size_t... Is>
struct build_indices<0, Is...>: indices<Is...> {
  //! work-around for Intel constexpr conversion bug
  typedef indices<Is...> base_indices;
};

/**
 * @param arr array to swap
 * @returns arr with X and Y values swapped
 */
template<size_t N, size_t... Is>
constexpr constexpr_array<ZEROPTV, N> swap_helper(const ZEROPTV (&arr)[N], indices<Is...>) {
  return {{ ZEROPTV(arr[Is][1], arr[Is][0])... }};
}
}  // namespace detail

/**
 * XY-plane surface gauss points for 3D tetrahedron (surface ID 4).
 * "Mirrored" version of the 2D triangle gauss points.
 */
template <int order>
struct TetItgPts<order, 4> {
  //! number of integration points
  static constexpr int n_itg_pts = TriItgPts<order, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = detail2::swap_helper(TriItgPts<order, 0>::itg_pts, detail2::build_indices<n_itg_pts>());
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = TriItgPts<order, 0>::weights;
};


#define DEF_TET_SURF(surf_id) \
template<int order> \
constexpr constexpr_array<ZEROPTV, TetItgPts<order, surf_id>::n_itg_pts> TetItgPts<order, surf_id>::itg_pts; \
template<int order> \
constexpr constexpr_array<double, TetItgPts<order, surf_id>::n_itg_pts> TetItgPts<order, surf_id>::weights;

DEF_TET_SURF(1)
DEF_TET_SURF(2)
DEF_TET_SURF(3)
DEF_TET_SURF(4)

}  // namespace TALYFEMLIB
