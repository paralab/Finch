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
 * Volume gauss points and weights for a 4D pentatope for N = 2 (5 points).
 * These points are taken from N=2 of this paper:
 * http://onlinelibrary.wiley.com/doi/10.1002/fld.1796/full
 * We don't have the analytic values for these, and it looks very
 * painful to calculate.
 */
template<>
struct PentItgPts<2, 0> {
  static constexpr int n_itg_pts = 5;  ///< number of integration points
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, 5> itg_pts = {{
    ZEROPTV(0.118350341907227374, 0.118350341907227374, 0.118350341907227374, 0.118350341907227374),
    ZEROPTV(0.526598632371090503, 0.118350341907227374, 0.118350341907227374, 0.118350341907227374),
    ZEROPTV(0.118350341907227374, 0.526598632371090503, 0.118350341907227374, 0.118350341907227374),
    ZEROPTV(0.118350341907227374, 0.118350341907227374, 0.526598632371090503, 0.118350341907227374),
    ZEROPTV(0.118350341907227374, 0.118350341907227374, 0.118350341907227374, 0.526598632371090503)
    
  }};
  ///! integration point weights
  static constexpr constexpr_array<double, 5> weights = {{
    1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0
  }};
};

//
// TODO(4D) Find quadrature points for higher order elements.
//
// template<>
// struct TetItgPts<3, 0> {
//   static constexpr int n_itg_pts = 10;  ///< number of integration points
//   ///! integration points
//   static constexpr constexpr_array<ZEROPTV, 10> itg_pts = {{
//     ZEROPTV(0.7784952948213300,     0.0738349017262234,     0.0738349017262234),
//     ZEROPTV(0.0738349017262234,     0.7784952948213300,     0.0738349017262234),
//     ZEROPTV(0.0738349017262234,     0.0738349017262234,     0.7784952948213300),
//     ZEROPTV(0.0738349017262234,     0.0738349017262234,     0.0738349017262234),
//     ZEROPTV(0.4062443438840510,     0.4062443438840510,     0.0937556561159491),
//     ZEROPTV(0.4062443438840510,     0.0937556561159491,     0.4062443438840510),
//     ZEROPTV(0.4062443438840510,     0.0937556561159491,     0.0937556561159491),
//     ZEROPTV(0.0937556561159491,     0.4062443438840510,     0.4062443438840510),
//     ZEROPTV(0.0937556561159491,     0.4062443438840510,     0.0937556561159491),
//     ZEROPTV(0.0937556561159491,     0.0937556561159491,     0.4062443438840510)
//   }};
//   ///! integration point weights
//   static constexpr constexpr_array<double, 10> weights = {{
//     0.0476331348432089,0.0476331348432089,0.0476331348432089,0.0476331348432089,0.1349112434378610,0.1349112434378610,0.1349112434378610,0.1349112434378610,0.1349112434378610,0.1349112434378610
//   }};
// };

// template<>
// struct TetItgPts<4, 0> : public TetItgPts<3, 0> {};

/**
 * Hypotenuse surface gauss points for 3D tetrahedron (surface ID 1).
 * Defined as a 4x4 matrix transformation of the 2D triangle gauss points.
 */
// template <int order>
// struct TetItgPts<order, 1> {  // whatever the other side is called (hypotenuse?)
//   //! number of integration points
//   static constexpr int n_itg_pts = TriItgPts<order, 0>::n_itg_pts;

//   /**
//    * Integration points.
//    * Matrix found by transforming the following set of points:
//    *  (0, 0, 0), (1, 0, 0), (0, 1, 0)
//    * to
//    *  (0, 0, 1), (1, 0, 0), (0, 1, 0)
//    * Using the Python script found here: http://stackoverflow.com/a/27547597
//    */
//   static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4 {{
//       1, 0, 1, 0,
//       0, 1, 1, 0,
//       -1, -1, 1, 1,
//       0, 0, 0, 1 }} * TriItgPts<order, 0>::itg_pts;
//   ///! integration point weights
//   static constexpr constexpr_array<double, n_itg_pts> weights = TriItgPts<order, 0>::weights;
// };

/**
 * YZ-plane surface gauss points for 3D tetrahedron (surface ID 2).
 * Defined as a 4x4 matrix transformation of the 2D triangle gauss points.
 */
// template <int order>
// struct TetItgPts<order, 2> {  // YZ-plane surface
//   //! number of integration points
//   static constexpr int n_itg_pts = TriItgPts<order, 0>::n_itg_pts;
//   ///! integration points
//   static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::rotateY(-90) * TriItgPts<order, 0>::itg_pts;
//   ///! integration point weights
//   static constexpr constexpr_array<double, n_itg_pts> weights = TriItgPts<order, 0>::weights;
// };

/**
 * XZ-plane surface gauss points for 3D tetrahedron (surface ID 3).
 * Defined as a 4x4 matrix transformation of the 2D triangle gauss points.
 */
// template <int order>
// struct TetItgPts<order, 3> {  // XZ-plane surface
//   //! number of integration points
//   static constexpr int n_itg_pts = TriItgPts<order, 0>::n_itg_pts;
//   ///! integration points
//   static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::rotateX(90) * TriItgPts<order, 0>::itg_pts;
//   ///! integration point weights
//   static constexpr constexpr_array<double, n_itg_pts> weights = TriItgPts<order, 0>::weights;
// };

/**
 * XY-plane surface gauss points for 3D tetrahedron (surface ID 4).
 * Identical to the 2D triangle gauss points.
 */
// template <int order>  // XY-plane surface
// struct TetItgPts<order, 4> : public TriItgPts<order, 0> {  // Y-plane surface, conveniently identical to triangle
// };


// #define DEF_TET_SURF(surf_id) \
// template<int order> \
// constexpr constexpr_array<ZEROPTV, TetItgPts<order, surf_id>::n_itg_pts> TetItgPts<order, surf_id>::itg_pts; \
// template<int order> \
// constexpr constexpr_array<double, TetItgPts<order, surf_id>::n_itg_pts> TetItgPts<order, surf_id>::weights;

// DEF_TET_SURF(1)
// DEF_TET_SURF(2)
// DEF_TET_SURF(3)
// no surface for 4, since it's just normal triangle

}  // namespace TALYFEMLIB
