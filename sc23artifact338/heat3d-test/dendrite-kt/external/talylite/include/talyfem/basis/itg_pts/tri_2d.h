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
 * N = 2 triangle gauss points (3 points).
 */
template<>
struct TriItgPts<2, 0> {
  static constexpr int n_itg_pts = 3;  ///< number of integration points
  ///! integration points
  static constexpr ZEROPTV itg_pts[3] = {
    /*ZEROPTV(0.1666666667, 0.7886751346, 0),
    ZEROPTV(0.6220084679, 0.2113248654, 0),
    ZEROPTV(0.4465819874e-1, 0.7886751346, 0),
    ZEROPTV(0.1666666667, 0.2113248654, 0)*/

    ZEROPTV(1.0 / 2.0, 0, 0),
    ZEROPTV(0, 1.0 / 2.0, 0),
    ZEROPTV(1.0 / 2.0, 1.0 / 2.0, 0),
  };
  ///! integration point weights
  static constexpr constexpr_array<double, 3> weights = {{
    1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0
  }};
};

/**
 * 7-point quadrature, PEREFECT for cubic elements.
 * Source: http://freeit.free.fr/Finite%20Element/Solin%20P.,%20Segeth%20K.,%20Dolezel%20I.,%20Higher-Order%20Finite%20Element%20Methods,%202004/C438X_04.pdf
 *
 * Note: this is currently being used by both quadratic and cubic triangles.
 * For some reason, it seems to give a much better result for quadratic triangles compared to the 6-point quadrature.
 * (the difference between 1e-2 error and 1e-6 error in the SSHT example).
 */
template<>
struct TriItgPts<3, 0> {
  static constexpr int n_itg_pts = 7;  ///< number of integration points
  ///! integration points
  static constexpr ZEROPTV itg_pts[7] = {
      ZEROPTV(0.333333333333334,0.333333333333334),
      ZEROPTV(0.470142064105115,0.470142064105115),
      ZEROPTV(0.470142064105115,0.059715871789770),
      ZEROPTV(0.059715871789770,0.470142064105115),
      ZEROPTV(0.101286507323456,0.101286507323456),
      ZEROPTV(0.101286507323456,0.797426985353087),
      ZEROPTV(0.797426985353087,0.101286507323456),
  };
  ///! integration point weights
  static constexpr constexpr_array<double, 7> weights = {{
     0.450000000000000/2,0.264788305577012/2,0.264788305577012/2,0.264788305577012/2,
     0.251878361089654/2,0.251878361089654/2,0.251878361089654/2
 }};
};

/**
 * 7-point quadrature, PEREFECT for cubic elements.
 * http://freeit.free.fr/Finite%20Element/Solin%20P.,%20Segeth%20K.,%20Dolezel%20I.,%20Higher-Order%20Finite%20Element%20Methods,%202004/C438X_04.pdf
 */
template<>
struct TriItgPts<4, 0> : public TriItgPts<3, 0> {};

/**
 * Left surface (X-) gauss points for triangle (surface ID -1).
 * Defined as a matrix transformation of the 1D line gauss points.
 */
template <int order>
struct TriItgPts<order, -1> {
  ///! number of integration points
  static constexpr int n_itg_pts = BoxItgPts<order, 1, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0.5, 0) * Matrix4::rotateZ(-90) * Matrix4::scale(0.5, 0, 0) * BoxItgPts<order, 1, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = array_mult(BoxItgPts<order, 1, 0>::weights, 0.5);
};

/**
 * Hypotenuse gauss points for triangle (surface ID +1).
 * Defined as a matrix transformation of the 1D line gauss points.
 */
template <int order>
struct TriItgPts<order, +1> {
  ///! number of integration points
  static constexpr int n_itg_pts = BoxItgPts<order, 1, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0.5, 0.5, 0)  // (0, 1/sqrt(2), 0)
                                                               * Matrix4::rotateZ(135) * Matrix4::scale(0.707106781186547524400844, 0, 0) * BoxItgPts<order, 1, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = array_mult(BoxItgPts<order, 1, 0>::weights, 0.707106781186547524400844);
};

/**
 * Bottom surface (Y-) gauss points for triangle (surface ID -2).
 * Defined as a matrix transformation of the 1D line gauss points.
 */
template <int order>
struct TriItgPts<order, -2> {
  ///! number of integration points
  static constexpr int n_itg_pts = BoxItgPts<order, 1, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0.5, 0, 0) * Matrix4::scale(0.5, 0, 0) * BoxItgPts<order, 1, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = array_mult(BoxItgPts<order, 1, 0>::weights, 0.5);
};

template <int order>
constexpr constexpr_array<ZEROPTV, TriItgPts<order, -1>::n_itg_pts> TriItgPts<order, -1>::itg_pts;
template <int order>
constexpr constexpr_array<double, TriItgPts<order, -1>::n_itg_pts> TriItgPts<order, -1>::weights;

template <int order>
constexpr constexpr_array<ZEROPTV, TriItgPts<order, +1>::n_itg_pts> TriItgPts<order, +1>::itg_pts;
template <int order>
constexpr constexpr_array<double, TriItgPts<order, +1>::n_itg_pts> TriItgPts<order, +1>::weights;

template <int order>
constexpr constexpr_array<ZEROPTV, TriItgPts<order, -2>::n_itg_pts> TriItgPts<order, -2>::itg_pts;
template <int order>
constexpr constexpr_array<double, TriItgPts<order, -2>::n_itg_pts> TriItgPts<order, -2>::weights;

}  // namespace TALYFEMLIB
