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

#include <talyfem/data_structures/matrix4.h>

namespace TALYFEMLIB {

// 1D surfaces (always one point)
/**
 * Left 1D surface gauss point (surface ID -1).
 */
template<int order>
struct BoxItgPts<order, 1, -1>  {
  static constexpr int n_itg_pts = 1;  ///< number of integration points
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = {{ ZEROPTV(-1, 0, 0) }};
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = {{ 1 }};
};

/**
 * Right 1D surface gauss point (surface ID +1).
 */
template<int order>
struct BoxItgPts<order, 1, +1> {
  static constexpr int n_itg_pts = 1;  ///< number of integration points
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = {{ ZEROPTV(+1, 0, 0) }};
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = {{ 1 }};
};

// 2D surfaces

/**
 * Left (X-) surface gauss points for 2D box (surface ID -1).
 */
template <int order>
struct BoxItgPts<order, 2, -1> {
  ///! number of integration points
  static constexpr int n_itg_pts = BoxItgPts<order, 1, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(-1, 0, 0) * Matrix4::rotateZ(-90) * BoxItgPts<order, 1, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 1, 0>::weights;
};

/**
 * Right (X+) surface gauss points for 2D box (surface ID +1).
 */
template <int order>
struct BoxItgPts<order, 2, +1> {
  ///! number of integration points
  static constexpr int n_itg_pts = BoxItgPts<order, 1, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(+1, 0, 0) * Matrix4::rotateZ(+90) * BoxItgPts<order, 1, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 1, 0>::weights;
};

/**
 * Bottom (Y-) surface gauss points for 2D box (surface ID -2).
 */
template <int order>
struct BoxItgPts<order, 2, -2> {
  ///! number of integration points
  static constexpr int n_itg_pts = BoxItgPts<order, 1, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, -1, 0) * Matrix4::rotateZ(0) * BoxItgPts<order, 1, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 1, 0>::weights;
};

/**
 * Bottom (Y+) surface gauss points for 2D box (surface ID +2).
 */
template <int order>
struct BoxItgPts<order, 2, +2> {
  ///! number of integration points
  static constexpr int n_itg_pts = BoxItgPts<order, 1, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, +1, 0) * Matrix4::rotateZ(0) * BoxItgPts<order, 1, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 1, 0>::weights;
};

// 3D surfaces

/**
 * Left (X-) surface gauss points for 3D hexahedron (surface ID -1).
 */
template<int order>
struct BoxItgPts<order, 3, -1> {
  ///! number of integration points
  static constexpr int n_itg_pts = BoxItgPts<order, 2, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(-1, 0, 0) * Matrix4::rotateY(-90) * BoxItgPts<order, 2, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 2, 0>::weights;
};

/**
 * Right (X+) surface gauss points for 3D hexahedron (surface ID +1).
 */
template<int order>
struct BoxItgPts<order, 3, +1> {
  ///! number of integration points
  static constexpr int n_itg_pts = BoxItgPts<order, 2, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(+1, 0, 0) * Matrix4::rotateY(90) * BoxItgPts<order, 2, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 2, 0>::weights;
};

/**
 * Bottom (Y-) surface gauss points for 3D hexahedron (surface ID -2).
 */
template<int order>
struct BoxItgPts<order, 3, -2> {
  ///! number of integration points
  static constexpr int n_itg_pts = BoxItgPts<order, 2, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, -1, 0) * Matrix4::rotateX(-90) * BoxItgPts<order, 2, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 2, 0>::weights;
};

/**
 * Top (Y+) surface gauss points for 3D hexahedron (surface ID +2).
 */
template<int order>
struct BoxItgPts<order, 3, +2> {
  ///! number of integration points
  static constexpr int n_itg_pts = BoxItgPts<order, 2, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, +1, 0) * Matrix4::rotateX(90) * BoxItgPts<order, 2, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 2, 0>::weights;
};

/**
 * Back (Z-) surface gauss points for 3D hexahedron (surface ID -3).
 */
template<int order>
struct BoxItgPts<order, 3, -3> {
  ///! number of integration points
  static constexpr int n_itg_pts = BoxItgPts<order, 2, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0, -1) * BoxItgPts<order, 2, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 2, 0>::weights;
};

/**
 * Front (Z+) surface gauss points for 3D hexahedron (surface ID +3).
 */
template<int order>
struct BoxItgPts<order, 3, +3> {
  ///! number of integration points
  static constexpr int n_itg_pts = BoxItgPts<order, 2, 0>::n_itg_pts;
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0, +1) * BoxItgPts<order, 2, 0>::itg_pts;
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 2, 0>::weights;
};
#ifdef ENABLE_4D

// 4D surfaces
// order = 2
/**
 * Left (X-) surface gauss points for 4D tesseract (surface ID -1).
 */
template<>
struct BoxItgPts<2, 4, -1> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 8;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(-1, 0, 0) * Matrix4::rotateY(-90) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 8> itg_pts = {{
                                                              ZEROPTV(-1.0, -0.5773502691896258, -0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(-1.0, 0.5773502691896258, -0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(-1.0, -0.5773502691896258, 0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(-1.0, 0.5773502691896258, 0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(-1.0, -0.5773502691896258, -0.5773502691896258, 0.5773502691896258),
                                                              ZEROPTV(-1.0, 0.5773502691896258, -0.5773502691896258, 0.5773502691896258),
                                                              ZEROPTV(-1.0, -0.5773502691896258, 0.5773502691896258, 0.5773502691896258),
                                                              ZEROPTV(-1.0, 0.5773502691896258, 0.5773502691896258, 0.5773502691896258)
                                                          }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 8> weights = {{
                                                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
                                                         }};
};

/**
 * Right (X+) surface gauss points for 4D tesseract (surface ID +1).
 */
template<>
struct BoxItgPts<2, 4, +1> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 8;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(+1, 0, 0) * Matrix4::rotateY(90) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 8> itg_pts = {{
                                                              ZEROPTV(+1.0, -0.5773502691896258, -0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(+1.0, 0.5773502691896258, -0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(+1.0, -0.5773502691896258, 0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(+1.0, 0.5773502691896258, 0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(+1.0, -0.5773502691896258, -0.5773502691896258, 0.5773502691896258),
                                                              ZEROPTV(+1.0, 0.5773502691896258, -0.5773502691896258, 0.5773502691896258),
                                                              ZEROPTV(+1.0, -0.5773502691896258, 0.5773502691896258, 0.5773502691896258),
                                                              ZEROPTV(+1.0, 0.5773502691896258, 0.5773502691896258, 0.5773502691896258)
                                                          }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 8> weights = {{
                                                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
                                                         }};
};

/**
 * Bottom (Y-) surface gauss points for 4D tesseract (surface ID -2).
 */
template<>
struct BoxItgPts<2, 4, -2> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 8;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, -1, 0) * Matrix4::rotateX(-90) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 8> itg_pts = {{
                                                              ZEROPTV(-0.5773502691896258, -1.0, -0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, -1.0, -0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(-0.5773502691896258, -1.0, 0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, -1.0, 0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(-0.5773502691896258, -1.0, -0.5773502691896258, 0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, -1.0, -0.5773502691896258, 0.5773502691896258),
                                                              ZEROPTV(-0.5773502691896258, -1.0, 0.5773502691896258, 0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, -1.0, 0.5773502691896258, 0.5773502691896258)
                                                          }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 8> weights = {{
                                                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
                                                         }};
};

/**
 * Top (Y+) surface gauss points for 4D tesseract (surface ID +2).
 */
template<>
struct BoxItgPts<2, 4, +2> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 8;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, +1, 0) * Matrix4::rotateX(90) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 8> itg_pts = {{
                                                              ZEROPTV(-0.5773502691896258, +1.0, -0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, +1.0, -0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(-0.5773502691896258, +1.0, 0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, +1.0, 0.5773502691896258, -0.5773502691896258),
                                                              ZEROPTV(-0.5773502691896258, +1.0, -0.5773502691896258, 0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, +1.0, -0.5773502691896258, 0.5773502691896258),
                                                              ZEROPTV(-0.5773502691896258, +1.0, 0.5773502691896258, 0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, +1.0, 0.5773502691896258, 0.5773502691896258)
                                                          }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 8> weights = {{
                                                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
                                                         }};
};

/**
 * Back (Z-) surface gauss points for 4D tesseract (surface ID -3).
 */
template<>
struct BoxItgPts<2, 4, -3> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 8;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0, -1) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 8> itg_pts = {{
                                                              ZEROPTV(-0.5773502691896258, -0.5773502691896258, -1.0, -0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, -0.5773502691896258, -1.0, -0.5773502691896258),
                                                              ZEROPTV(-0.5773502691896258, 0.5773502691896258, -1.0, -0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, 0.5773502691896258, -1.0, -0.5773502691896258),
                                                              ZEROPTV(-0.5773502691896258, -0.5773502691896258, -1.0, 0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, -0.5773502691896258, -1.0, 0.5773502691896258),
                                                              ZEROPTV(-0.5773502691896258, 0.5773502691896258, -1.0, 0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, 0.5773502691896258, -1.0, 0.5773502691896258)
                                                          }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 8> weights = {{
                                                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
                                                         }};
};

/**
 * Front (Z+) surface gauss points for 4D tesseract (surface ID +3).
 */
template<>
struct BoxItgPts<2, 4, +3> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 8;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0, +1) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 8> itg_pts = {{
                                                              ZEROPTV(-0.5773502691896258, -0.5773502691896258, +1.0, -0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, -0.5773502691896258, +1.0, -0.5773502691896258),
                                                              ZEROPTV(-0.5773502691896258, 0.5773502691896258, +1.0, -0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, 0.5773502691896258, +1.0, -0.5773502691896258),
                                                              ZEROPTV(-0.5773502691896258, -0.5773502691896258, +1.0, 0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, -0.5773502691896258, +1.0, 0.5773502691896258),
                                                              ZEROPTV(-0.5773502691896258, 0.5773502691896258, +1.0, 0.5773502691896258),
                                                              ZEROPTV(0.5773502691896258, 0.5773502691896258, +1.0, 0.5773502691896258)
                                                          }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 8> weights = {{
                                                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
                                                         }};
};

/**
 * Back (T-) surface gauss points for 4D tesseract (surface ID -4).
 */
template<>
struct BoxItgPts<2, 4, -4> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 8;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0, -1) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 8> itg_pts = {{
                                                              ZEROPTV(-0.5773502691896258, -0.5773502691896258, -0.5773502691896258, -1.0),
                                                              ZEROPTV(0.5773502691896258, -0.5773502691896258, -0.5773502691896258, -1.0),
                                                              ZEROPTV(-0.5773502691896258, 0.5773502691896258, -0.5773502691896258, -1.0),
                                                              ZEROPTV(0.5773502691896258, 0.5773502691896258, -0.5773502691896258, -1.0),
                                                              ZEROPTV(-0.5773502691896258, -0.5773502691896258, 0.5773502691896258, -1.0),
                                                              ZEROPTV(0.5773502691896258, -0.5773502691896258, 0.5773502691896258, -1.0),
                                                              ZEROPTV(-0.5773502691896258, 0.5773502691896258, 0.5773502691896258, -1.0),
                                                              ZEROPTV(0.5773502691896258, 0.5773502691896258, 0.5773502691896258, -1.0)
                                                          }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 8> weights = {{
                                                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
                                                         }};
};

/**
 * Front (T+) surface gauss points for 4D tesseract (surface ID +4).
 */
template<>
struct BoxItgPts<2, 4, +4> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 8;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0, +1) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 8> itg_pts = {{
                                                              ZEROPTV(-0.5773502691896258, -0.5773502691896258, -0.5773502691896258, +1.0),
                                                              ZEROPTV(0.5773502691896258, -0.5773502691896258, -0.5773502691896258, +1.0),
                                                              ZEROPTV(-0.5773502691896258, 0.5773502691896258, -0.5773502691896258, +1.0),
                                                              ZEROPTV(0.5773502691896258, 0.5773502691896258, -0.5773502691896258, +1.0),
                                                              ZEROPTV(-0.5773502691896258, -0.5773502691896258, 0.5773502691896258, +1.0),
                                                              ZEROPTV(0.5773502691896258, -0.5773502691896258, 0.5773502691896258, +1.0),
                                                              ZEROPTV(-0.5773502691896258, 0.5773502691896258, 0.5773502691896258, +1.0),
                                                              ZEROPTV(0.5773502691896258, 0.5773502691896258, 0.5773502691896258, +1.0)
                                                          }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 8> weights = {{
                                                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
                                                         }};
};

// order = 3
/**
 * Left (X-) surface gauss points for 4D tesseract (surface ID -1).
 */
template<>
struct BoxItgPts<3, 4, -1> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 27;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(-1, 0, 0) * Matrix4::rotateY(-90) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 27> itg_pts = {{
                                                               ZEROPTV(-1.0, -0.7745966692414834, -0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(-1.0, 0, -0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(-1.0, 0.7745966692414834, -0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(-1.0, -0.7745966692414834, 0, -0.7745966692414834),
                                                               ZEROPTV(-1.0, 0, 0, -0.7745966692414834),
                                                               ZEROPTV(-1.0, 0.7745966692414834, 0, -0.7745966692414834),
                                                               ZEROPTV(-1.0, -0.7745966692414834, 0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(-1.0, 0, 0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(-1.0, 0.7745966692414834, 0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(-1.0, -0.7745966692414834, -0.7745966692414834, 0),
                                                               ZEROPTV(-1.0, 0, -0.7745966692414834, 0),
                                                               ZEROPTV(-1.0, 0.7745966692414834, -0.7745966692414834, 0),
                                                               ZEROPTV(-1.0, -0.7745966692414834, 0, 0),
                                                               ZEROPTV(-1.0, 0, 0, 0),
                                                               ZEROPTV(-1.0, 0.7745966692414834, 0, 0),
                                                               ZEROPTV(-1.0, -0.7745966692414834, 0.7745966692414834, 0),
                                                               ZEROPTV(-1.0, 0, 0.7745966692414834, 0),
                                                               ZEROPTV(-1.0, 0.7745966692414834, 0.7745966692414834, 0),
                                                               ZEROPTV(-1.0, -0.7745966692414834, -0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(-1.0, 0, -0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(-1.0, 0.7745966692414834, -0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(-1.0, -0.7745966692414834, 0, 0.7745966692414834),
                                                               ZEROPTV(-1.0, 0, 0, 0.7745966692414834),
                                                               ZEROPTV(-1.0, 0.7745966692414834, 0, 0.7745966692414834),
                                                               ZEROPTV(-1.0, -0.7745966692414834, 0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(-1.0, 0, 0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(-1.0, 0.7745966692414834, 0.7745966692414834, 0.7745966692414834)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 27> weights = {{
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.7023319615912208,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666
                                                          }};
};

/**
 * Right (X+) surface gauss points for 4D tesseract (surface ID +1).
 */
template<>
struct BoxItgPts<3, 4, +1> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 27;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(+1, 0, 0) * Matrix4::rotateY(90) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 27> itg_pts = {{
                                                               ZEROPTV(+1.0, -0.7745966692414834, -0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(+1.0, 0, -0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(+1.0, 0.7745966692414834, -0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(+1.0, -0.7745966692414834, 0, -0.7745966692414834),
                                                               ZEROPTV(+1.0, 0, 0, -0.7745966692414834),
                                                               ZEROPTV(+1.0, 0.7745966692414834, 0, -0.7745966692414834),
                                                               ZEROPTV(+1.0, -0.7745966692414834, 0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(+1.0, 0, 0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(+1.0, 0.7745966692414834, 0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(+1.0, -0.7745966692414834, -0.7745966692414834, 0),
                                                               ZEROPTV(+1.0, 0, -0.7745966692414834, 0),
                                                               ZEROPTV(+1.0, 0.7745966692414834, -0.7745966692414834, 0),
                                                               ZEROPTV(+1.0, -0.7745966692414834, 0, 0),
                                                               ZEROPTV(+1.0, 0, 0, 0),
                                                               ZEROPTV(+1.0, 0.7745966692414834, 0, 0),
                                                               ZEROPTV(+1.0, -0.7745966692414834, 0.7745966692414834, 0),
                                                               ZEROPTV(+1.0, 0, 0.7745966692414834, 0),
                                                               ZEROPTV(+1.0, 0.7745966692414834, 0.7745966692414834, 0),
                                                               ZEROPTV(+1.0, -0.7745966692414834, -0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(+1.0, 0, -0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(+1.0, 0.7745966692414834, -0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(+1.0, -0.7745966692414834, 0, 0.7745966692414834),
                                                               ZEROPTV(+1.0, 0, 0, 0.7745966692414834),
                                                               ZEROPTV(+1.0, 0.7745966692414834, 0, 0.7745966692414834),
                                                               ZEROPTV(+1.0, -0.7745966692414834, 0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(+1.0, 0, 0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(+1.0, 0.7745966692414834, 0.7745966692414834, 0.7745966692414834)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 27> weights = {{
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.7023319615912208,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666
                                                          }};
};

/**
 * Bottom (Y-) surface gauss points for 4D tesseract (surface ID -2).
 */
template<>
struct BoxItgPts<3, 4, -2> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 27;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, -1, 0) * Matrix4::rotateX(-90) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 27> itg_pts = {{
                                                               ZEROPTV(-0.7745966692414834, -1.0, -0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(0, -1.0, -0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, -1.0, -0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, -1.0, 0, -0.7745966692414834),
                                                               ZEROPTV(0, -1.0, 0, -0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, -1.0, 0, -0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, -1.0, 0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(0, -1.0, 0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, -1.0, 0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, -1.0, -0.7745966692414834, 0),
                                                               ZEROPTV(0, -1.0, -0.7745966692414834, 0),
                                                               ZEROPTV(0.7745966692414834, -1.0, -0.7745966692414834, 0),
                                                               ZEROPTV(-0.7745966692414834, -1.0, 0, 0),
                                                               ZEROPTV(0, -1.0, 0, 0),
                                                               ZEROPTV(0.7745966692414834, -1.0, 0, 0),
                                                               ZEROPTV(-0.7745966692414834, -1.0, 0.7745966692414834, 0),
                                                               ZEROPTV(0, -1.0, 0.7745966692414834, 0),
                                                               ZEROPTV(0.7745966692414834, -1.0, 0.7745966692414834, 0),
                                                               ZEROPTV(-0.7745966692414834, -1.0, -0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(0, -1.0, -0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, -1.0, -0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, -1.0, 0, 0.7745966692414834),
                                                               ZEROPTV(0, -1.0, 0, 0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, -1.0, 0, 0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, -1.0, 0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(0, -1.0, 0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, -1.0, 0.7745966692414834, 0.7745966692414834)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 27> weights = {{
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.7023319615912208,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666
                                                          }};
};

/**
 * Top (Y+) surface gauss points for 4D tesseract (surface ID +2).
 */
template<>
struct BoxItgPts<3, 4, +2> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 27;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, +1, 0) * Matrix4::rotateX(90) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 27> itg_pts = {{
                                                               ZEROPTV(-0.7745966692414834, +1.0, -0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(0, +1.0, -0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, +1.0, -0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, +1.0, 0, -0.7745966692414834),
                                                               ZEROPTV(0, +1.0, 0, -0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, +1.0, 0, -0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, +1.0, 0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(0, +1.0, 0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, +1.0, 0.7745966692414834, -0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, +1.0, -0.7745966692414834, 0),
                                                               ZEROPTV(0, +1.0, -0.7745966692414834, 0),
                                                               ZEROPTV(0.7745966692414834, +1.0, -0.7745966692414834, 0),
                                                               ZEROPTV(-0.7745966692414834, +1.0, 0, 0),
                                                               ZEROPTV(0, +1.0, 0, 0),
                                                               ZEROPTV(0.7745966692414834, +1.0, 0, 0),
                                                               ZEROPTV(-0.7745966692414834, +1.0, 0.7745966692414834, 0),
                                                               ZEROPTV(0, +1.0, 0.7745966692414834, 0),
                                                               ZEROPTV(0.7745966692414834, +1.0, 0.7745966692414834, 0),
                                                               ZEROPTV(-0.7745966692414834, +1.0, -0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(0, +1.0, -0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, +1.0, -0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, +1.0, 0, 0.7745966692414834),
                                                               ZEROPTV(0, +1.0, 0, 0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, +1.0, 0, 0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, +1.0, 0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(0, +1.0, 0.7745966692414834, 0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, +1.0, 0.7745966692414834, 0.7745966692414834)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 27> weights = {{
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.7023319615912208,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666
                                                          }};
};

/**
 * Back (Z-) surface gauss points for 4D tesseract (surface ID -3).
 */
template<>
struct BoxItgPts<3, 4, -3> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 27;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0, -1) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 27> itg_pts = {{
                                                               ZEROPTV(-0.7745966692414834, -0.7745966692414834, -1.0, -0.7745966692414834),
                                                               ZEROPTV(0, -0.7745966692414834, -1.0, -0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, -0.7745966692414834, -1.0, -0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, 0, -1.0, -0.7745966692414834),
                                                               ZEROPTV(0, 0, -1.0, -0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, 0, -1.0, -0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, 0.7745966692414834, -1.0, -0.7745966692414834),
                                                               ZEROPTV(0, 0.7745966692414834, -1.0, -0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, 0.7745966692414834, -1.0, -0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, -0.7745966692414834, -1.0, 0),
                                                               ZEROPTV(0, -0.7745966692414834, -1.0, 0),
                                                               ZEROPTV(0.7745966692414834, -0.7745966692414834, -1.0, 0),
                                                               ZEROPTV(-0.7745966692414834, 0, -1.0, 0),
                                                               ZEROPTV(0, 0, -1.0, 0),
                                                               ZEROPTV(0.7745966692414834, 0, -1.0, 0),
                                                               ZEROPTV(-0.7745966692414834, 0.7745966692414834, -1.0, 0),
                                                               ZEROPTV(0, 0.7745966692414834, -1.0, 0),
                                                               ZEROPTV(0.7745966692414834, 0.7745966692414834, -1.0, 0),
                                                               ZEROPTV(-0.7745966692414834, -0.7745966692414834, -1.0, 0.7745966692414834),
                                                               ZEROPTV(0, -0.7745966692414834, -1.0, 0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, -0.7745966692414834, -1.0, 0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, 0, -1.0, 0.7745966692414834),
                                                               ZEROPTV(0, 0, -1.0, 0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, 0, -1.0, 0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, 0.7745966692414834, -1.0, 0.7745966692414834),
                                                               ZEROPTV(0, 0.7745966692414834, -1.0, 0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, 0.7745966692414834, -1.0, 0.7745966692414834)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 27> weights = {{
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.7023319615912208,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666
                                                          }};
};

/**
 * Front (Z+) surface gauss points for 4D tesseract (surface ID +3).
 */
template<>
struct BoxItgPts<3, 4, +3> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 27;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0, +1) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 27> itg_pts = {{
                                                               ZEROPTV(-0.7745966692414834, -0.7745966692414834, +1.0, -0.7745966692414834),
                                                               ZEROPTV(0, -0.7745966692414834, +1.0, -0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, -0.7745966692414834, +1.0, -0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, 0, +1.0, -0.7745966692414834),
                                                               ZEROPTV(0, 0, +1.0, -0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, 0, +1.0, -0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, 0.7745966692414834, +1.0, -0.7745966692414834),
                                                               ZEROPTV(0, 0.7745966692414834, +1.0, -0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, 0.7745966692414834, +1.0, -0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, -0.7745966692414834, +1.0, 0),
                                                               ZEROPTV(0, -0.7745966692414834, +1.0, 0),
                                                               ZEROPTV(0.7745966692414834, -0.7745966692414834, +1.0, 0),
                                                               ZEROPTV(-0.7745966692414834, 0, +1.0, 0),
                                                               ZEROPTV(0, 0, +1.0, 0),
                                                               ZEROPTV(0.7745966692414834, 0, +1.0, 0),
                                                               ZEROPTV(-0.7745966692414834, 0.7745966692414834, +1.0, 0),
                                                               ZEROPTV(0, 0.7745966692414834, +1.0, 0),
                                                               ZEROPTV(0.7745966692414834, 0.7745966692414834, +1.0, 0),
                                                               ZEROPTV(-0.7745966692414834, -0.7745966692414834, +1.0, 0.7745966692414834),
                                                               ZEROPTV(0, -0.7745966692414834, +1.0, 0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, -0.7745966692414834, +1.0, 0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, 0, +1.0, 0.7745966692414834),
                                                               ZEROPTV(0, 0, +1.0, 0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, 0, +1.0, 0.7745966692414834),
                                                               ZEROPTV(-0.7745966692414834, 0.7745966692414834, +1.0, 0.7745966692414834),
                                                               ZEROPTV(0, 0.7745966692414834, +1.0, 0.7745966692414834),
                                                               ZEROPTV(0.7745966692414834, 0.7745966692414834, +1.0, 0.7745966692414834)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 27> weights = {{
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.7023319615912208,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666
                                                          }};
};

/**
 * Back (T-) surface gauss points for 4D tesseract (surface ID -4).
 */
template<>
struct BoxItgPts<3, 4, -4> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 27;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0, -1) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 27> itg_pts = {{
                                                               ZEROPTV(-0.7745966692414834, -0.7745966692414834, -0.7745966692414834, -1.0),
                                                               ZEROPTV(0, -0.7745966692414834, -0.7745966692414834, -1.0),
                                                               ZEROPTV(0.7745966692414834, -0.7745966692414834, -0.7745966692414834, -1.0),
                                                               ZEROPTV(-0.7745966692414834, 0, -0.7745966692414834, -1.0),
                                                               ZEROPTV(0, 0, -0.7745966692414834, -1.0),
                                                               ZEROPTV(0.7745966692414834, 0, -0.7745966692414834, -1.0),
                                                               ZEROPTV(-0.7745966692414834, 0.7745966692414834, -0.7745966692414834, -1.0),
                                                               ZEROPTV(0, 0.7745966692414834, -0.7745966692414834, -1.0),
                                                               ZEROPTV(0.7745966692414834, 0.7745966692414834, -0.7745966692414834, -1.0),
                                                               ZEROPTV(-0.7745966692414834, -0.7745966692414834, 0, -1.0),
                                                               ZEROPTV(0, -0.7745966692414834, 0, -1.0),
                                                               ZEROPTV(0.7745966692414834, -0.7745966692414834, 0, -1.0),
                                                               ZEROPTV(-0.7745966692414834, 0, 0, -1.0),
                                                               ZEROPTV(0, 0, 0, -1.0),
                                                               ZEROPTV(0.7745966692414834, 0, 0, -1.0),
                                                               ZEROPTV(-0.7745966692414834, 0.7745966692414834, 0, -1.0),
                                                               ZEROPTV(0, 0.7745966692414834, 0, -1.0),
                                                               ZEROPTV(0.7745966692414834, 0.7745966692414834, 0, -1.0),
                                                               ZEROPTV(-0.7745966692414834, -0.7745966692414834, 0.7745966692414834, -1.0),
                                                               ZEROPTV(0, -0.7745966692414834, 0.7745966692414834, -1.0),
                                                               ZEROPTV(0.7745966692414834, -0.7745966692414834, 0.7745966692414834, -1.0),
                                                               ZEROPTV(-0.7745966692414834, 0, 0.7745966692414834, -1.0),
                                                               ZEROPTV(0, 0, 0.7745966692414834, -1.0),
                                                               ZEROPTV(0.7745966692414834, 0, 0.7745966692414834, -1.0),
                                                               ZEROPTV(-0.7745966692414834, 0.7745966692414834, 0.7745966692414834, -1.0),
                                                               ZEROPTV(0, 0.7745966692414834, 0.7745966692414834, -1.0),
                                                               ZEROPTV(0.7745966692414834, 0.7745966692414834, 0.7745966692414834, -1.0)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 27> weights = {{
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.7023319615912208,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666
                                                          }};
};

/**
 * Front (T+) surface gauss points for 4D tesseract (surface ID +4).
 */
template<>
struct BoxItgPts<3, 4, +4> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 27;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0, +1) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 27> itg_pts = {{
                                                               ZEROPTV(-0.7745966692414834, -0.7745966692414834, -0.7745966692414834, +1.0),
                                                               ZEROPTV(0, -0.7745966692414834, -0.7745966692414834, +1.0),
                                                               ZEROPTV(0.7745966692414834, -0.7745966692414834, -0.7745966692414834, +1.0),
                                                               ZEROPTV(-0.7745966692414834, 0, -0.7745966692414834, +1.0),
                                                               ZEROPTV(0, 0, -0.7745966692414834, +1.0),
                                                               ZEROPTV(0.7745966692414834, 0, -0.7745966692414834, +1.0),
                                                               ZEROPTV(-0.7745966692414834, 0.7745966692414834, -0.7745966692414834, +1.0),
                                                               ZEROPTV(0, 0.7745966692414834, -0.7745966692414834, +1.0),
                                                               ZEROPTV(0.7745966692414834, 0.7745966692414834, -0.7745966692414834, +1.0),
                                                               ZEROPTV(-0.7745966692414834, -0.7745966692414834, 0, +1.0),
                                                               ZEROPTV(0, -0.7745966692414834, 0, +1.0),
                                                               ZEROPTV(0.7745966692414834, -0.7745966692414834, 0, +1.0),
                                                               ZEROPTV(-0.7745966692414834, 0, 0, +1.0),
                                                               ZEROPTV(0, 0, 0, +1.0),
                                                               ZEROPTV(0.7745966692414834, 0, 0, +1.0),
                                                               ZEROPTV(-0.7745966692414834, 0.7745966692414834, 0, +1.0),
                                                               ZEROPTV(0, 0.7745966692414834, 0, +1.0),
                                                               ZEROPTV(0.7745966692414834, 0.7745966692414834, 0, +1.0),
                                                               ZEROPTV(-0.7745966692414834, -0.7745966692414834, 0.7745966692414834, +1.0),
                                                               ZEROPTV(0, -0.7745966692414834, 0.7745966692414834, +1.0),
                                                               ZEROPTV(0.7745966692414834, -0.7745966692414834, 0.7745966692414834, +1.0),
                                                               ZEROPTV(-0.7745966692414834, 0, 0.7745966692414834, +1.0),
                                                               ZEROPTV(0, 0, 0.7745966692414834, +1.0),
                                                               ZEROPTV(0.7745966692414834, 0, 0.7745966692414834, +1.0),
                                                               ZEROPTV(-0.7745966692414834, 0.7745966692414834, 0.7745966692414834, +1.0),
                                                               ZEROPTV(0, 0.7745966692414834, 0.7745966692414834, +1.0),
                                                               ZEROPTV(0.7745966692414834, 0.7745966692414834, 0.7745966692414834, +1.0)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 27> weights = {{
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.7023319615912208,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.438957475994513,
                                                              0.27434842249657065,
                                                              0.17146776406035666,
                                                              0.27434842249657065,
                                                              0.17146776406035666
                                                          }};
};

// order = 4
/**
 * Left (X-) surface gauss points for 4D tesseract (surface ID -1).
 */
template<>
struct BoxItgPts<4, 4, -1> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 64;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(-1, 0, 0) * Matrix4::rotateY(-90) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 64> itg_pts = {{
                                                               ZEROPTV(-1.0, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.3399810435848562, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.3399810435848562, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.8611363115940526, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.8611363115940526, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.3399810435848562, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.3399810435848562, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.8611363115940526, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.8611363115940526, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.3399810435848562, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.3399810435848562, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.8611363115940526, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.8611363115940526, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.3399810435848562, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.3399810435848562, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.8611363115940526, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.8611363115940526, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.3399810435848562, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.3399810435848562, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.8611363115940526, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.8611363115940526, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.3399810435848562, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.3399810435848562, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.8611363115940526, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.8611363115940526, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.3399810435848562, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.3399810435848562, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.8611363115940526, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.8611363115940526, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.3399810435848562, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.3399810435848562, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.8611363115940526, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.8611363115940526, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.3399810435848562, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.3399810435848562, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.8611363115940526, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.8611363115940526, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.3399810435848562, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.3399810435848562, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.8611363115940526, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.8611363115940526, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.3399810435848562, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.3399810435848562, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.8611363115940526, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.8611363115940526, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.3399810435848562, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.3399810435848562, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-1.0, 0.8611363115940526, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-1.0, -0.8611363115940526, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.3399810435848562, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.3399810435848562, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.8611363115940526, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.8611363115940526, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.3399810435848562, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.3399810435848562, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.8611363115940526, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.8611363115940526, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.3399810435848562, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.3399810435848562, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.8611363115940526, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.8611363115940526, 0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(-1.0, -0.3399810435848562, 0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.3399810435848562, 0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(-1.0, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 64> weights = {{
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147
                                                          }};
};

/**
 * Right (X+) surface gauss points for 4D tesseract (surface ID +1).
 */
template<>
struct BoxItgPts<4, 4, +1> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 64;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(+1, 0, 0) * Matrix4::rotateY(90) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 64> itg_pts = {{
                                                               ZEROPTV(+1.0, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.3399810435848562, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.3399810435848562, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.8611363115940526, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.8611363115940526, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.3399810435848562, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.3399810435848562, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.8611363115940526, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.8611363115940526, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.3399810435848562, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.3399810435848562, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.8611363115940526, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.8611363115940526, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.3399810435848562, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.3399810435848562, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.8611363115940526, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.8611363115940526, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.3399810435848562, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.3399810435848562, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.8611363115940526, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.8611363115940526, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.3399810435848562, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.3399810435848562, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.8611363115940526, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.8611363115940526, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.3399810435848562, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.3399810435848562, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.8611363115940526, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.8611363115940526, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.3399810435848562, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.3399810435848562, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.8611363115940526, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.8611363115940526, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.3399810435848562, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.3399810435848562, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.8611363115940526, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.8611363115940526, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.3399810435848562, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.3399810435848562, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.8611363115940526, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.8611363115940526, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.3399810435848562, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.3399810435848562, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.8611363115940526, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.8611363115940526, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.3399810435848562, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.3399810435848562, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(+1.0, 0.8611363115940526, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(+1.0, -0.8611363115940526, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.3399810435848562, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.3399810435848562, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.8611363115940526, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.8611363115940526, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.3399810435848562, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.3399810435848562, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.8611363115940526, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.8611363115940526, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.3399810435848562, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.3399810435848562, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.8611363115940526, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.8611363115940526, 0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(+1.0, -0.3399810435848562, 0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.3399810435848562, 0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(+1.0, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 64> weights = {{
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147
                                                          }};
};

/**
 * Bottom (Y-) surface gauss points for 4D tesseract (surface ID -2).
 */
template<>
struct BoxItgPts<4, 4, -2> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 64;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, -1, 0) * Matrix4::rotateX(-90) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 64> itg_pts = {{
                                                               ZEROPTV(-0.8611363115940526, -1.0, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -1.0, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -1.0, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -1.0, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, -1.0, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -1.0, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -1.0, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -1.0, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, -1.0, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -1.0, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -1.0, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -1.0, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, -1.0, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -1.0, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -1.0, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -1.0, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, -1.0, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -1.0, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -1.0, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -1.0, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -1.0, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -1.0, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -1.0, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -1.0, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -1.0, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -1.0, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -1.0, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -1.0, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -1.0, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -1.0, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -1.0, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -1.0, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -1.0, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -1.0, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -1.0, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -1.0, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -1.0, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -1.0, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -1.0, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -1.0, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -1.0, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -1.0, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -1.0, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -1.0, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -1.0, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -1.0, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -1.0, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -1.0, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -1.0, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -1.0, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -1.0, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -1.0, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, -1.0, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -1.0, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -1.0, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -1.0, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, -1.0, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -1.0, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -1.0, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -1.0, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, -1.0, 0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -1.0, 0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -1.0, 0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -1.0, 0.8611363115940526, 0.8611363115940526)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 64> weights = {{
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147
                                                          }};
};

/**
 * Top (Y+) surface gauss points for 4D tesseract (surface ID +2).
 */
template<>
struct BoxItgPts<4, 4, +2> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 64;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, +1, 0) * Matrix4::rotateX(90) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 64> itg_pts = {{
                                                               ZEROPTV(-0.8611363115940526, +1.0, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, +1.0, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, +1.0, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, +1.0, -0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, +1.0, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, +1.0, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, +1.0, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, +1.0, -0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, +1.0, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, +1.0, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, +1.0, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, +1.0, 0.3399810435848562, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, +1.0, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, +1.0, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, +1.0, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, +1.0, 0.8611363115940526, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, +1.0, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, +1.0, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, +1.0, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, +1.0, -0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, +1.0, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, +1.0, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, +1.0, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, +1.0, -0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, +1.0, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, +1.0, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, +1.0, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, +1.0, 0.3399810435848562, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, +1.0, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, +1.0, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, +1.0, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, +1.0, 0.8611363115940526, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, +1.0, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, +1.0, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, +1.0, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, +1.0, -0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, +1.0, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, +1.0, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, +1.0, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, +1.0, -0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, +1.0, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, +1.0, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, +1.0, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, +1.0, 0.3399810435848562, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, +1.0, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, +1.0, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, +1.0, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, +1.0, 0.8611363115940526, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, +1.0, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, +1.0, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, +1.0, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, +1.0, -0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, +1.0, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, +1.0, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, +1.0, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, +1.0, -0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, +1.0, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, +1.0, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, +1.0, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, +1.0, 0.3399810435848562, 0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, +1.0, 0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, +1.0, 0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, +1.0, 0.8611363115940526, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, +1.0, 0.8611363115940526, 0.8611363115940526)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 64> weights = {{
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147
                                                          }};
};

/**
 * Back (Z-) surface gauss points for 4D tesseract (surface ID -3).
 */
template<>
struct BoxItgPts<4, 4, -3> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 64;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0, -1) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 64> itg_pts = {{
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, -1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, -1.0, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, -1.0, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, -1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, -1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, -1.0, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, -1.0, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, -1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, -1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, -1.0, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, -1.0, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, -1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, -1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, -1.0, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, -1.0, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, -1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, -1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, -1.0, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, -1.0, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, -1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, -1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, -1.0, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, -1.0, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, -1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, -1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, -1.0, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, -1.0, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, -1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, -1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, -1.0, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, -1.0, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, -1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, -1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, -1.0, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, -1.0, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, -1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, -1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, -1.0, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, -1.0, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, -1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, -1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, -1.0, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, -1.0, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, -1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, -1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, -1.0, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, -1.0, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, -1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, -1.0, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, -1.0, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, -1.0, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, -1.0, 0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, -1.0, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, -1.0, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, -1.0, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, -1.0, 0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, -1.0, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, -1.0, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, -1.0, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, -1.0, 0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, -1.0, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, -1.0, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, -1.0, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, -1.0, 0.8611363115940526)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 64> weights = {{
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147
                                                          }};
};

/**
 * Front (Z+) surface gauss points for 4D tesseract (surface ID +3).
 */
template<>
struct BoxItgPts<4, 4, +3> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 64;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0, +1) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 64> itg_pts = {{
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, +1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, +1.0, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, +1.0, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, +1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, +1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, +1.0, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, +1.0, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, +1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, +1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, +1.0, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, +1.0, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, +1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, +1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, +1.0, -0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, +1.0, -0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, +1.0, -0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, +1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, +1.0, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, +1.0, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, +1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, +1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, +1.0, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, +1.0, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, +1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, +1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, +1.0, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, +1.0, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, +1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, +1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, +1.0, -0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, +1.0, -0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, +1.0, -0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, +1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, +1.0, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, +1.0, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, +1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, +1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, +1.0, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, +1.0, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, +1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, +1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, +1.0, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, +1.0, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, +1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, +1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, +1.0, 0.3399810435848562),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, +1.0, 0.3399810435848562),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, +1.0, 0.3399810435848562),
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, +1.0, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, +1.0, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, +1.0, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, +1.0, 0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, +1.0, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, +1.0, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, +1.0, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, +1.0, 0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, +1.0, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, +1.0, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, +1.0, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, +1.0, 0.8611363115940526),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, +1.0, 0.8611363115940526),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, +1.0, 0.8611363115940526),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, +1.0, 0.8611363115940526),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, +1.0, 0.8611363115940526)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 64> weights = {{
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147
                                                          }};
};

/**
 * Back (T-) surface gauss points for 4D tesseract (surface ID -4).
 */
template<>
struct BoxItgPts<4, 4, -4> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 64;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0, -1) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 64> itg_pts = {{
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, -0.8611363115940526, -1.0),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, -0.8611363115940526, -1.0),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, -0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, -0.8611363115940526, -1.0),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, -0.8611363115940526, -1.0),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, -0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, -0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, -0.8611363115940526, -1.0),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, -0.8611363115940526, -1.0),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, -0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, -0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, -0.8611363115940526, -1.0),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, -0.8611363115940526, -1.0),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, -0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, -0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, -0.3399810435848562, -1.0),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, -0.3399810435848562, -1.0),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, -0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, -0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, -0.3399810435848562, -1.0),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, -0.3399810435848562, -1.0),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, -0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, -0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, -0.3399810435848562, -1.0),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, -0.3399810435848562, -1.0),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, -0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, -0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, -0.3399810435848562, -1.0),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, -0.3399810435848562, -1.0),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, -0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, 0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, 0.3399810435848562, -1.0),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, 0.3399810435848562, -1.0),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, 0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, 0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, 0.3399810435848562, -1.0),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, 0.3399810435848562, -1.0),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, 0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, 0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, 0.3399810435848562, -1.0),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, 0.3399810435848562, -1.0),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, 0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, 0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, 0.3399810435848562, -1.0),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, 0.3399810435848562, -1.0),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, 0.3399810435848562, -1.0),
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, 0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, 0.8611363115940526, -1.0),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, 0.8611363115940526, -1.0),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, 0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, 0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, 0.8611363115940526, -1.0),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, 0.8611363115940526, -1.0),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, 0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, 0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, 0.8611363115940526, -1.0),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, 0.8611363115940526, -1.0),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, 0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, 0.8611363115940526, -1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, 0.8611363115940526, -1.0),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, 0.8611363115940526, -1.0),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, 0.8611363115940526, -1.0)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 64> weights = {{
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147
                                                          }};
};

/**
 * Front (T+) surface gauss points for 4D tesseract (surface ID +4).
 */
template<>
struct BoxItgPts<4, 4, +4> {
  ///! number of integration points
  // static constexpr int n_itg_pts = BoxItgPts<order, 3, 0>::n_itg_pts;
  static constexpr int n_itg_pts = 64;
  ///! integration points
  // static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = Matrix4::translate(0, 0, +1) * BoxItgPts<order, 3, 0>::itg_pts;
  static constexpr constexpr_array<ZEROPTV, 64> itg_pts = {{
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, -0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, -0.8611363115940526, +1.0),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, -0.8611363115940526, +1.0),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, -0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, -0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, -0.8611363115940526, +1.0),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, -0.8611363115940526, +1.0),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, -0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, -0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, -0.8611363115940526, +1.0),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, -0.8611363115940526, +1.0),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, -0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, -0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, -0.8611363115940526, +1.0),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, -0.8611363115940526, +1.0),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, -0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, -0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, -0.3399810435848562, +1.0),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, -0.3399810435848562, +1.0),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, -0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, -0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, -0.3399810435848562, +1.0),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, -0.3399810435848562, +1.0),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, -0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, -0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, -0.3399810435848562, +1.0),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, -0.3399810435848562, +1.0),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, -0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, -0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, -0.3399810435848562, +1.0),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, -0.3399810435848562, +1.0),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, -0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, 0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, 0.3399810435848562, +1.0),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, 0.3399810435848562, +1.0),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, 0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, 0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, 0.3399810435848562, +1.0),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, 0.3399810435848562, +1.0),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, 0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, 0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, 0.3399810435848562, +1.0),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, 0.3399810435848562, +1.0),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, 0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, 0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, 0.3399810435848562, +1.0),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, 0.3399810435848562, +1.0),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, 0.3399810435848562, +1.0),
                                                               ZEROPTV(-0.8611363115940526, -0.8611363115940526, 0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.8611363115940526, 0.8611363115940526, +1.0),
                                                               ZEROPTV(0.3399810435848562, -0.8611363115940526, 0.8611363115940526, +1.0),
                                                               ZEROPTV(0.8611363115940526, -0.8611363115940526, 0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.8611363115940526, -0.3399810435848562, 0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.3399810435848562, -0.3399810435848562, 0.8611363115940526, +1.0),
                                                               ZEROPTV(0.3399810435848562, -0.3399810435848562, 0.8611363115940526, +1.0),
                                                               ZEROPTV(0.8611363115940526, -0.3399810435848562, 0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.3399810435848562, 0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.3399810435848562, 0.8611363115940526, +1.0),
                                                               ZEROPTV(0.3399810435848562, 0.3399810435848562, 0.8611363115940526, +1.0),
                                                               ZEROPTV(0.8611363115940526, 0.3399810435848562, 0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.8611363115940526, 0.8611363115940526, 0.8611363115940526, +1.0),
                                                               ZEROPTV(-0.3399810435848562, 0.8611363115940526, 0.8611363115940526, +1.0),
                                                               ZEROPTV(0.3399810435848562, 0.8611363115940526, 0.8611363115940526, +1.0),
                                                               ZEROPTV(0.8611363115940526, 0.8611363115940526, 0.8611363115940526, +1.0)
                                                           }};
  ///! integration point weights
  // static constexpr constexpr_array<double, n_itg_pts> weights = BoxItgPts<order, 3, 0>::weights;
  static constexpr constexpr_array<double, 64> weights = {{
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.2773529669539129,
                                                              0.2773529669539129,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.1479403360567813,
                                                              0.1479403360567813,
                                                              0.07891151579507055,
                                                              0.04209147749053147,
                                                              0.07891151579507055,
                                                              0.07891151579507055,
                                                              0.04209147749053147
                                                          }};
};

#endif
#define DEF_BOX_SURF(nsd, surf_id) \
template<int order> \
constexpr constexpr_array<ZEROPTV, BoxItgPts<order, nsd, surf_id>::n_itg_pts> BoxItgPts<order, nsd, surf_id>::itg_pts; \
template<int order> \
constexpr constexpr_array<double, BoxItgPts<order, nsd, surf_id>::n_itg_pts> BoxItgPts<order, nsd, surf_id>::weights;

DEF_BOX_SURF(1, -1)
DEF_BOX_SURF(1, +1)

DEF_BOX_SURF(2, -1)
DEF_BOX_SURF(2, +1)
DEF_BOX_SURF(2, -2)
DEF_BOX_SURF(2, +2)

DEF_BOX_SURF(3, -1)
DEF_BOX_SURF(3, +1)
DEF_BOX_SURF(3, -2)
DEF_BOX_SURF(3, +2)
DEF_BOX_SURF(3, -3)
DEF_BOX_SURF(3, +3)

}  // namespace TALYFEMLIB
