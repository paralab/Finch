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

#include <limits>

#include <talyfem/data_structures/constexpr_array.h>
#include <talyfem/grid/zeroptv.h>

namespace TALYFEMLIB {

/**
 * A compile-time 4x4 double matrix.
 * Used for calculating gauss points that can be described as transformations
 * of existing gauss points (typically surfaces).
 * Heavily dependant on constexpr.
 */
struct Matrix4 {
  static const int n_rows = 4;  ///< number of rows
  static const int n_cols = 4;  ///< number of columns
  typedef double T;  ///< type of storage

  T data[n_rows*n_cols];  ///< data  // NOLINT(runtime/arrays)

  /**
   * Get an entry in the matrix.
   * @param row row
   * @param col column
   * @returns the value of (row, col)
   */
  constexpr T get(int row, int col = 0) const {
    return data[row * n_cols + col];
  }

  // w coordinate assumed to be 1
  /**
   * Multiply a point by the matrix.
   * @param vec vector/point to multiply
   * @returns the transformed point
   */
  constexpr ZEROPTV operator*(const ZEROPTV& vec) const {
    return ZEROPTV(
      vec[0] * get(0, 0) + vec[1] * get(0, 1) + vec[2] * get(0, 2) + 1 * get(0, 3),   // NOLINT(whitespace/line_length)
      vec[0] * get(1, 0) + vec[1] * get(1, 1) + vec[2] * get(1, 2) + 1 * get(1, 3),   // NOLINT(whitespace/line_length)
      vec[0] * get(2, 0) + vec[1] * get(2, 1) + vec[2] * get(2, 2) + 1 * get(2, 3));  // NOLINT(whitespace/line_length)
  }

  /**
   * Multiply two matrices together.
   * @param rhs the other matrix
   * @returns the product matrix
   */
  constexpr Matrix4 operator*(const Matrix4& rhs) const {
    return {{
      get(0, 0)*rhs.get(0, 0) + get(0, 1)*rhs.get(1, 0) + get(0, 2)*rhs.get(2, 0) + get(0, 3)*rhs.get(3, 0),  // NOLINT(whitespace/line_length)
        get(0, 0)*rhs.get(0, 1) + get(0, 1)*rhs.get(1, 1) + get(0, 2)*rhs.get(2, 1) + get(0, 3)*rhs.get(3, 1),  // NOLINT(whitespace/line_length)
          get(0, 0)*rhs.get(0, 2) + get(0, 1)*rhs.get(1, 2) + get(0, 2)*rhs.get(2, 2) + get(0, 3)*rhs.get(3, 2),  // NOLINT(whitespace/line_length)
            get(0, 0)*rhs.get(0, 3) + get(0, 1)*rhs.get(1, 3) + get(0, 2)*rhs.get(2, 3) + get(0, 3)*rhs.get(3, 3),  // NOLINT(whitespace/line_length)

      get(1, 0)*rhs.get(0, 0) + get(1, 1)*rhs.get(1, 0) + get(1, 2)*rhs.get(2, 0) + get(1, 3)*rhs.get(3, 0),  // NOLINT(whitespace/line_length)
        get(1, 0)*rhs.get(0, 1) + get(1, 1)*rhs.get(1, 1) + get(1, 2)*rhs.get(2, 1) + get(1, 3)*rhs.get(3, 1),  // NOLINT(whitespace/line_length)
          get(1, 0)*rhs.get(0, 2) + get(1, 1)*rhs.get(1, 2) + get(1, 2)*rhs.get(2, 2) + get(1, 3)*rhs.get(3, 2),  // NOLINT(whitespace/line_length)
            get(1, 0)*rhs.get(0, 3) + get(1, 1)*rhs.get(1, 3) + get(1, 2)*rhs.get(2, 3) + get(1, 3)*rhs.get(3, 3),  // NOLINT(whitespace/line_length)

      get(2, 0)*rhs.get(0, 0) + get(2, 1)*rhs.get(1, 0) + get(2, 2)*rhs.get(2, 0) + get(2, 3)*rhs.get(3, 0),  // NOLINT(whitespace/line_length)
        get(2, 0)*rhs.get(0, 1) + get(2, 1)*rhs.get(1, 1) + get(2, 2)*rhs.get(2, 1) + get(2, 3)*rhs.get(3, 1),  // NOLINT(whitespace/line_length)
          get(2, 0)*rhs.get(0, 2) + get(2, 1)*rhs.get(1, 2) + get(2, 2)*rhs.get(2, 2) + get(2, 3)*rhs.get(3, 2),  // NOLINT(whitespace/line_length)
            get(2, 0)*rhs.get(0, 3) + get(2, 1)*rhs.get(1, 3) + get(2, 2)*rhs.get(2, 3) + get(2, 3)*rhs.get(3, 3),  // NOLINT(whitespace/line_length)

      get(3, 0)*rhs.get(0, 0) + get(3, 1)*rhs.get(1, 0) + get(3, 2)*rhs.get(2, 0) + get(3, 3)*rhs.get(3, 0),  // NOLINT(whitespace/line_length)
        get(3, 0)*rhs.get(0, 1) + get(3, 1)*rhs.get(1, 1) + get(3, 2)*rhs.get(2, 1) + get(3, 3)*rhs.get(3, 1),  // NOLINT(whitespace/line_length)
          get(3, 0)*rhs.get(0, 2) + get(3, 1)*rhs.get(1, 2) + get(3, 2)*rhs.get(2, 2) + get(3, 3)*rhs.get(3, 2),  // NOLINT(whitespace/line_length)
            get(3, 0)*rhs.get(0, 3) + get(3, 1)*rhs.get(1, 3) + get(3, 2)*rhs.get(2, 3) + get(3, 3)*rhs.get(3, 3),  // NOLINT(whitespace/line_length)
    }};
  }

  /**
   * Multiply an array of points.
   * @param arr array of points to multiply
   * @returns a constexpr_array of transformed points
   */
  template <size_t count>
  constexpr constexpr_array<ZEROPTV, count> operator*(const ZEROPTV (&arr)[count]) const {  // NOLINT(whitespace/line_length)
#ifdef __ICC
    return array_mul_helper(arr, build_indices<count>::base_indices());
#else
    return array_mul_helper(arr, build_indices<count>());
#endif
  }

  /**
   * Multiply a constexpr_array of points.
   * @param arr array to multiply
   * @returns a constexpr_array of transformed points
   */
  template <size_t count>
  constexpr constexpr_array<ZEROPTV, count> operator*(const constexpr_array<ZEROPTV, count>& arr) const {  // NOLINT(whitespace/line_length)
#ifdef __ICC
    return array_mul_helper(arr, build_indices<count>::base_indices());
#else
    return array_mul_helper(arr, build_indices<count>());
#endif
  }

  /**
   * Build a rotation matrix around the X axis.
   * @param angle angle to rotate (in degrees)
   * @returns the rotation matrix
   */
  static constexpr Matrix4 rotateX(double angle) {
    return {{
      1, 0, 0, 0,
      0, cos(angle), -sin(angle), 0,
      0, sin(angle), cos(angle), 0,
      0, 0, 0, 1
    }};
  }

  /**
   * Build a rotation matrix around the Y axis.
   * @param angle angle to rotate (in degrees)
   * @returns the rotation matrix
   */
  static constexpr Matrix4 rotateY(double angle) {
    return {{
      cos(angle), 0, sin(angle), 0,
      0, 1, 0, 0,
      -sin(angle), 0, cos(angle), 0,
      0, 0, 0, 1
    }};
  }

  /**
   * Build a rotation matrix around the Z axis.
   * @param angle angle to rotate (in degrees)
   * @returns the rotation matrix
   */
  static constexpr Matrix4 rotateZ(double angle) {
    return {{
      cos(angle), -sin(angle), 0, 0,
      sin(angle), cos(angle), 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1
    }};
  }

  /**
   * Build a translation matrix.
   * @param x x translation
   * @param y y translation
   * @param z z translation
   * @returns the translation matrix
   */
  static constexpr Matrix4 translate(double x, double y, double z) {
    return {{
      1, 0, 0, x,
      0, 1, 0, y,
      0, 0, 1, z,
      0, 0, 0, 1
    }};
  }

  /**
   * Build a scale matrix.
   * @param sx x scale
   * @param sy y scale
   * @param sz z scale
   * @returns the scale matrix
   */
  static constexpr Matrix4 scale(double sx, double sy, double sz) {
    return {{
      sx, 0, 0, 0,
      0, sy, 0, 0,
      0, 0, sz, 0,
      0, 0, 0, 1
    }};
  }

 private:
  // could use compile-time taylor series expansion to calculate arbitrary
  // values, but I am not that smart and only need a few values
  /**
   * Get the sine of an angle in degrees.
   * Only works for a limited set of angles.
   * @param angle angle in degrees
   * @returns sin(angle)
   */
  static constexpr double sin(double angle) {
    return (angle == 90 ? 1 :
            angle == 0 ? 0 :
            angle == -90 ? -1 :
            angle == 180 ? 0 :
            angle == 45 ? 0.7071067811865475244008443621 :
            angle == -45 ? -0.7071067811865475244008443621 :
            angle == 135 ? 0.707106781186547524400844 :
            std::numeric_limits<T>::signaling_NaN());
  }

  /**
   * Get the cosine of an angle in degrees.
   * Only works for a limited set of angles.
   * @param angle angle in degrees
   * @returns cos(angle)
   */
  static constexpr double cos(double angle) {
    return (angle == 90 ? 0 :
            angle == 0 ? 1 :
            angle == -90 ? 0 :
            angle == 180 ? -1 :
            angle == 45 ? 0.7071067811865475244008443621 :
            angle == -45 ? 0.7071067811865475244008443621 :
            angle == 135 ? -0.7071067811865475244008443621 :
            std::numeric_limits<T>::signaling_NaN());
  }

  // http://loungecpp.wikidot.com/tips-and-tricks%3aindices
  /**
   * Compile-time recursion helper.
   * Just holds template information.
   */
  template <std::size_t... Is>
  struct indices {};

  /**
   * Compile-time recursion helper.
   * Just holds template information.
   */
  template <std::size_t N, std::size_t... Is>
  struct build_indices: build_indices<N-1, N-1, Is...> {
  };

  /**
   * Compile-time recursion helper.
   * Just holds template information.
   */
  template <std::size_t... Is>
  struct build_indices<0, Is...>: indices<Is...> {
    // work-around for Intel constexpr conversion bug
    typedef indices<Is...> base_indices;
  };

  /**
   * Compile-time recursion helper.
   * @param vecs array of points
   * @returns transformed points
   */
  template<size_t size, size_t... Is>
  constexpr constexpr_array<ZEROPTV, size> array_mul_helper(const ZEROPTV (&vecs)[size], indices<Is...>) const {  // NOLINT(whitespace/line_length)
    return {{ (*this * vecs[Is])... }};
  }

  /**
   * Compile-time recursion helper.
   * @param vecs constexpr_array of points
   * @returns transformed points
   */
  template<size_t size, size_t... Is>
  constexpr constexpr_array<ZEROPTV, size> array_mul_helper(const constexpr_array<ZEROPTV, size> vecs, indices<Is...>) const {  // NOLINT(whitespace/line_length)
    return {{ (*this * vecs[Is])... }};
  }
};

}  // namespace TALYFEMLIB
