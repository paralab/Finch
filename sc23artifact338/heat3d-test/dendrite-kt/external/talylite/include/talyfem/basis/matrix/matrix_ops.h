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

#include <type_traits>

#include "matrix.h"

namespace MatrixUtils {

// determinant
/**
 * @returns the determinant of the matrix
 */
template <typename MatType>
constexpr typename std::enable_if<MatType::n == 1, double>::type calc_det(const MatType& mat) {
  return mat(0, 0);
}

template <typename MatType>
constexpr typename std::enable_if<MatType::n == 2, double>::type calc_det(const MatType& mat) {
  return mat(0, 0) * mat(1, 1) - mat(1, 0) * mat(0, 1);
}

template <typename MatType>
constexpr typename std::enable_if<MatType::n == 3, double>::type calc_det(const MatType& mat) {
  return mat(0, 0) * calc_det(mat.template create_view<0, 0>()) - mat(0, 1) * calc_det(mat.template create_view<0, 1>())
      + mat(0, 2) * calc_det(mat.template create_view<0, 2>());
}

template <typename MatType>
constexpr typename std::enable_if<MatType::n == 4, double>::type calc_det(const MatType& mat) {
  return mat(0, 0) * calc_det(mat.template create_view<0, 0>()) - mat(0, 1) * calc_det(mat.template create_view<0, 1>())
      + mat(0, 2) * calc_det(mat.template create_view<0, 2>()) - mat(0, 3) * calc_det(mat.template create_view<0, 3>());
}

template <typename MatType>
constexpr typename std::enable_if<MatType::n == 5, double>::type calc_det(const MatType& mat) {
  return mat(0, 0) * calc_det(mat.template create_view<0, 0>()) - mat(0, 1) * calc_det(mat.template create_view<0, 1>())
      + mat(0, 2) * calc_det(mat.template create_view<0, 2>()) - mat(0, 3) * calc_det(mat.template create_view<0, 3>())
      + mat(0, 4) * calc_det(mat.template create_view<0, 4>());
}

// helper for passing in raw 2D arrays
/**
 * @returns the determinant of a 2D array (matrix)
 */
template <int n>
constexpr double calc_det(const double (&mat)[n][n]) {
  return calc_det(Matrix(mat));
}

/**
 * Calculates the cofactor matrix.
 * @param mat matrix to calculate the cofactor of
 * @param[out] cof_out cofactor matrix (output)
 */
template <typename MatType>
typename std::enable_if<MatType::n == 1>::type calc_cof_t(const MatType& mat, double (&cof_out)[MatType::n][MatType::n]) {
  cof_out[0][0] = 1.0;
}

template <typename MatType>
typename std::enable_if<MatType::n == 2>::type calc_cof_t(const MatType& mat, double (&cof_out)[MatType::n][MatType::n]) {
  cof_out[0][0] = calc_det(mat.template create_view<0, 0>()); //mat(1, 1)
  cof_out[0][1] = -calc_det(mat.template create_view<0, 1>()); //-mat(0, 1);
  cof_out[1][0] = -calc_det(mat.template create_view<1, 0>()); //-mat(1, 0);
  cof_out[1][1] = calc_det(mat.template create_view<1, 1>()); //mat(0, 0);
}

template <typename MatType>
typename std::enable_if<MatType::n == 3>::type calc_cof_t(const MatType& mat, double (&cof_out)[MatType::n][MatType::n]) {
  cof_out[0][0] = calc_det(mat.template create_view<0, 0>());  // ok
  cof_out[0][1] = -calc_det(mat.template create_view<0, 1>());
  cof_out[0][2] = calc_det(mat.template create_view<0, 2>());

  cof_out[1][0] = -calc_det(mat.template create_view<1, 0>());
  cof_out[1][1] = calc_det(mat.template create_view<1, 1>());  // ok
  cof_out[1][2] = -calc_det(mat.template create_view<1, 2>());

  cof_out[2][0] = calc_det(mat.template create_view<2, 0>());
  cof_out[2][1] = -calc_det(mat.template create_view<2, 1>());
  cof_out[2][2] = calc_det(mat.template create_view<2, 2>());  // ok
}

template <typename MatType>
typename std::enable_if<MatType::n == 4>::type calc_cof_t(const MatType& mat, double (&cof_out)[MatType::n][MatType::n]) {
  cof_out[0][0] = calc_det(mat.template create_view<0, 0>());
  cof_out[0][1] = -calc_det(mat.template create_view<0, 1>());
  cof_out[0][2] = calc_det(mat.template create_view<0, 2>());
  cof_out[0][3] = -calc_det(mat.template create_view<0, 3>());

  cof_out[1][0] = -calc_det(mat.template create_view<1, 0>());
  cof_out[1][1] = calc_det(mat.template create_view<1, 1>());
  cof_out[1][2] = -calc_det(mat.template create_view<1, 2>());
  cof_out[1][3] = calc_det(mat.template create_view<1, 3>());

  cof_out[2][0] = calc_det(mat.template create_view<2, 0>());
  cof_out[2][1] = -calc_det(mat.template create_view<2, 1>());
  cof_out[2][2] = calc_det(mat.template create_view<2, 2>());
  cof_out[2][3] = -calc_det(mat.template create_view<2, 3>());

  cof_out[3][0] = -calc_det(mat.template create_view<3, 0>());
  cof_out[3][1] = calc_det(mat.template create_view<3, 1>());
  cof_out[3][2] = -calc_det(mat.template create_view<3, 2>());
  cof_out[3][3] = calc_det(mat.template create_view<3, 3>());
}

template <typename MatType>
typename std::enable_if<MatType::n == 5>::type calc_cof_t(const MatType& mat, double (&cof_out)[MatType::n][MatType::n]) {
  cof_out[0][0] = calc_det(mat.template create_view<0, 0>());
  cof_out[0][1] = -calc_det(mat.template create_view<0, 1>());
  cof_out[0][2] = calc_det(mat.template create_view<0, 2>());
  cof_out[0][3] = -calc_det(mat.template create_view<0, 3>());
  cof_out[0][4] = calc_det(mat.template create_view<0, 4>());

  cof_out[1][0] = -calc_det(mat.template create_view<1, 0>());
  cof_out[1][1] = calc_det(mat.template create_view<1, 1>());
  cof_out[1][2] = -calc_det(mat.template create_view<1, 2>());
  cof_out[1][3] = calc_det(mat.template create_view<1, 3>());
  cof_out[1][4] = -calc_det(mat.template create_view<1, 4>());

  cof_out[2][0] = calc_det(mat.template create_view<2, 0>());
  cof_out[2][1] = -calc_det(mat.template create_view<2, 1>());
  cof_out[2][2] = calc_det(mat.template create_view<2, 2>());
  cof_out[2][3] = -calc_det(mat.template create_view<2, 3>());
  cof_out[2][4] = calc_det(mat.template create_view<2, 4>());

  cof_out[3][0] = -calc_det(mat.template create_view<3, 0>());
  cof_out[3][1] = calc_det(mat.template create_view<3, 1>());
  cof_out[3][2] = -calc_det(mat.template create_view<3, 2>());
  cof_out[3][3] = calc_det(mat.template create_view<3, 3>());
  cof_out[3][4] = -calc_det(mat.template create_view<3, 4>());

  cof_out[4][0] = calc_det(mat.template create_view<4, 0>());
  cof_out[4][1] = -calc_det(mat.template create_view<4, 1>());
  cof_out[4][2] = calc_det(mat.template create_view<4, 2>());
  cof_out[4][3] = -calc_det(mat.template create_view<4, 3>());
  cof_out[4][4] = calc_det(mat.template create_view<4, 4>());
}

template <typename MatType>
typename std::enable_if<MatType::n == 6>::type calc_cof_t(const MatType& mat, double (&cof_out)[MatType::n][MatType::n]) {
  cof_out[0][0] = calc_det(mat.template create_view<0, 0>());
  cof_out[0][1] = -calc_det(mat.template create_view<0, 1>());
  cof_out[0][2] = calc_det(mat.template create_view<0, 2>());
  cof_out[0][3] = -calc_det(mat.template create_view<0, 3>());
  cof_out[0][4] = calc_det(mat.template create_view<0, 4>());
  cof_out[0][5] = -calc_det(mat.template create_view<0, 5>());

  cof_out[1][0] = -calc_det(mat.template create_view<1, 0>());
  cof_out[1][1] = calc_det(mat.template create_view<1, 1>());
  cof_out[1][2] = -calc_det(mat.template create_view<1, 2>());
  cof_out[1][3] = calc_det(mat.template create_view<1, 3>());
  cof_out[1][4] = -calc_det(mat.template create_view<1, 4>());
  cof_out[1][5] = calc_det(mat.template create_view<1, 5>());

  cof_out[2][0] = calc_det(mat.template create_view<2, 0>());
  cof_out[2][1] = -calc_det(mat.template create_view<2, 1>());
  cof_out[2][2] = calc_det(mat.template create_view<2, 2>());
  cof_out[2][3] = -calc_det(mat.template create_view<2, 3>());
  cof_out[2][4] = calc_det(mat.template create_view<2, 4>());
  cof_out[2][5] = -calc_det(mat.template create_view<2, 5>());

  cof_out[3][0] = -calc_det(mat.template create_view<3, 0>());
  cof_out[3][1] = calc_det(mat.template create_view<3, 1>());
  cof_out[3][2] = -calc_det(mat.template create_view<3, 2>());
  cof_out[3][3] = calc_det(mat.template create_view<3, 3>());
  cof_out[3][4] = -calc_det(mat.template create_view<3, 4>());
  cof_out[3][5] = calc_det(mat.template create_view<3, 5>());

  cof_out[4][0] = calc_det(mat.template create_view<4, 0>());
  cof_out[4][1] = -calc_det(mat.template create_view<4, 1>());
  cof_out[4][2] = calc_det(mat.template create_view<4, 2>());
  cof_out[4][3] = -calc_det(mat.template create_view<4, 3>());
  cof_out[4][4] = calc_det(mat.template create_view<4, 4>());
  cof_out[4][5] = -calc_det(mat.template create_view<4, 5>());

  cof_out[5][0] = -calc_det(mat.template create_view<5, 0>());
  cof_out[5][1] = calc_det(mat.template create_view<5, 1>());
  cof_out[5][2] = -calc_det(mat.template create_view<5, 2>());
  cof_out[5][3] = calc_det(mat.template create_view<5, 3>());
  cof_out[5][4] = -calc_det(mat.template create_view<5, 4>());
  cof_out[5][5] = calc_det(mat.template create_view<5, 5>());
}

// helper for passing in raw 2D arrays
/**
 * Calculates the cofactor matrix for a 2D array (matrix).
 * @param[in] mat matrix to calculate the cofactor of
 * @param[out] cof_out cofactor matrix output
 */
template <int n>
void calc_cof_t(const double (&mat)[n][n], double (&cof_out)[n][n]) {
  calc_cof_t(Matrix(mat), cof_out);
}

// calculate jacobian using cofactor matrix
// (since the cofactor matrix is needed elsewhere (dN), we can reuse it here)
/**
 * Calculate the determinant of the jacobian of a matrix using its previously
 * calculated cofactor matrix to reduce calculations.
 * @param[in] dxde matrix to calculate the determinant of
 * @param[in] cof cofactor matrix of dxde
 * @returns determinant of the jacobian of dxde
 */
template <int n>
constexpr typename std::enable_if<n == 1, double>::type calc_jacobian_w_cof(const double (&dxde)[n][n], const double (&cof)[n][n]) {
  return dxde[0][0];
}

template <int n>
constexpr typename std::enable_if<n == 2, double>::type calc_jacobian_w_cof(const double (&dxde)[n][n], const double (&cof)[n][n]) {
  return dxde[0][0] * dxde[1][1] - dxde[0][1] * dxde[1][0];
}

template <int n>
constexpr typename std::enable_if<n == 3, double>::type calc_jacobian_w_cof(const double (&dxde)[n][n], const double (&cof)[n][n]) {
  return dxde[0][0] * cof[0][0] + dxde[0][1] * cof[0][1]
      + dxde[0][2] * cof[0][2];
}

template <int n>
constexpr typename std::enable_if<n == 4, double>::type calc_jacobian_w_cof(const double (&dxde)[n][n], const double (&cof)[n][n]) {
  return dxde[0][0] * cof[0][0] + dxde[0][1] * cof[0][1]
      + dxde[0][2] * cof[0][2] + dxde[0][3] * cof[0][3];
}

template <int n>
constexpr typename std::enable_if<n == 6, double>::type calc_jacobian_w_cof(const double (&dxde)[n][n], const double (&cof)[n][n]) {
  return dxde[0][0] * cof[0][0] + dxde[0][1] * cof[0][1]
      + dxde[0][2] * cof[0][2] + dxde[0][3] * cof[0][3]
      + dxde[0][4] * cof[0][4] + dxde[0][5] * cof[0][5];
}

}  // namespace MatrixUtils
