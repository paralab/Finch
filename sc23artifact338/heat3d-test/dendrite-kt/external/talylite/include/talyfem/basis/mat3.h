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
#include <assert.h>
#include <type_traits>  // for std::enable_if
#include <array>  // for std::array
#include <cmath>

namespace TALYFEMLIB {

/**
 * Compile-time sized square matrix.
 */
template <int n>
class Matrix {
 public:
  std::array<double, n*n> data;  ///< matrix data

  /**
   * @returns the identity matrix (diagonal is 1s)
   */
  static Matrix<n> identity() {
    Matrix<n> mat;
    for (int row = 0; row < n; row++) {
      for (int col = 0; col < n; col++) {
        mat(row, col) = (row == col) ? 1.0 : 0.0;
      }
    }
    return mat;
  }

  /**
   * Access a cell in the matrix.
   * @param i row
   * @param j column
   * @returns the value
   */
  double& operator()(int i, int j) {
    assert(i >= 0 && i < n);
    assert(j >= 0 && j < n);
    return data[i*n + j];
  }

  /**
   * Access a cell in the matrix (const).
   * @param i row
   * @param j column
   * @returns the value
   */
  double operator()(int i, int j) const {
    assert(i >= 0 && i < n);
    assert(j >= 0 && j < n);
    return data[i*n + j];
  }

  /**
   * Matrix-matrix multiplication.
   * @param rhs matrix to multiply by
   * @returns the matrix product
   */
  Matrix<n> operator*(const Matrix<n>& rhs) const {
    const Matrix<n>& lhs = *this;
    Matrix<n> out;
    for (unsigned int i = 0; i < n; i++) {  // row in lhs
      for (unsigned int j = 0; j < n; j++) {  // col in rhs
        double sum = 0.0;
        for (unsigned int k = 0; k < n; k++) {  // row in rhs
          sum += lhs(i, k) * rhs(k, j);
        }
        out(i, j) = sum;
      }
    }
    return out;
  }

  /**
   * @returns the transpose of this matrix
   */
  Matrix<n> transpose() const {
    Matrix<n> out;
    for (unsigned int i = 0; i < n; i++) {
      for (unsigned int j = 0; j < n; j++) {
        out(i, j) = (*this)(j, i);
      }
    }
    return out;
  }

  /**
   * matrix-matrix multiplication
   * @param rhs  matrix to multiply by
   * @returns the matrix product
   */
  template <int rhs_m>
  std::array<std::array<double, n>, rhs_m> operator*(
      const double(&rhs)[n][rhs_m]) const {
    const Matrix<n>& lhs = *this;
    std::array< std::array<double, n>, rhs_m> out;
    for (unsigned int i = 0; i < n; i++) {  // row in lhs
      for (unsigned int j = 0; j < rhs_m; j++) {  // col in rhs
        double sum = 0.0;
        for (unsigned int k = 0; k < n; k++) {  // row in rhs
          sum += lhs(i, k) * rhs(k, j);
        }
        out[i][j] = sum;
      }
    }
    return out;
  }

  /**
   * Builds a matrix that rotates a 3D point about axis by cos(angle).
   * Angle should be cos(angle_in_radians).
   * Too lazy to figure out a dimension-independent implementation for this
   * @param axis axis to rotate about
   * @param angle cos(angle) value (you must take the cosine yourself!)
   * @returns matrix that rotates about axis by angle radians
   */
  static Matrix<n> axisAngle(
      ZEROPTV axis, double angle) {
    assert(n == 3);
    const double c = angle;
    const double s = sin(acos(angle));
    const double ic = 1 - c;

    Matrix<n> m;
    m(0, 0) = axis.x() * axis.x() * ic + c;
    m(0, 1) = axis.x() * axis.y() * ic - axis.z() * s;
    m(0, 2) = axis.x() * axis.z() * ic + axis.y() * s;

    m(1, 0) = axis.y() * axis.x() * ic + axis.z() * s;
    m(1, 1) = axis.y() * axis.y() * ic + c;
    m(1, 2) = axis.y() * axis.z() * ic - axis.x() * s;

    m(2, 0) = axis.z() * axis.x() * ic - axis.y() * s;
    m(2, 1) = axis.z() * axis.y() * ic + axis.x() * s;
    m(2, 2) = axis.z() * axis.z() * ic + c;

    return m;
  }

  /**
   * Same algorithm as matrix-matrix multiplication, but j is assumed to be 0
   * (since rhs is 1x3).
   * @param rhs vector to multiply by
   * @returns matrix-vector product
   */
  ZEROPTV operator*(const ZEROPTV& rhs) const {
    assert(n == 3);

    const Matrix<n>& lhs = *this;
    ZEROPTV out;
    for (unsigned int i = 0; i < n; i++) {  // row in lhs
      double sum = 0.0;
      for (unsigned int k = 0; k < 3; k++) {  // row in rhs
        sum += lhs(i, k) * rhs(k);
      }
      out(i) = sum;
    }
    return out;
  }
};

/**
 * Helper for multiplying a Matrix by a 2D array.
 * @param lhs left hand side matrix
 * @param rhs right hand side matrix
 * @returns matrix-matrix product
 */
template <uint64_t lhs_n, uint64_t lhs_m>
std::array<std::array<double, lhs_m>, lhs_n> operator*(
    const double(&lhs)[lhs_n][lhs_m], const Matrix<lhs_m>& rhs) {
  static const int rhs_n = lhs_m;
  static const int rhs_m = rhs_n;
  assert(lhs_m == rhs_n);

  std::array< std::array<double, rhs_m>, lhs_n> out;
  for (unsigned int i = 0; i < lhs_n; i++) {  // row in lhs
    for (unsigned int j = 0; j < rhs_m; j++) {  // col in rhs
      double sum = 0.0;
      for (unsigned int k = 0; k < lhs_m; k++) {  // row in rhs
        sum += lhs[i][k] * rhs(k, j);
      }
      out[i][j] = sum;
    }
  }
  return out;
}

//! Convenience typedef for 3x3 matrix
typedef Matrix<3> Matrix3;

}  // namespace TALYFEMLIB
