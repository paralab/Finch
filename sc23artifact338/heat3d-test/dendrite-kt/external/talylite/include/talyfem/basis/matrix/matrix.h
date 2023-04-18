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

namespace MatrixUtils {

/**
 * A "view" into a matrix that eliminates one row and one column.
 * Views have a parent (which can also eliminate one row and one column),
 * and "crossed out" rows/columns are cumulative.
 * This construct is useful when calculating the determinant and cofactor
 * matrices. Data is never duplicated, and since the "skipped" rows/cols are
 * known at compile-time, operations usually optimize very well.
 */
template <typename ParentMatrixType, int parent_n, int skip_row, int skip_col>
struct MatrixView
{
  static_assert((skip_row >= 0 && skip_col >= 0), "Skips must be positive");

  static const int n = ParentMatrixType::n - 1;  ///< width/height of matrix
  typedef MatrixView<ParentMatrixType, parent_n, skip_row, skip_col> this_type;
  ///< type of this

  /**
   * @param parent_init parent
   */
  constexpr MatrixView(const ParentMatrixType& parent_init) : parent(parent_init) {}

  /**
   * @param i row
   * @param j column
   * @returns row i, col j of the matrix
   */
  double operator()(int i, int j) const {
    const int real_i = i + (i >= skip_row ? 1 : 0);
    const int real_j = j + (j >= skip_col ? 1 : 0);
    return parent(real_i, real_j);
  }

  /**
   * Create a child view of this matrix with i_skip, j_skip crossed out.
   * @returns child matrix (of size n - 1)
   */
  template <int i_skip, int j_skip>
  constexpr MatrixView<this_type, n, i_skip, j_skip> create_view() const {
    return MatrixView<this_type, n, i_skip, j_skip>(*this);
  }

 private:
  const ParentMatrixType parent;  ///< parent matrix
};

/**
 * "Base case" for a MatrixView, with no rows/columns crossed out.
 */
template <int data_dim>
struct MatrixView<void, data_dim, -1, -1> {

  static const int n = data_dim;  ///< width/height of matrix
  typedef MatrixView<void, data_dim, -1, -1> this_type;
  ///< type of this

  /**
   * @param data_init matrix data
   */
  constexpr MatrixView(const double (&data_init)[data_dim][data_dim]) : data(data_init) {}

  /**
   * @param i row
   * @param j column
   * @returns row i, column j of the matrix
   */
  double operator()(int i, int j) const {
    return data[i][j];
  }

  /**
   * @returns a "child" of this matrix with row i_skip, column j_skip removed
   */
  template <int i_skip, int j_skip>
  constexpr MatrixView<this_type, n, i_skip, j_skip> create_view() const {
    return MatrixView<this_type, n, i_skip, j_skip>(*this);
  }

 private:
  const double (&data)[data_dim][data_dim];  ///< our data
};

/**
 * @param data matrix data
 * @returns a "base" MatrixView for a 2D array of data
 */
template <int n>
MatrixView<void, n, -1, -1> Matrix(const double (&data)[n][n]) {
  return MatrixView<void, n, -1, -1>(data);
}

}  // namespace MatrixUtils
