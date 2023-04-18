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
#ifndef COMMON_INDEXER3D_H_
#define COMMON_INDEXER3D_H_

namespace TALYFEMLIB {

/**
 * A simple indexer class to convert 3 grid indices into a single index.
 *
 * Sample usage is:
 * int nx = 3;
 * int ny = 2;
 * int nz = 5;
 * Indexer3D indexer(nx, ny, nz);
 * for (int i = 0; i < nx; i++) {
 *   for (int j = 0; j < ny; j++) {
 *     for (int k = 0; k < nz; k++) {
 *       int index = indexer(i, j, k);
 *     }
 *   }
 * }
 */
class Indexer3D {
 public:
  /**
   * Create the indexer with the given dimensions
   *
   * @param dim_x_value number of items in x direction
   * @param dim_y_value number of items in y direction
   * @param dim_z_value number of items in z direction
   */
  Indexer3D(int dim_x_value, int dim_y_value, int dim_z_value)
      : dim_x_(dim_x_value),
        dim_y_(dim_y_value),
        dim_z_(dim_z_value) { }

  ~Indexer3D() { }

  /**
   * Returns the index calculated with the given values
   *
   * @param x_index x value
   * @param y_index y value
   * @param z_index z value
   * @return the index calculated from teh given values
   */
  int operator()(int x_index, int y_index, int z_index) const {
    return x_index + y_index * dim_x_ + z_index * dim_x_ * dim_y_;
  }

  /**
   * Returns the number of items in the x direction
   */
  inline int dim_x() { return dim_x_; }

  /**
   * Returns the number of items in the y direction
   */
  inline int dim_y() { return dim_y_; }

  /**
   * Returns the number of items in the z direction
   */
  inline int dim_z() { return dim_z_; }

 private:
  const int dim_x_;  ///< x dimension of grid
  const int dim_y_;  ///< y dimension of grid
  const int dim_z_;  ///< y dimension of grid (not used)
};


}  // namespace TALYFEMLIB

#endif  // COMMON_INDEXER3D_H_
