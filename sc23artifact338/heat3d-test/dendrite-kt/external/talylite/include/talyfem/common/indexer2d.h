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
#ifndef COMMON_INDEXER2D_H_
#define COMMON_INDEXER2D_H_

namespace TALYFEMLIB {

/**
 * A simple indexer class to convert 2 grid indices into a single index.
 *
 * Sample usage is:
 * int nx = 3;
 * int ny = 2;
 * Indexer2D indexer(nx, ny);
 * for (int i = 0; i < nx; i++) {
 *   for (int j = 0; j < ny; j++) {
 *     int index = indexer(i, j);
 *   }
 * }
 */
class Indexer2D {
 public:
  /**
   * Create the indexer with the given dimensions
   *
   * @param dim_x_value number of items in x direction
   * @param dim_y_value number of items in y direction
   */
  Indexer2D(int dim_x_value, int dim_y_value)
      : dim_x_(dim_x_value),
        dim_y_(dim_y_value) { }

  ~Indexer2D() { }

  /**
   * Returns the index calculated with the given values
   *
   * @param x_index x value
   * @param y_index y value
   * @return the index calculated from teh given values
   */
  int operator()(int x_index, int y_index) const {
    return x_index + y_index * dim_x_;
  }

  /**
   * Returns the number of items in the x direction
   */
  inline int dim_x() { return dim_x_; }

  /**
   * Returns the number of items in the y direction
   */
  inline int dim_y() { return dim_y_; }

 private:
  const int dim_x_;  // x dimension of grid
  const int dim_y_;  // y dimension of grid (not used)
};

}  // namespace TALYFEMLIB

#endif  // COMMON_INDEXER2D_H_
