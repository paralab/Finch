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

#include <talyfem/fem/preallocator.h>
#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/node.h>

namespace TALYFEMLIB {

/**
 * Original implementation of the preallocator.
 * Faster than the perfect preallocator (and requires much less
 * inter-process communication), but also much less accurate
 * (overestimate significantly in some scenarios).
 */
class PreallocatorOriginal : public Preallocator {
 public:
  PreallocatorOriginal();
  virtual ~PreallocatorOriginal();

  /**
   * This equation specifies the grid and number of freedom per node.
   * @param pGrid The pointer of the grid
   * @param ndof The number of freedom per node.
   */
  virtual void redim(GRID* pGrid, int ndof);

  virtual void PresetPeriodicData(int n_per_bounds, int n_per_vars,
                                  int n_total_vars);

  virtual const PetscInt* get_dnz() const;
  virtual const PetscInt* get_onz() const;

  /**
   * Frees the memory allocated by this class.
   */
  void Destroy();

  virtual void calc_preallocation(Mat& mat, PetscInt bs);

 protected:
  GRID* p_grid_;  ///< pointer to grid object to preallocate for

  int n_dof_;  ///< The number of degrees of freedom per node.

  PetscInt *dnz;  ///< Array that stores the number of nonzeroes in the diagonal
                  ///< part of matrix Ag.

  PetscInt *onz;  ///< Array that stores the number of nonzeroes in the
                  ///< off-diagonal part of matrix Ag.

  // Preset periodic information. These values are set prior to initPerBC in
  // order to tell Preallocate about the periodic boundaries.
  int prealloc_n_periodic_vars;  ///< number of periodic variables per node
  int prealloc_n_periodic_bounds;  ///< number of periodic boundaries
  int prealloc_n_total_vars;  ///< number of total variables per node
};

}  // namespace TALYFEMLIB
