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
#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petsc.h>
#else
#include <petsc.h>
#endif

namespace TALYFEMLIB {

/**
 * Base class for preallocators to determine memory requirement for A matrix.
 *
 * This class is responsible for determining how much memory is required in the
 * A matrix so that it can be correctly created. Inaccurate preallocation will
 * lead to either poor performance or wasted memory. See the PETSc user manual
 * for more details about the preallocation process.
 *
 * Subclasses need to implement the following:
 *  calc_preallocation - to do the actual preallocation calculation
 *  get_dnz - to return the diagonal non-zero counts for the matrix
 *  get_onz - to return the off diagonal non-zero counts for the matrix
 *
 */
class Preallocator {
 public:
  virtual ~Preallocator() {}

  /**
   * Calculate diagonal non-zeros and off-diagonal non-zeros.
   * This call is collective - that is, you should call this function
   * on every process simultaneously.
   * After this function is called, get_dnz() and get_onz() will return
   * the diagonal/off-diagonal non-zeros for this process for each row,
   * which are then passed into MatMPIAIJSetPreallocation() or similar.
   * @param mat Matrix we are preallocating for. Used to check which rows
   *            this process owns.
   * @param blockSize Block size when using block compressed row format.
   */
  virtual void calc_preallocation(Mat& mat, PetscInt blockSize = 1) = 0;

  /**
   * Return the diagonal non-zeros in the global PETSc matrix for the
   * matrix rows owned by this process.
   * Only valid after calc_preallocation() has been called.
   * @returns diagonal non-zeros by row
   *          (get_dnz()[0] == number of diagonal non-zeros for the first row
   *          owned by this process).
   */
  virtual const PetscInt* get_dnz() const = 0;

  /**
   * Return the off-diagonal non-zeros in the global PETSc matrix for the
   * matrix rows owned by this process.
   * Only valid after calc_preallocation() has been called.
   * @returns off-diagonal non-zeros by row
   *          (get_onz()[0] == number of off-diagonal non-zeros for the first row
   *          owned by this process).
   */
  virtual const PetscInt* get_onz() const = 0;


  /**
   * Tell the solver information about the periodic struction of the system
   *
   * Since the periodic bounary conditions are not set until after the
   * Preallocate function is called, Preallocate is not aware of the system's
   * periodic structure. Without that information, it cannot accurately
   * determine memory usage. This function should be called prior to
   * redim in order to ensure Preallocate knows about the periodic bounds.
   *
   * @param n_per_bounds number of periodic boundaries
   * @param n_per_vars number of periodic variables per node
   * @param n_total_vars number of total variables per node
   * @deprecated This function will be removed from the interface soon!
   */
  virtual void PresetPeriodicData(int n_per_bounds, int n_per_vars,
                                  int n_total_vars) {}
};

}  // namespace TALYFEMLIB

