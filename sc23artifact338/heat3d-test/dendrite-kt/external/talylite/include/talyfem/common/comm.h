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
#ifndef COMMON_COMM_H_
#define COMMON_COMM_H_

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petsc.h>
#else
#include <petsc.h>
#endif


namespace TALYFEMLIB {

/**
 * Get the maximum maxV across all processes.
 * This is a collective call (all processors must call simultaneously).
 * @param maxV this process's maximum value
 * @returns maximum maxV across all processes
 */
double maxValueInAllProc(double maxV);

/**
 * Get the minimum minV across all processes.
 * This is a collective call (all processors must call simultaneously).
 * @param minV this process's minimum value
 * @returns minimum minV across all processes
 */
double minValueInAllProc(double minV);

/**
 * Calculates, in place, the sum of double value(s) from all processors.
 *
 * This is a collective operationa and all processors will have the result.
 *
 * @param[inout] value array of values to sum over processors
 * @param length number of items to sum (length of array)
 */
void DoGlobalSum(double &value, int length);

/**
 * Calculates, in place, the sum of integer value(s) from all processors.
 *
 * This is a collective operationa and all processors will have the result.
 *
 * @param[inout] value array of values to sum over processors
 * @param length number of items to sum (length of array)
 */
void DoGlobalSum(int &value, int length);

/**
 * Used to set timeout value.
 * Deprecated.
 * @param seconds number of seconds to time out after
 */
void SetCalculationSeconds(int seconds);

/**
 * Used to check if the timeout value has been exceeded.
 * Deprecated.
 * @returns 1 if calculation has timed out, 0 otherwise
 */
int IsCalculationTimeOut();

}  // namespace TALYFEMLIB

#endif  // COMMON_COMM_H_
