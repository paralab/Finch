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
#include <talyfem/common/comm.h>

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petscksp.h>
#import <PETSc/petscvec.h>
#import <PETSc/petsc.h>
#else
#include <petscksp.h>
#include <petscvec.h>
#include <petsc.h>
#endif

#include <talyfem/utils/utils.h>


namespace TALYFEMLIB {
static int calseconds = 3600 * 24 * 7;
static time_t tm_start;

// static void DefaultSet()
// {
//   time(&tm_start);
// }

double maxValueInAllProc(double maxV) {
  Vec maxVg;
  VecCreateMPI(PETSC_COMM_WORLD, 1, PETSC_DETERMINE, &maxVg);
  PetscMPIInt rank = GetMPIRank();
  PetscInt pos = rank;
  VecSetValues(maxVg, 1, &pos, &maxV, INSERT_VALUES);
  VecAssemblyBegin(maxVg);
  VecAssemblyEnd(maxVg);
  PetscInt index;
  VecMax(maxVg, &index, &maxV);
  VecDestroy(&maxVg);

  return maxV;
}

double minValueInAllProc(double minV) {
  return -maxValueInAllProc(-minV);
}

// Do a global sum of the given variable *in place*
void DoGlobalSum(double &value, int length) {
  MPI_Allreduce(MPI_IN_PLACE, &value, length, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
}

// Do a global sum of the given variable *in place*
void DoGlobalSum(int &value, int length) {
  MPI_Allreduce(MPI_IN_PLACE, &value, length, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
}

void SetCalculationSeconds(int seconds) {
  calseconds = seconds;
  time(&tm_start);
}

int IsCalculationTimeOut() {
  int rank = GetMPIRank();

  double tmp = 0;
  if (rank == 0) {
    time_t tm1;
    time(&tm1);
    if (static_cast<int>(difftime(tm1, tm_start) > calseconds)) {
      tmp = 1;
    }
  }

  tmp = maxValueInAllProc(tmp);
  if (tmp > 0.5)
    return 1;
  return 0;
}

}  // namespace TALYFEMLIB
