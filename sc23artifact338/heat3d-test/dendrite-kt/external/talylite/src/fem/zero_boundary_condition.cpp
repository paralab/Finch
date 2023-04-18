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
#include <talyfem/fem/zero_boundary_condition.h>


namespace TALYFEMLIB {

void ZeroBoundaryCondition::Apply(Mat& Ag, Vec& bg,
                                  const ZEROARRAY<GlobalVarIdx>* pXMap,
                                  bool recalc_matrix) const {
  GlobalVarIdx petsc_index = PetscIndex(pXMap);
  if (recalc_matrix) {
    MatSetValue(Ag, petsc_index, petsc_index, 1.0, INSERT_VALUES);
  }
  VecSetValue(bg, petsc_index, 0.0, INSERT_VALUES);
}

void ZeroBoundaryCondition::ApplyEssBCToSolution(
    ZEROARRAY<double>& solution) const {
  solution(index_) = 0.0;
}

}  // namespace TALYFEMLIB
