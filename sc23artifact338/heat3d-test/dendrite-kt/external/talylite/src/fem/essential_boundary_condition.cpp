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
#include <talyfem/fem/essential_boundary_condition.h>

namespace TALYFEMLIB {

EssentialBoundaryCondition::EssentialBoundaryCondition(LocalVarIdx index,
                                                       double matrix_value,
                                                       double value)
    : BaseBoundaryCondition(index),
      essential_value_(value),
      matrix_value_(matrix_value) {
}

void EssentialBoundaryCondition::Apply(Mat& Ag, Vec& bg,
                                       const ZEROARRAY<GlobalVarIdx>* pXMap,
                                       bool recalc_matrix) const {
  GlobalVarIdx petsc_index = PetscIndex(pXMap);
  if (recalc_matrix) {
    MatSetValue(Ag, petsc_index, petsc_index, matrix_value_, INSERT_VALUES);
  }
  VecSetValue(bg, petsc_index, essential_value_, INSERT_VALUES);
}

void EssentialBoundaryCondition::ApplyEssBCToSolution(
    ZEROARRAY<double>& solution) const {
  solution(index_) = essential_value_ / matrix_value_;
}

}  // namespace TALYFEMLIB
