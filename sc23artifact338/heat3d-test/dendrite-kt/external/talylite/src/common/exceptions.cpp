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
#include <talyfem/common/exceptions.h>
#include <talyfem/utils/utils.h>

namespace TALYFEMLIB {

void TALYException::print() const {
  PrintError(what());
}

PetscException::PetscException() {
  msg_ = "A PETSc error occured!";
}

PetscErrorCode PetscExceptionErrorHandler(MPI_Comm comm, int line,
    const char* func, const char* file, PetscErrorCode code,
    PetscErrorType type, const char* msg, void* ctx) {
  PetscTraceBackErrorHandler(comm, line, func, file, code, type, msg, ctx);
  throw PetscException();
}

}  // namespace TALYFEMLIB

