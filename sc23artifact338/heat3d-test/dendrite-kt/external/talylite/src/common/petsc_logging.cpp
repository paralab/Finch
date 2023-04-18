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
#include <talyfem/common/petsc_logging.h>

#include <assert.h>

#include <talyfem/utils/utils.h>

// UNUSED(x) is to prevent "variable unused" warnings when asserts are off
#define UNUSED(x) (void)(x)

// this should only be used for paranoia-tier error checking, where an error
// is incredibly unlikely and also won't affect the program if it is missed
#define CHKERRASSERT(err) UNUSED(err); assert(err == 0);

namespace TALYFEMLIB {

PetscStageLogger::PetscStageLogger(const char* name) {
  PetscErrorCode err = PetscLogStageRegister(name, &stage_);
  CHKERRASSERT(err);

  PrintInfo("Entering stage ", name, ".");
  err = PetscLogStagePush(stage_);
  CHKERRASSERT(err);
}

PetscStageLogger::~PetscStageLogger() {
  PetscErrorCode err = PetscLogStagePop();
  CHKERRASSERT(err);
}



PetscEventLogger::PetscEventLogger(const char* name, PetscClassId class_id) {
  PetscErrorCode err = PetscLogEventRegister(name, class_id, &event_);
  CHKERRASSERT(err);

  err = PetscLogEventBegin(event_, 0, 0, 0, 0);
  CHKERRASSERT(err);
}

PetscEventLogger::~PetscEventLogger() {
  PetscErrorCode err = PetscLogEventEnd(event_, 0, 0, 0, 0);
  CHKERRASSERT(err);
}

}  // namespace TALYFEMLIB

