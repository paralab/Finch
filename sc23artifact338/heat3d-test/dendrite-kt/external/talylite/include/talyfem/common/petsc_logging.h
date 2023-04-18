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


// C++-style wrappers for PETSc logging.
// Stages/events are pushed/began in the constructor and
// popped/ended in the destructor.

namespace TALYFEMLIB {

/**
 * C++-style wrapper for PETSc stage logging.
 * Use this to time "stages" of your program (assembly, post-processing, etc).
 *
 * Calls PetscLogStagePush() when constructed and PetscLogStagePop() in the
 * destructor. You can use this to log stages based on scope.
 *
 * Example:
 *
 * void calculateStuff() {
 *   PetscStageLogger("calculateStuff");
 *   // [calculations intensify]
 * }
 *
 * The stage will automatically end when calculateStuff finishes.
 * You can see stage timing information in the PETSc logging output.
 */
class PetscStageLogger {
 public:
  /**
   * Push a new PETSc stage.
   * Events will be grouped into this stage until it goes out of scope.
   * @param name Name of the stage.
   */
  explicit PetscStageLogger(const char* name);

  /**
   * Pops this stage.
   */
  ~PetscStageLogger();

 private:
  PetscLogStage stage_;  ///< PETSc stage ID
};

/**
 * C++-style wrapper for PETSc event logging.
 * Use this for events that happen multiple times (like CEquation::Integrands).
 *
 * Usage is the same as PetscStageLogger.
 */
class PetscEventLogger {
 public:
  /**
   * Start a new PETSc event.
   * Events can be used to gather timing information.
   * @param name Name of the event.
   * @param class_id PETSc class of the event, optional.
   */
  explicit PetscEventLogger(const char* name, PetscClassId class_id = 0);

  /**
   * End this PETSc event.
   */
  ~PetscEventLogger();

 private:
  PetscLogEvent event_;  ///< PETSc event ID
};

}  // namespace TALYFEMLIB
