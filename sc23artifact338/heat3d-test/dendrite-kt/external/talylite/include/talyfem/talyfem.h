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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petscksp.h>
#import <PETSc/petscsnes.h>
#else
#include "petscksp.h"
#include "petscsnes.h"
#endif

#include <talyfem/utils/utils.h>
#include <talyfem/utils/test_utils.h>
#include <talyfem/utils/timers.h>
#include <talyfem/math/math.h>
#include <talyfem/input_data/input_data.h>
#include <talyfem/common/comm.h>
#include <talyfem/common/exceptions.h>
#include <talyfem/common/pack_comm.h>
#include <talyfem/common/petsc_logging.h>
#include <talyfem/grid/grid-types.h>
#include <talyfem/grid/nodedata.h>
#include <talyfem/domain_decomposition/mesh_partition.h>
#include <talyfem/integrator/function_integrator.h>
#include <talyfem/integrator/value_function.h>
#include <talyfem/fem/cequation.h>
#include <talyfem/fem/nonlinear_equation.h>
#include <talyfem/fem/preallocator_perfect.h>
#include <talyfem/stabilizer/tezduyar_upwind.h>
#include <talyfem/file_io/file_io.h>

#pragma GCC diagnostic pop

using namespace TALYFEMLIB;
