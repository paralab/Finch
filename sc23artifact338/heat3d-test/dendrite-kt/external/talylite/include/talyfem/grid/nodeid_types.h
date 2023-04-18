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

#include <vector>

// Index into the node array for an element.
typedef int ElemNodeID;

// Index of node on the local processor.
typedef int LocalNodeID;

// Global physical ID.
typedef PetscInt PhysicalNodeID;

// Index of node in PETSc data structures.
typedef PetscInt SolutionNodeID;

// Local variable index for local structures
typedef int LocalVarIdx;

// Local variable index for petsc structures
typedef PetscInt GlobalVarIdx;

inline std::vector<LocalNodeID> phys_to_local(
                                const std::vector<PhysicalNodeID>& phys) {
#ifndef PETSC_USE_64BIT_INDICES
  return phys;
#else
  std::vector<LocalNodeID> lcl(phys.size());
  for (unsigned int i = 0; i < phys.size(); i++) {
    lcl[i] = static_cast<LocalNodeID>(phys[i]);
  }
  return lcl;
#endif
}
