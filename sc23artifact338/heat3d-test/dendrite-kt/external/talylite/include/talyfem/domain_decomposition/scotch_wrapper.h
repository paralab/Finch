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

/**
 * Scotch provides a ParMETIS compatibility library, libptscotchparmetis, which
 * has a wrapper for ParMETIS_V3_PartKway (same inputs, same outputs, just uses
 * Scotch to do the decomposition).
 *
 * We originally used ParMETIS_V3_Part*Mesh*Kway, which takes in a mesh,
 * converts it to a dual graph, and then runs the Kway partitioning algorithm.
 * As far as I can tell, Scotch does not provide a similar function, and it's
 * nontrivial to implement ourselves.
 *
 * The solution is to link with both ParMETIS and Scotch's ParMETIS
 * compatibility layer. We use ParMETIS's Mesh2Dual function to build
 * the graph, and Scotch's ParMETIS compatibility function (declared
 * in this file) to do the partitioning.
 *
 * We can't directly include both the real ParMETIS and scotch's ParMETIS
 * compatibility layer because both files have the same name. That's why
 * this file exists - to allow us to include both (and check that types
 * are roughly compatibile between the two).
 *
 * This file will appear blank unless the ENABLE_SCOTCH flag is specified.
 *
 * This file is only included by .cpp files in the library, so it is okay
 * if user code does not have a matching "ENABLE_SCOTCH" flag.
 */

#ifdef ENABLE_SCOTCH

#include <mpi.h>  // for MPI_Comm
#include <scotch.h>  // for SCOTCH_Num
#include <ptscotch.h>

static_assert(sizeof(idx_t) == sizeof(SCOTCH_Num),
              "ParMETIS idx_t and SCOTCH_Num are incompatible!");
extern "C" {

// declaration should be identical to the one in scotch's parmetis.h

void SCOTCH_ParMETIS_V3_PartKway(
const SCOTCH_Num * const    vtxdist,
SCOTCH_Num * const          xadj,
SCOTCH_Num * const          adjncy,
SCOTCH_Num * const          vwgt,
SCOTCH_Num * const          adjwgt,
const SCOTCH_Num * const    wgtflag,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    ncon,                 /* Not used */
const SCOTCH_Num * const    nparts,
const float * const         tpwgts,
const float * const         ubvec,                /* Not used */
const SCOTCH_Num * const    options,              /* Not used */
SCOTCH_Num * const          edgecut,
SCOTCH_Num * const          part,
MPI_Comm *                  comm);

}

#endif
