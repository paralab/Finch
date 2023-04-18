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
#ifndef GRID_CREMOTENEIGHBOR_H_
#define GRID_CREMOTENEIGHBOR_H_

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petscsys.h>
#else
#include <petscsys.h>
#endif



namespace TALYFEMLIB {


/**
 * Structure used for date transfer about neighboring elements across
 * domain boundary. It stores data needed when the MPI_Alltoallv()
 * function is called.
 *
 * Example of usage:
 *
 *     for(int i = 0; i < pGrid->remoteNeighbors.ngbrno(); i++) {
 *       counts[pGrid->remoteNeighbors.remote_rank(i)] += 1;
 *     }
 *     int counter = 0;
 *     for(int i = 0; i < pGrid->size; i++) {
 *       displs[i] += counter;
 *       counter += counts[i];
 *     }
 *     MPI_Alltoallv(localFaceData, counts, displs, type,
 *                   remoteFaceData, counts, displs, type,
 *                   PETSC_COMM_WORLD);
 */
class CRemoteNeighbor {
 public:
  static const PetscInt nNgbrVar;  ///< Number of fields in element neighbor
                                   ///< description

  /**
   * Fills object with data.
   * @param nsr number of elements' neighbors that are in other domain
   */
  void redim(PetscInt nsr) {
    free();
    nShrNgbrs = nsr;
    PetscMalloc(nShrNgbrs * nNgbrVar * sizeof(PetscInt), &shrNgbrs);
  }

  /**
   * Cleans the allocated memory.
   */
  void free() {
    PetscFree(shrNgbrs);
    nShrNgbrs = 0;
  }

  /**
   * Gives access to data (it is used in communication across domains).
   * This is just an array organized as follows:
   *   - global ID of an element in a current rank,
   *   - local ID of an element in a current rank,
   *   - ID of a current rank,
   *   - which face of an element is actually shared,
   *   - global ID of a neighbor element in a remote rank,
   *   - local ID of a neighbor element in a remote rank,
   *   - ID of a remote rank (where current element neighbor is),
   *   - which face the remote neighbor is actually shared.
   * @return pointer to array organized in fields of `nNgbrVar` size
   */
  PetscInt* data() {
    return shrNgbrs;
  }

  /**
   * Returns ID of an element in global numeration.
   * @param idx position of element in remote neighbors list
   * @return global ID of an element
   */
  inline PetscInt local_globalID(PetscInt idx) const {
    return shrNgbrs[idx * nNgbrVar];
  }

  /**
   * Returns ID of an element in local numeration.
   * @param idx position of element in remote neighbors list
   * @return local ID of an element
   */
  inline PetscInt local_localID(PetscInt idx) const {
    return shrNgbrs[idx * nNgbrVar + 1];
  }

  /**
   * Returns processor rank, where a "local" element is.
   * @param idx position of element in remote neighbors list
   * @return ID of the rank
   */
  inline PetscInt local_rank(PetscInt idx) const {
    return shrNgbrs[idx * nNgbrVar + 2];
  }

  /**
   * Returns information which face of a element is actually shared.
   * @param idx position of element in remote neighbors list
   * @return which face is shared
   */
  inline PetscInt local_face(PetscInt idx) const {
    return shrNgbrs[idx * nNgbrVar + 3];
  }

  /**
   * Returns ID of a neighbor element in global numeration.
   * @param idx position of element in remote neighbors list
   * @return global ID of an element in other rank
   */
  inline PetscInt remote_globalID(PetscInt idx) const {
    return shrNgbrs[idx * nNgbrVar + 4];
  }

  /**
   * Returns ID of a neigbor element in local numeration.
   * @param idx position of element in remote neighbors list
   * @return local ID of an element in other rank
   */
  inline PetscInt remote_localID(PetscInt idx) const {
    return shrNgbrs[idx * nNgbrVar + 5];
  }

  /**
   * Returns processor rank, where a "remote" (neighbor) element is.
   * @param idx position of element in remote neighbors list
   * @return ID of a rank, where the neighbor is
   */
  inline PetscInt remote_rank(PetscInt idx) const {
    return shrNgbrs[idx * nNgbrVar + 6];
  }

  /**
   * Returns information which face the remote neighbor is actually sharing.
   * @param idx position of element in remote neighbors list
   * @return which face the neighbor from other rank is sharing
   */
  inline PetscInt remote_face(PetscInt idx) const {
    return shrNgbrs[idx * nNgbrVar + 7];
  }

  /**
   * Returns number of remote neighbours.
   * @return number of remote neighbours
   */
  inline PetscInt ngbrno() const {
    return nShrNgbrs;
  }

  CRemoteNeighbor()
      : nShrNgbrs(0),
        shrNgbrs(NULL) {
  }
  ~CRemoteNeighbor() {
    free();
  }

 private:
  PetscInt nShrNgbrs;  ///< number of elements' neighbors that are in
                       ///< other domain
  PetscInt *shrNgbrs;  ///< access to raw data. See `CRemoteNeighbor::data()`
};

}  // namespace TALYFEMLIB

#endif  // GRID_CREMOTENEIGHBOR_H_
