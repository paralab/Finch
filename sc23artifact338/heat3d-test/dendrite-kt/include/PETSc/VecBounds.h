//
// Created by maksbh on 12/14/20.
//

#ifndef DENDRITEKT_VECBOUNDS_H
#define DENDRITEKT_VECBOUNDS_H

#include <DataTypes.h>

/**
 * @brief Returns minimum and maximum values of a vector.
 * Must be called only from active comm.
 */

static constexpr PetscScalar MAX_VEC_VALUE = std::numeric_limits<PetscScalar>::infinity();
static constexpr PetscScalar MIN_VEC_VALUE = -std::numeric_limits<PetscScalar>::infinity();

template<DENDRITE_UINT ndof>
class VecBounds {
 public:

  static void getMaxAndMinimumValues(const Vec & v, const MPI_Comm & comm, PetscScalar *getMax, PetscScalar *getMin);
  static void setMaxAndMinimumValues(Vec & v,  const MPI_Comm & comm, const PetscScalar *setMax = nullptr,const PetscScalar *setMin = nullptr);
};

template<DENDRITE_UINT ndof>
void VecBounds<ndof>::getMaxAndMinimumValues(const Vec &v,const MPI_Comm & comm, PetscScalar *getMax, PetscScalar *getMin) {
  DENDRITE_REAL min_[ndof];
  DENDRITE_REAL max_[ndof];

  std::fill(std::begin(min_),std::end(min_),MAX_VEC_VALUE);
  std::fill(std::begin(max_),std::end(max_),MIN_VEC_VALUE);

  const PetscScalar *array;
  VecGetArrayRead(v, &array);
  PetscInt localSize;
  VecGetLocalSize(v, &localSize);
  for (DENDRITE_UINT i = 0; i < localSize; i += ndof) {
    for (DENDRITE_UINT dof = 0; dof < ndof; dof++) {
      if (min_[dof] > array[i + dof]) {
        min_[dof] = array[i + dof];
      }
      if (max_[dof] < array[i + dof]) {
        max_[dof] = array[i + dof];
      }
    }
  }
  VecRestoreArrayRead(v, &array);
  MPI_Allreduce(max_, getMax, ndof, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(min_, getMin, ndof, MPI_DOUBLE, MPI_MIN, comm);

}


template<DENDRITE_UINT ndof>
void VecBounds<ndof>::setMaxAndMinimumValues(Vec &v, const MPI_Comm & comm, const PetscScalar *setMax, const PetscScalar *setMin) {
  DENDRITE_REAL min_[ndof];
  DENDRITE_REAL max_[ndof];
  std::fill(std::begin(min_),std::end(min_),MIN_VEC_VALUE);
  std::fill(std::begin(max_),std::end(max_),MAX_VEC_VALUE);
  if(setMax != nullptr){
    std::memcpy(max_,setMax,sizeof(PetscScalar)*ndof);
  }
  if(setMin != nullptr){
    std::memcpy(min_,setMin,sizeof(PetscScalar)*ndof);
  }
  PetscScalar *localArray;
  VecGetArray(v, &localArray);
  PetscInt localSize;
  VecGetLocalSize(v, &localSize);
  for (DENDRITE_UINT i = 0; i < localSize; i += ndof) {
    for (DENDRITE_UINT dof = 0; dof < ndof; dof++) {
      if (localArray[i + dof] < min_[dof]) {
        localArray[i + dof] = min_[dof];
      }
      if (localArray[i + dof] > max_[dof]) {
        localArray[i + dof] = max_[dof];
      }
    }
  }
  VecRestoreArray(v, &localArray);
}

#endif //DENDRITEKT_VECBOUNDS_H
