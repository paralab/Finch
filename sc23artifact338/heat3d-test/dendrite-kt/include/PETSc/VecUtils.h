//
// Created by maksbh on 6/4/21.
//

#ifndef DENDRITEKT_VECUTILS_H
#define DENDRITEKT_VECUTILS_H

#include <petscvec.h>

struct VecUtils {
public:
  const int ndof;
  Vec & petscVec;
  const PetscScalar *val;

  inline const int getDof() const {
    return ndof;
  }

  VecUtils(Vec &v, int _ndof)
    : petscVec(v), ndof(_ndof) {
  }
};

#endif //DENDRITEKT_VECUTILS_H
