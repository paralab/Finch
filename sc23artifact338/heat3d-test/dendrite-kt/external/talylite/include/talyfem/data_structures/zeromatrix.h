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
#ifndef DATASTRUCTURE_ZEROMATRIX_H_
#define DATASTRUCTURE_ZEROMATRIX_H_

#include <assert.h>

#include <cstdio>  // for std::printf

#include <talyfem/data_structures/zeroarray.h>
#include <talyfem/common/exceptions.h>

namespace TALYFEMLIB {
/**
 * Class of matrix data structure
 */
template<class T>
class ZEROMATRIX {
 public:
  ZEROMATRIX()
      : nx_(0),
        ny_(0),
        data_(NULL) {
  }

  /**
   * Creates a matrix with the given dimensions.
   *
   * @param new_nx the number of items in the x direction
   * @param new_ny the number of items in the y direction
   */
  ZEROMATRIX(int new_nx, int new_ny)
      : nx_(0),
        ny_(0),
        data_(NULL) {
    redim(new_nx, new_ny);
  }

  virtual ~ZEROMATRIX() {
    cleanup();
  }

  /**
   * Clears all data from the object
   */
  void cleanup() {
    if (data_) {
      delete[] data_;
      data_ = NULL;
    }
  }

  /**
   * Resizes the object to the given size
   *
   * If the new dimensions are the same as the current dimensions, nothing is
   * done.
   *
   * @param new_nx The new x dimension of the object
   * @param new_ny The new y dimension of the object
   */
  virtual void redim(int new_nx, int new_ny) {
    if (new_nx == nx_ && new_ny == ny_) {
      // size remains the same, don't do anything
      return;
    }
    cleanup();
    nx_ = new_nx;
    ny_ = new_ny;
    data_ = new T[nx_ * ny_];
  }

  /**
   * Accesses the item at the given index
   *
   * @param i the i index of item
   * @param j the j index of item
   * @return reference to item in the array
   */
  T& operator () (int i, int j) {
    assert(i >= 0 && i < nx_);
    assert(j >= 0 && j < ny_);
    return data_[i * ny_ + j];
  }

  /**
   * Same thing, just optionally const.
   * @param i the i index of the item
   * @param j the j index of the item
   * @return const reference to item in the array
   */
  const T& operator () (int i, int j) const {
    assert(i >= 0 && i < nx_);
    assert(j >= 0 && j < ny_);
    return data_[i * ny_ + j];
  }

  /**
   * Assignment operator
   *
   * @param A Matrix to assign to the object
   * @return object after assignment
   */
  ZEROMATRIX& operator=(const ZEROMATRIX& A) {
    redim(A.nx_, A.ny_);
    for (int i = 0; i < nx_ * ny_; i++) {
      data_[i] = A.data_[i];
    }
    return *this;
  }

  /**
   * Returns the raw data in this array.
   *
   * @return Pointer to the data in this array.
   */
  inline const T* data() const {
    return data_;
  }

  /**
   * Returns the raw data in this array.
   *
   * @return Pointer to the data in this array.
   */
  inline T* data_ptr() {
    return data_;
  }

  /**
   * Returns the size of the matrix in the x direction
   *
   * @return size of the matrix in the x direction
   */
  inline int nx() const {
    return nx_;
  }

  /**
   * Returns the size of the matrix in the y direction
   *
   * @return size of the matrix in the y direction
   */
  inline int ny() const {
    return ny_;
  }

 protected:
  int nx_;  ///< size of the matrix in the x direction
  int ny_;  ///< size of the matrix in the y direction
  T* data_;  ///< array to store the matrix data
};

/**
 * Class of matrix data structure (including operations, solving Ax=b)
 */
template<class T>
class ZeroMatrix : public ZEROMATRIX<T> {
 public:
  using ZEROMATRIX<T>::nx_;
  using ZEROMATRIX<T>::ny_;
  using ZEROMATRIX<T>::data_;
  using ZEROMATRIX<T>::cleanup;
  using ZEROMATRIX<T>::redim;
  using ZEROMATRIX<T>::operator();

  virtual ~ZeroMatrix() { }

  /**
   * Sets every value in the matrix to data.
   *
   * @param value Set every value to this.
   */
  void fill(const T& value) {
    for (int i = 0; i < nx_ * ny_; i++) {
      data_[i] = value;
    }
  }

  /**
   * Sets this Matrix as the product of the given matrices
   *
   * The matrices A and B are multiplied, after optionally being transposed.
   * The resulting matrix is stored in the object and the object is returned.
   *
   * @param A first matrix to multiply
   * @param B seconds matrix to multiply
   * @param transposeA whether to transpose A in multiplication
   * @param transposeB whether to transpose B in multiplication
   * @return reference to this object
   */
  ZeroMatrix& prod(const ZeroMatrix& A, const ZeroMatrix& B,
                   bool transposeA = false, bool transposeB = false) {
    int I, J, J1, K;

    if (!transposeA) {
      I = A.nx();
      J = A.ny();
    } else {
      I = A.ny();
      J = A.nx();
    }

    if (!transposeB) {
      J1 = B.nx();
      K = B.ny();
    } else {
      J1 = B.ny();
      K = B.nx();
    }

    if (J != J1) {
      throw TALYException() << "Error multiplying matrices (I: " << I
          << ", J: " << J << ", J1: " << J1 << ", K: " << K << ")";
    }

    redim(I, K);
    for (int i = 0; i < I; i++) {
      for (int k = 0; k < K; k++) {
        (*this)(i, k) = 0;
        for (int j = 0; j < J; j++) {
          const T& a_ij = transposeA ? A(j, i) : A(i, j);
          const T& b_jk = transposeB ? B(k, j) : B(j, k);
          (*this)(i, k) += a_ij * b_jk;
        }
      }
    }
    return *this;
  }

  /**
   * Scales the items in the matrix by the given value
   *
   * @param value the value by which to scale the matrix items
   * @return a reference to object after scaling
   */
  ZeroMatrix& ratio(const T& value) {
    for (int i = 0; i < nx_ * ny_; i++) {
      data_[i] *= value;
    }
    return *this;
  }

  /**
   * Adds a matrix to the current object
   *
   * Mathematically, this operation is (for all i,j): x(i,j) = x(i,j) + A(i,j)
   *
   * @param A matrix to add to the object
   * @return a reference to object after addition
   */
  ZeroMatrix& plus(const ZeroMatrix& A) {
    assert(A.nx_ == nx_ && A.ny_ == ny_);
    for (int i = 0; i < nx_ * ny_; i++) {
      data_[i] += A.data_[i];
    }
    return *this;
  }

  /**
   * Adds a constant to all items in the matrix
   *
   * Mathematically, this operation is (for all i,j): x(i,j) = x(i,j) + C
   *
   * @param constant value to add to each item in the matrix
   * @return a reference to object after addition
   */
  ZeroMatrix& plus(const T& constant) {
    for (int i = 0; i < nx_ * ny_; i++) {
      data_[i] += constant;
    }
    return *this;
  }

  /**
   * Transposes the matrix, returning a reference to the object
   *
   * @return reference to the object
   */
  ZeroMatrix& transpose() {
    ZeroMatrix tmpM;
    tmpM.redim(ny_, nx_);
    for (int i = 0; i < nx_; i++) {
      for (int j = 0; j < ny_; j++) {
        tmpM(j, i) = (*this)(i, j);
      }
    }
    *this = tmpM;
    return *this;
  }

  /**
   * ???
   *
   * @param trialStress_dev ???
   * @return ???
   */
  ZeroMatrix& getDevatoricPart(ZeroMatrix& trialStress_dev) const {
    trialStress_dev.redim(nx_, ny_);
    const T self_trace = this->trace();
    for (int i = 0; i < nx_; i++) {
      for (int j = 0; j < ny_; j++) {
        trialStress_dev(i, j) = (*this)(i, j);
      }
      trialStress_dev(i, i) -= self_trace / nx_;
    }
    return trialStress_dev;
  }

  /**
   * Returns the trace of the matrix
   *
   * @return the trace of the matrix
   */
  const T trace() const {
    T sum = 0;
    for (int i = 0; i < nx_; i++) {
      sum += (*this)(i, i);
    }
    return sum;
  }

  /**
   * Returns the sum of the squares of every item in the matrix
   *
   * @return sum of the squares of every item in the matrix
   */
  T square() const {
    T sum = 0;
    for (int i = 0; i < nx_; i++) {
      for (int j = 0; j < ny_; j++) {
        sum += (*this)(i, j) * (*this)(i, j);
      }
    }
    return sum;
  }

  /**
   * Returns the root of the sum of the squares of every item in the matrix
   *
   * @return square root sum of the squares of every item in the matrix
   */
  T sqr() const {
    return sqrt(square());
  }

  /**
   * Sets this matrix to the identity matrix
   */
  void setIdentity() {
    fill(0);
    for (int i = 0; i < nx_; i++) {
      (*this)(i, i) = 1;
    }
  }

  /**
   * ??? - I don't think this has been tested
   */
  void factLU() {
    int n = nx_;

    ZeroMatrix& L = *this;
    ZeroMatrix& A = *this;

    for (int k = 0; k < n - 1; k++) {
      for (int i = k + 1; i < n; i++) {
        L(i, k) = A(i, k) / A(k, k);
      }
      // calculate         (*)
      //                   (*)
      //        L'*U'=A' - (*){* * * * *}
      //                   (*)
      //                   (*)
      for (int j = k + 1; j < n; j++) {
        for (int i = k + 1; i < n; i++) {
          A(i, j) -= L(i, k) * A(k, j);
        }
      }
    }
  }

  /**
   * ??? - I don't think this has been tested
   *
   * @param b ???
   * @param x ???
   */
  void forwBackLU(const ZEROARRAY<T>& b, ZEROARRAY<T>& x) {
    ZeroMatrix& L = *this;
    ZeroMatrix& U = *this;
    ZEROARRAY<T>& y = x;

    int n = nx_;

    for (int i = 0; i < n; i++) {
      T temp = 0.;
      for (int j = 0; j < i; j++) {
        temp += L(i, j) * y(j);
      }
      y(i) = (b(i) - temp);
    }

    for (int i = n; i >= 0; i--) {
      T temp = 0.;
      for (int j = i + 1; j < n; j++) {
        temp += U(i, j) * x(j);
      }
      x(i) = (y(i) - temp) / U(i, i);
    }
  }

  /**
   * Calculates the determinant of the matrix
   *
   * @return the determinant of the matrix
   * @note: I don't think this has been tested
   */
  double det() const {
    if (nx_ == 1)
      return (*this)(0, 0);

    double sum = 0;
    for (int i = 0; i < nx_; i++) {
      int sgn = (i % 2 == 1) ? +1 : -1;

      double val = (*this)(0, i);
      ZeroMatrix<T> son;
      son.redim(nx_ - 1, nx_ - 1);
      for (int a = 0; a < i - 1; a++) {
        for (int b = 0; b < nx_ - 1; b++) {
          son(b, a) = (*this)(b + 1, a);
        }
      }
      for (int a = i; a < nx_ - 1; a++) {
        for (int b = 0; b < nx_ - 1; b++) {
          son(b, a) = (*this)(b + 1, a + 1);
        }
      }
      sum += sgn * val * son.det();
    }
    return sum;
  }
};

}  // namespace TALYFEMLIB

#endif  // DATASTRUCTURE_ZEROMATRIX_H_
