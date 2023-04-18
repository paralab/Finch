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
#ifndef DATASTRUCTURE_ZEROARRAY_H_
#define DATASTRUCTURE_ZEROARRAY_H_

#include <string.h>  // for memcpy/memset
#include <assert.h>

#include <iostream>  // for std::cout

namespace TALYFEMLIB {

/**
 * A zero-indexed array class.
 */
template<class Data>
class ZEROARRAY {
 public:
  ZEROARRAY()
      : size_(0),
        data_(NULL) {
  }

  virtual ~ZEROARRAY() {
    cleanup();
  }

  /**
   * Clears all data from the object
   */
  virtual void cleanup() {
    if (data_) {
      delete[] data_;
      data_ = NULL;
      size_ = 0;
    }
  }

  /**
   * Copy constructor.
   *
   * @param A Array to copy from.
   */
  ZEROARRAY(const ZEROARRAY& A)
    : ZEROARRAY() {
    redim(A.size_);
    for (int i = 0; i < size_; i++) {
      data_[i] = A.data_[i];
    }
  }

  /**
   * Assignment operator
   *
   * @param A Array to assign to the object
   * @return object after assignment
   */
  ZEROARRAY& operator=(const ZEROARRAY& A) {
    // if we're doing "A = A;" do nothing
    if (this == &A)
      return *this;

    redim(A.size_);
    for (int i = 0; i < size_; i++) {
      data_[i] = A.data_[i];
    }
    return *this;
  }

  /**
   * Resizes the object to the given size
   *
   * If the new dimension is the same as the current dimensions, nothing is
   * done.
   *
   * @param new_size The new size of the object
   */
  void redim(int new_size) {
    if (new_size == size_) {
      // size remains the same, don't do anything
      return;
    }
    cleanup();
    size_ = new_size;
    data_ = new Data[new_size];
  }

  /**
   * Fills the array with a numeric sequence of the given length
   *
   * This stores a sequence 0, 1, 2, 3, ... in the array. The number of items
   * is given by the argument. The starting value defaults to zero but can be
   * changed as needed. This function does not make sense for data types that
   * are not numeric.
   *
   * @param length number of terms in the sequence
   * @param start the starting term in the sequence
   */
  void fill_sequence(int length, int start = 0) {
    redim(length);
    for (int i = 0; i < size_; i++) {
      data_[i] = i + start;
    }
  }

  /**
   * Fills the array with data from a C array
   *
   * @param new_data pointer to data to copy
   * @param length number of terms in the sequence
   */
  void fill_from_array(Data *new_data, int length) {
    redim(length);
    for (int i = 0; i < size_; i++) {
      data_[i] = new_data[i];
    }
  }

  /**
   * Fills the array with data from a C array
   *
   * @param new_data pointer to data to copy
   * @param length number of terms in the sequence
   */
  void fill_from_array(const Data *new_data, int length) {
    redim(length);
    for (int i = 0; i < size_; i++) {
      data_[i] = new_data[i];
    }
  }

  /**
   * Accesses the item at the given index
   *
   * @param i index of item
   * @return reference to item in the array
   */
  inline Data& operator()(int i) {
    assert(i >= 0 && i < size_);
    return data_[i];
  }

  /**
   * Const accessor.
   * @param i index of item
   * @return a const reference to item in the array
   */
  inline const Data& operator()(int i) const {
    assert(i >= 0 && i < size_);
    return data_[i];
  }

  /**
   * Sets every value in this array to data.
   *
   * @param new_data_value Set every value to this.
   */
  void fill(const Data new_data_value) {
    for (int i = 0; i < size_; i++) {
      data_[i] = new_data_value;
    }
  }

  /**
   * Adds an array, optionally scaled by a constant, to the object
   *
   * Mathematically, this operation is (for all i): a(i) = a(i) + x * b(i)
   *
   * @param other_array array to add to the object
   * @param scale_ratio value to scale added array by
   */
  void add(const ZEROARRAY<Data>& other_array, Data scale_ratio) {
    for (int i = 0; i < size_; i++) {
      data_[i] += other_array.data_[i] * scale_ratio;
    }
  }

  /**
   * Returns the norm (square root of the sum of the squares) of the array
   *
   * @return the norm of the array
   */
  Data norm() const {
    Data sum = 0;
    for (int i = 0; i < size_; i++) {
      sum += data_[i] * data_[i];
    }
    return sqrt(sum);
  }

  /**
   * Prints the array contents to the terminal
   *
   * The contents are formatted in the specified number of columns.
   *
   * @param n_columns number of columns to use while printing
   */
  void print(int n_columns = 2) const {
    for (int i = 0; i < size_; i++) {
      std::cout << data_[i] << " ";
      if (i % n_columns == (n_columns - 1)) {
        std::cout << "\n";
      }
    }
  }

  /**
   * Adds the given value to the end of the array.
   *
   * This will extend the length of the data array by one.
   *
   * @param value the value to append to the array
   * @note: do not use this function frequently, its efficiency is bad
   */
  void appendData(const Data value) {
    Data* new_data = new Data[size_ + 1];
    for (int i = 0; i < size_; i++) {
      new_data[i] = data_[i];
    }
    new_data[size_] = value;

    delete[] data_;
    data_ = new_data;
    size_++;
  }

  /**
   * Checks if given value is in the array
   *
   * @param value the value to search for
   * @return true if the value is in the array
   */
  bool contains(const Data value) const {
    for (int i = 0; i < size_; i++) {
      if (data_[i] == value) {
        return true;
      }
    }
    return false;
  }

  /**
   * Returns the number of elements in this array (0 if uninitialized)
   *
   * @return The number of elements in this array.
   */
  inline int size() const {
    return size_;
  }

  /**
   * Returns the raw data in this array.
   *
   * @return Pointer to the data in this array.
   */
  inline Data* data() const {
    return data_;
  }

  /**
   * Returns index of first instance of the given item or -1 if it is not found.
   *
   * @param value The data item to find
   * @return index of first instance of item in array or -1 if not found.
   */
  int find(const Data value) const {
    for (int i = 0; i < size_; i++) {
      if (data_[i] == value) {
        return i;
      }
    }
    return -1;
  }

  /**
   * Scales the items in the array by the given value
   *
   * @param value the value by which to scale the array items
   */
  void ratio(const Data value) {
    for (int i = 0; i < size_; i++) {
      Data& v = data_[i];
      v = v * value;
    }
  }

 private:
  int size_;  ///< number of data items in the array
  Data* data_;  ///< array of data
};

/**
 * Boolean specialization for ARRAY class.
 *
 * Represents each value using a single bit to save memory.
 */
template<>
class ZEROARRAY<bool> {
 public:
  ZEROARRAY() {
    data_ = NULL;
    size_ = 0;
  }

  virtual ~ZEROARRAY() {
    cleanup();
  }

  /**
   * Copy constructor.
   *
   * @param A Array to copy from.
   */
  ZEROARRAY(const ZEROARRAY& A) : ZEROARRAY() {
    redim(A.size_);
    memcpy(data_, A.data_, storage_size(size_) * sizeof(BoolStorage));
  }

  /**
   * Assignment operator
   *
   * @param A Array to assign to the object
   * @return object after assignment
   */
  ZEROARRAY& operator=(const ZEROARRAY& A) {
    // if we're doing "A = A;" do nothing
    if (this == &A)
      return *this;

    redim(A.size_);
    memcpy(data_, A.data_, storage_size(size_) * sizeof(BoolStorage));
    return *this;
  }

  /**
   * Frees the memory used by the object.
   */
  void cleanup() {
    if (data_) {
      delete[] data_;
      data_ = NULL;
      size_ = 0;
    }
  }

  /**
   * Resizes the object
   *
   * @param new_size The new size of the object
   */
  void redim(int new_size) {
    cleanup();
    size_ = new_size;
    data_ = new BoolStorage[storage_size(new_size)];
  }

  /**
   * Accesses the item at the given index
   *
   * @param i index of item
   * @return reference to item in the array
   */
  inline bool get(int i) const {
    assert(i >= 0 && i < size_);
    return data_[i / sizeof(BoolStorage)] & (1 << (i % sizeof(BoolStorage)));
  }

  /**
   * Accesses the item at the given index
   *
   * @param i index of item
   * @param value value to set at this index
   * @return reference to item in the array
   */
  inline void set(int i, bool value) {
    assert(i >= 0 && i < size_);

    // An example. sizeof(BoolStorage) == 16. size_ == 18.
    // Since ceil(size_ / sizeof(BoolStorage)) == 2, data_ has 2 elements.
    // Initially, only the item at i == 14 is true.
    // |0100000000000000|0000000000000000|

    // set(2, true):
    // 2 / sizeof(BoolStorage) == 0, so we use the 0th element in data_.
    // 2 % sizeof(BoolStorage) == 2, so we shift left by 2 in set/get.
    //         |0100000000000000|0000000000000000|
    // | (OR)  |0000000000000100| (1 << 2)
    //         -----------------------------------
    //         |0100000000000100|0000000000000000|

    // get(2):
    //         |0100000000000100|0000000000000000|
    // & (AND) |0000000000000100| (1 << 2)
    //         -----------------------------------
    //         |0000000000000100| (which is nonzero, so C interprets as 'true')

    // set(2, false):
    //         |0100000000000100|0000000000000000|
    // & (AND) |1111111111111011| ~(1 << 2)
    //         -----------------------------------
    //         |0100000000000000|0000000000000000|

    if (value) {
      data_[i / sizeof(BoolStorage)] |= (1 << (i % sizeof(BoolStorage)));
    } else {
      data_[i / sizeof(BoolStorage)] &= ~(1 << (i % sizeof(BoolStorage)));
    }
  }

  /**
   * Returns the number of elements in this array (0 if uninitialized)
   *
   * @return The number of elements in this array.
   */
  inline int size() const {
    return size_;
  }

  /**
   * Sets every value in this array to data.
   *
   * @param value Set every value to this.
   */
  void fill(bool value) {
    memset(data_, value ? ~0 : 0, storage_size(size_) * sizeof(BoolStorage));
  }

 private:
  typedef unsigned int BoolStorage;

  /**
   * Returns the smallest number of BoolStorages
   * that can hold bool_count booleans.
   * @param bool_count :  Number of boolean values to store.
   * @return The number of BoolStorages needed to store bool_count bools.
   */
  inline static size_t storage_size(int bool_count) {
    return (bool_count / sizeof(BoolStorage)) +
            (bool_count % sizeof(BoolStorage) == 0 ? 0 : 1);
  }

  int size_;
  BoolStorage* data_;
};


}  // namespace TALYFEMLIB

#endif  // DATASTRUCTURE_ZEROARRAY_H_
