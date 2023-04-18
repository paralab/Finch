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
#ifndef UTILS_UTILS_H_
#define UTILS_UTILS_H_

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/mpi.h>
#import <PETSc/petsc.h>
#else
#include <mpi.h>
#include <petsc.h>
#endif

#include <iostream>
#include <sstream>
#include <string>

#include <talyfem/utils/macros.h>
#include <talyfem/utils/reproducibility.h>

// Used to mark a function as deprecated. Deprecated functions are obsolete and
// should not be used. They are still included in the library to avoid breaking
// existing code, but should not be used in new code. They may be removed in
// the future.
// Marking a function as deprecated will cause the compiler to issue a
// warning during compilation if that function is used in a code. The warnings
// can be fixed by changing the code to use the new function that replaces
// the deprecated one. When marking a function deprecated, there should be a
// comment above the function stating what to use instead.
// Example of marking function deprecated:
//   not deprecated:
//     void DoSomething(int count, double value, int err) { ... }
//
//   deprecated:
//     // Deprecated in v6.9.7, instead use:
//     // int DoSomething(int count, double value)
//     void DEPRECATED(DoSomething(int count, double value, int err)) { ... }
//
// #ifdef __GNUC__
// #define DEPRECATED(func) func __attribute__ ((deprecated))
// #elif defined(_MSC_VER)
// #define DEPRECATED(func) __declspec(deprecated) func
// #elif __INTEL_COMPILER
// #define DEPRECATED(func) func __attribute__ ((deprecated))
// #else
// #pragma message("***WARNING: need to implement DEPRECATED for this compiler")
// #define DEPRECATED(func) func
// #endif

namespace TALYFEMLIB {

// Store the flags for whether to print types of output. These are global
// variables kept in a separate namespace to avoid polluting the global
// namespace. If a value is set to true, that type of output will be done.
// These can be changed at any point in the program, but the intent is that
// they will be set in the InputData class based on input values.
namespace GLOBALS {
extern bool gPrintStat;
extern bool gPrintInfo;
extern bool gPrintWarn;
extern bool gPrintLog;
extern bool gPrintTime;
extern bool gPrintEOL;  // whether to print the end of line character
extern int gRankOfPrinter;
}

/**
 * Set which process can print using the logging functions (PrintInfo, etc).
 * @param new_rank rank of the new printer
 */
void set_RankOfPrinter(const int new_rank);

/**
 * Set whether or not the logging functions automatically end with a newline.
 * @param new_value if logging functions should automatically print a newline
 */
void set_gPrintEOL(const bool new_value);

/**
 * Print the "\r" character to std::cout.
 * This "rewinds" the print cursor to the start of the current line.
 */
void PrintStatusRewind();

/**
 * Returns the MPI rank of this process
 *
 * @param is_serial a bool specifying whether process is serial
 * @return MPI rank of this process
 */
static inline int GetMPIRank(bool is_serial) {
  int rank;
  if (is_serial) {
    MPI_Comm_rank(PETSC_COMM_SELF, &rank);
  } else {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  }
  return rank;
}

/**
 * Returns the total number of MPI processes
 *
 * @param is_serial a bool specifying whether process is serial
 * @return number of MPI processes
 */
static inline int GetMPISize(bool is_serial) {
  int size;
  if (is_serial) {
    MPI_Comm_size(PETSC_COMM_SELF, &size);
  } else {
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
  }
  return size;
}

/**
 * Returns the MPI rank of this process
 *
 * @return MPI rank of this process
 */
static inline int GetMPIRank() {
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  return rank;
}

/**
 * Returns the total number of MPI processes
 *
 * @return number of MPI processes
 */
static inline int GetMPISize() {
  int size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  return size;
}


// Include variadic print functions.
#include <talyfem/utils/print_utils_variadic.h>


/**
 * Converts a value to a std::string, returning the result.
 *
 * @param val Value which will be converted to a string
 * @returns val as a string
 */
template<typename T>
static std::string ToString(const T& val) {
  std::stringstream out;
  out << val;
  return out.str();
}

/**
 * This is used to count the number of bits set to 1 in an integer.
 *
 * This is actually a well-documented problem; it is often called the
 * "Hamming Weight" or popcount
 * Check out this wikipedia article for more info:
 * http://en.wikipedia.org/wiki/Hamming_weight
 * This particular implementation is taken from here:
 * http://bisqwit.iki.fi/source/misc/bitcounting/ (WP3 - Nifty Revised)
 *
 * @param n integer to count the number of bits set of
 * @returns the number of bits set in n
 */
template<typename T>
int popcount(T n) {
  static const unsigned int TEST_BITS = sizeof(T) * 8;

  static const T m1 = (~(T) 0) / 3;
  static const T m2 = (~(T) 0) / 5;
  static const T m4 = (~(T) 0) / 17;
  static const T h01 = (~(T) 0) / 255;

  n -= (n >> 1) & m1;
  n = (n & m2) + ((n >> 2) & m2);
  n = (n + (n >> 4)) & m4;

  return (n * h01) >> (TEST_BITS - 8);
}

/**
 * Prints a simple ASCII bar chart (to PrintStatus) for an array of values.
 * @param values values to print
 * @param count number of values
 * @param label label for the chart
 */
template <typename T>
void PrintBarChart(T* values, unsigned int count, const std::string& label) {
  T max = 0;
  unsigned int max_idx = 0;
  for (unsigned int i = 0; i < count; i++) {
    if (values[i] > max) {
      max = values[i];
      max_idx = i;
    }
  }

  std::stringstream ss;
  ss << label << " (chart):\n\n";

  const int lines = 32;
  const long double scale = static_cast<long double>(max) / lines;
  for (int line = lines; line > 0; line--) {
    ss << "| ";
    for (unsigned int i = 0; i < count; i++) {
      if (values[i] >= scale * line)
        ss << (i == max_idx ? "@" : "#");
      else
        ss << "-";
    }
    ss << "\n";
  }
  for (unsigned int i = 0; i < count + 3; i++) {
    ss << "=";
  }
  ss << "\n\n";
  ss << "(maximum: " << max << ", at " << max_idx << ")";

  PrintStatus(ss.str());
}

// compatibility layer for PETSc 3.7's new PetscOptions* functions (maps to the 3.6- versions)
// used to keep the tests that use PetscOptionsGet* compatible with both 3.6 and 3.7+
// can be removed once we can fully deprecate 3.6 (in however many years...)
// intentionally in the TALYFEMLIB namespace to make it clear these are not the "real" petsc functions
#if PETSC_VERSION_LT(3, 7, 0)
inline PetscErrorCode PetscOptionsHasName(void*,const char pre[],const char name[],PetscBool  *set) { return PetscOptionsHasName(pre, name, set); }
inline PetscErrorCode PetscOptionsGetInt(void*,const char pre[],const char name[],PetscInt *ivalue,PetscBool  *set) { return PetscOptionsGetInt(pre, name, ivalue, set); }
inline PetscErrorCode PetscOptionsGetBool(void*,const char pre[],const char name[],PetscBool  *ivalue,PetscBool  *set) { return PetscOptionsGetBool(pre, name, ivalue, set); }
inline PetscErrorCode PetscOptionsGetScalar(void*,const char pre[],const char name[],PetscScalar *dvalue,PetscBool  *set) { return PetscOptionsGetScalar(pre, name, dvalue, set); }
inline PetscErrorCode PetscOptionsGetString(void*,const char pre[],const char name[],char str[],size_t len,PetscBool  *set) { return PetscOptionsGetString(pre, name, str, len, set); }
#define PetscOptionsSetValue(dummy, iname, value) PetscOptionsSetValue(iname, value)
#define PetscOptionsClearValue(dummy, iname) PetscOptionsClearValue(iname)
#endif

}  // namespace TALYFEMLIB

#endif  // UTILS_UTILS_H_
