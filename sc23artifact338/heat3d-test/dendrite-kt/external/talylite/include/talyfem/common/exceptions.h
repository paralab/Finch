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
#import <PETSc/petscsys.h>
#else
#include <petscsys.h>
#endif

#include <string>
#include <exception>
#include <sstream>

namespace TALYFEMLIB {

/**
 * Base exception class for all TalyFEM errors.
 *
 * You can write messages C++ stream-style like so:
 *
 * if (error_condition) {
 *   throw TALYException() << "Bad things happened! More information: " << error_info;
 * }
 *
 * Use this to include relevant information in your error messages!
 *
 * If you want to throw a new type of exception that users should be
 * able to catch, please subclass this.
 */
class TALYException : public std::exception {
 public:
  virtual ~TALYException() throw() {}

  /**
   * Return the message associated with this exception, if any.
   */
  virtual const char* what() const throw() { return msg_.c_str(); }

  /**
   * Stream input for message.
   * @param e Exception object (left-hand side of <<)
   * @param obj Value to append to message (right-hand side of <<)
   * @returns TALYException with new error message
   */
  template<typename T>
  friend TALYException operator<<(TALYException e, const T& obj);

  /**
   * Print this exception's message using TalyFEM's PrintError function.
   */
  void print() const;

 protected:
  std::string msg_;  ///< The message string
};

template<typename T>
TALYException operator<<(TALYException e, const T& obj) {
  std::stringstream ss;
  ss << e.msg_ << obj;
  e.msg_ = ss.str();
  return e;
}

/**
 * Exception thrown when an option as not been implemented.
 */
class NotImplementedException : public TALYException {};

/**
 * Exception thrown when a file was not loaded successfully.
 */
class FileIOException : public TALYException {};

/**
 * Exception thrown when a PETSc error occurs and PetscExceptionErrorHandler
 * is set as the current PETSc error handler.
 * Typically used in cases where PETSc errors should abort the program
 * (for example, when there is an error during assembly).
 */
class PetscException : public TALYException {
 public:
  PetscException();
};

/**
 * PETSc error handler that throws an exception instead of just printing
 * an error message and continuing with the program.
 */
PetscErrorCode PetscExceptionErrorHandler(MPI_Comm comm, int line,
    const char* func, const char* file, PetscErrorCode code,
    PetscErrorType type, const char* msg, void* ctx);

}  // namespace TALYFEMLIB

