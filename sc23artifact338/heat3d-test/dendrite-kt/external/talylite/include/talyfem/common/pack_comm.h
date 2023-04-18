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
#ifndef COMMON_PACKCOMM_H_
#define COMMON_PACKCOMM_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace TALYFEMLIB {

/**
 * Get a line from the file.
 *
 * @param[in] fp a file pointer
 * @param[out] str contents of a line from the file
 */
void getLine(FILE* fp, char* str);

/**
 * Get value of a parameter from file.
 *
 * @param[in] fp a file pointer
 * @param[in] varname name of the parameter
 * @param[out] var value of the parameter
 * @param divideChar character string separating name and value
 * @return the parameter value (same as var)
 */
char* getParameter(FILE* fp, const char* varname, char* var,
                   const char* divideChar = "=");

/**
 * trim the spaces of a string (both sides).
 *
 * @param[in] str string to be trimed
 * @param[out] s string after trimes
 */
void trimSpace(char* str, char* s);

/**
 * divide a string into several substrings
 *
 * @param[in] str origin string
 * @param[out] position positions of substrings
 * @param[in] divideChar characters to divide string
 * @return number of substrings
 */
int divideParameters(char* str, char* position[], const char* divideChar);

}  // namespace TALYFEMLIB

#endif  // COMMON_PACKCOMM_H_
