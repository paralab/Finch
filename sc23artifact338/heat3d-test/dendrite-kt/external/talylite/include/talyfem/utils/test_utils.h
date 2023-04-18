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
#ifndef UTILS_TEST_UTILS_H_
#define UTILS_TEST_UTILS_H_

#include <iostream>  // std::cout, std::endl
#include <string>

#include <talyfem/grid/zeroptv.h>
#include <talyfem/utils/utils.h>  // for GetMPIRank()


namespace TALYFEMLIB {

/**
 * Prints the results (pass/fail) of this program run
 *
 * This is used when using the code as a test for the library. Code that passes
 * will normally print "[TEST-PASSED]" is green and code that fails will
 * normally print "[TEST-FAILED]" in red. The colors can be reversed by setting
 * should_fail to true. In this case, we are expecting the test to fail so we
 * get [TEST-PASSED] is red to indicate a problem and [TEST-FAILED] in green
 * to indicate that it is the expected result.
 *
 * @param err_string print that is printed after the test results (not colored)
 * @param pass whether the test passed
 * @param should_fail whether the test should fail
 */
void PrintResults(std::string err_string, bool pass, bool should_fail);

/**
 * Returns the relative error (% error / 100.0) between the two values.
 *
 * If the values are the same (even zero), the error is zero. Otherwise, the
 * error is defined relative to the expected value.
 * Note: this the value is undefined if expected == 0 && actual != 0
 *
 * @param expected the expected value that the error is based on
 * @param actual the actual value we have
 * @return the error value
 */
double calc_value_error(double expected, double actual);

/**
 * Returns the relative error (% error / 100.0) between two points.
 *
 * This is the sum of the errors between the x, y, and z pairs.
 *
 * @param expected the point value that the error is based on
 * @param actual the actual point we have
 * @return the error value
 */
double calc_value_error(const ZEROPTV& expected, const ZEROPTV& actual);

/**
 * Returns the relative error (% error / 100.0) between two integers.
 *
 * If the values are the same, 0.0 is returned. Otherwise, the error is
 * determined by converting the values to doubles and calculating the error of
 * those values.
 *
 * @param expected the expected integer value that the error is based on
 * @param actual the actual integer we have
 * @return the error value
 */
double calc_value_error(int expected, int actual);

}  // namespace TALYFEMLIB

#endif  // UTILS_TEST_UTILS_H_
