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
#include <talyfem/utils/test_utils.h>

#include <string>

namespace TALYFEMLIB {

void PrintResults(std::string err_string, bool pass, bool should_fail) {
  if (GetMPIRank() == 0) {  // only print from rank 0
    const char* pass_color = "\033[92m";
    const char* fail_color = "\033[91m";
    const char* end_color = "\033[0m";
    if (pass) {
      std::cout << (should_fail ? fail_color : pass_color) << "[TEST-PASSED]"
          << (should_fail ? " (should fail) " : " ") << end_color
          << err_string << std::endl;
    } else {
      std::cout << (should_fail ? pass_color : fail_color) << "[TEST-FAILED]"
          << (should_fail ? " (should fail) " : " ") << end_color
          << err_string << std::endl;
    }
  }
}

double calc_value_error(double expected, double actual) {
  // for the calc_value_error(0.0, 0.0) case
  if (expected == actual) {
    return 0.0;
  }

  return fabs(expected - actual) / fabs(expected);
}

double calc_value_error(const ZEROPTV& expected, const ZEROPTV& actual) {
  return calc_value_error(expected.x(), actual.x()) +
         calc_value_error(expected.y(), actual.y()) +
         calc_value_error(expected.z(), actual.z());
}

double calc_value_error(int expected, int actual) {
  // If the values are the same, the error is zero. We do this check to avoid
  // doing floating point calculations unless needed
  if (expected == actual) {
    return 0.0;
  }

  const double double_expected = static_cast<double>(expected);
  const double double_actual = static_cast<double>(actual);
  return calc_value_error(double_expected, double_actual);
}

}  // namespace TALYFEMLIB
