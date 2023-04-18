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

#include <stabilizertest.h>
#include <globals.h>


#define TRY_RECORD(sec, subsec, expected, actual) \
try { \
  record(sec, subsec, expected, actual); \
} catch (NotImplementedException& e) { \
  missing(sec, subsec); \
}

void SUPGTest::process_point(const FEMElm& fe, int itg_point) {
  std::stringstream subsec;
  subsec << "Itg Pt " << itg_point;

  TezduyarUpwindFE stabilizer;
  stabilizer.calcSUPGWeightingFunction(fe, velocity(), nu());
  for (int i = 0; i < fe.nbf(); i++) {
    std::stringstream ss1;
    ss1 << "SUPG(" << i << ")";

    std::stringstream ss2;
    ss2 << "GP " << itg_point;

    TRY_RECORD(ss1.str(), ss2.str(), SUPG(itg_point, i), stabilizer.SUPG(i));
  }
}

void SUPGTest::report(std::ostream& stream, int* errors_out,
                       int* missing_out, bool all) const {
  int errors = 0;
  int missing_count = 0;
  for (auto it = checks_.begin(); it != checks_.end(); it++) {
    bool interesting = all;

    std::stringstream temp;
    temp << it->first << "\n";
    for (auto pt = it->second.begin(); pt != it->second.end(); pt++) {
      if (pt->second.implemented) {
        const bool passed = pt->second.error <= ERROR_THRESHOLD;
        if (all || !passed) {
          temp << "   " << (passed ? "" : FAIL_COLOR) << pt->first
               << ": " << pt->second.msg << ", err: " << pt->second.error
               << END_COLOR << "\n";
          interesting = true;
        }

        if (!passed)
          errors++;
      } else {
        temp << "   " << MISSING_COLOR << pt->first << ": " << pt->second.msg
             << END_COLOR << "\n";
        missing_count++;
        interesting = true;
      }
    }

    if (interesting)
      stream << temp.str();
  }

  if (errors_out)
    *errors_out = errors;
  if (missing_out)
    *missing_out = missing_count;
}


int SUPGTest::basis_rel_order() {
  return 0;
}

int SUPGTest::surface_id() {
  return 0;
}

ZEROPTV SUPGTest::velocity() {
  return ZEROPTV(1, 0, 0);
}

double SUPGTest::nu() {
  return 0.01;
}

double SUPGTest::calc_error(double expected, double actual) {
  // for the calc_error(0, 0) case
  if (expected == actual)
    return 0;

  // if both expected and actual are "very close" to zero
  // (defined as both within 1/2 error tolerance of 0),
  // return no error.
  // This fixes error blowing up to inf when expected is 0 and actual is 1e-17.
  static const double half_error = (ERROR_THRESHOLD / 2.0);
  if (fabs(expected) < half_error && fabs(actual) < half_error)
    return 0;

  return fabs(expected - actual) / fabs(expected);
}

double SUPGTest::calc_error(const ZEROPTV& expected, const ZEROPTV& actual) {
  return calc_error(expected.x(), actual.x()) +
         calc_error(expected.y(), actual.y()) +
         calc_error(expected.z(), actual.z());
}

