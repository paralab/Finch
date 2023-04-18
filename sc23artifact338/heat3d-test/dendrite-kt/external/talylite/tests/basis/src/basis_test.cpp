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

#include <basistest.h>
#include <globals.h>


#define TRY_RECORD(sec, subsec, expected, actual) \
try { \
  record(sec, subsec, expected, actual); \
} catch (NotImplementedException& e) { \
  missing(sec, subsec); \
}

int BasisTest::num_nodes() {
  int mesh_order = basis_get_mesh_order(basis_function());
  return get_nodes_in_element(elm_type(), mesh_order);
}

void BasisTest::process_point(const FEMElm& fe, int itg_point) {
  std::stringstream subsec;
  subsec << "Itg Pt " << itg_point;
  TRY_RECORD("Gauss Points", subsec.str(), gauss_point(itg_point),
         fe.itg_pt());

  TRY_RECORD("|Jacobian|", subsec.str(), jacobian_det(itg_point),
             fabs(fe.jacc()))
  TRY_RECORD("|Jacobian| * weight", subsec.str(),
             jacobian_det_times_weight(itg_point), fe.detJxW())

  // take surface ID into account for nsd - this is important, so we use
  // the right array/matrix dimensions when checking surface integration
  int iso_nsd = nsd() - (surface_id() != 0 ? 1 : 0);
  int iso_n_nodes = fe.nbf();

  for (int i = 0; i < iso_n_nodes; i++) {
    subsec.str(""); subsec.clear();
    subsec << "GP " << i;

    std::stringstream ss;
    ss << "N at Integration Point " << itg_point;
    TRY_RECORD(ss.str(), subsec.str(), N(itg_point, i), fe.N(i))

    ss.str(""); ss.clear();
    ss << "dNde at Integration Point " << itg_point;
    try {
      ZEROPTV test_dNde(dNde(itg_point, i, 0),
                        iso_nsd > 1 ? dNde(itg_point, i, 1) : 0,
                        iso_nsd > 2 ? dNde(itg_point, i, 2) : 0);

      ZEROPTV actual_dNde(fe.dNde(i, 0),
                          iso_nsd > 1 ? fe.dNde(i, 1) : 0,
                          iso_nsd > 2 ? fe.dNde(i, 2) : 0);

      TRY_RECORD(ss.str(), subsec.str(), test_dNde, actual_dNde)
    } catch (NotImplementedException& e) {
      missing(ss.str(), subsec.str());
    }
  }

  std::stringstream ss;
  ss.str(""); ss.clear();
  ss << "dXde at Integration Point " << itg_point;
  for (int i = 0; i < iso_nsd; i++) {
    for (int j = 0; j < iso_nsd; j++) {
      subsec.str(""); subsec.clear();
      subsec << "dXde(" << i << ", " << j << ")";
      TRY_RECORD(ss.str(), subsec.str(), dXde(itg_point, i, j), fe.dXde(i, j))
    }
  }

  ss.str(""); ss.clear();
  ss << "cof at Integration Point " << itg_point;
  for (int i = 0; i < iso_nsd; i++) {
    for (int j = 0; j < iso_nsd; j++) {
      subsec.str(""); subsec.clear();
      subsec << "cof(" << i << ", " << j << ")";
      TRY_RECORD(ss.str(), subsec.str(), cof(itg_point, i, j), fe.cof(i, j))
    }
  }

  ss.str(""); ss.clear();
  ss << "dN at Integration Point " << itg_point;
  for (int i = 0; i < iso_n_nodes; i++) {
    for (int axis = 0; axis < iso_nsd; axis++) {
      subsec.str(""); subsec.clear();
      subsec << "dN(" << i << ", " << axis << ")";
      TRY_RECORD(ss.str(), subsec.str(), dN(itg_point, i, axis), fe.dN(i, axis))
    }
  }

  ss.str(""); ss.clear();
  ss << "d2Nde at Integration Point " << itg_point;
  int d2_dim = iso_nsd * (iso_nsd + 1) / 2;
  int nbf = iso_n_nodes;
  for (int i = 0; i < nbf; i++) {
    for (int axis = 0; axis < d2_dim; axis++) {
      subsec.str(""); subsec.clear();
      subsec << "d2Nde(" << i << ", " << axis << ")";
      TRY_RECORD(ss.str(), subsec.str(), d2Nde(itg_point, i, axis), fe.d2Nde(i, axis))
    }
  }

  ss.str(""); ss.clear();
  ss << "d2N at Integration Point " << itg_point;
  for (int i = 0; i < nbf; i++) {
    for (int axis = 0; axis < d2_dim; axis++) {
      subsec.str(""); subsec.clear();
      subsec << "d2N(" << i << ", " << axis << ")";
      TRY_RECORD(ss.str(), subsec.str(), d2N(itg_point, i, axis), fe.d2N(i, axis))
    }
  }
}

void BasisTest::report(std::ostream& stream, int* errors_out,
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

int BasisTest::surface_id() {
  return 0;
}

double BasisTest::calc_error(double expected, double actual) {
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

double BasisTest::calc_error(const ZEROPTV& expected, const ZEROPTV& actual) {
  return calc_error(expected.x(), actual.x()) +
         calc_error(expected.y(), actual.y()) +
         calc_error(expected.z(), actual.z());
}

