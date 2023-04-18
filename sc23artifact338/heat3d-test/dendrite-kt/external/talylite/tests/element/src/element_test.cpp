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

#include <element_test.h>

#include <vector>
#include <utility>  // std::pair

#include <globals.h>


int ElementTest::num_nodes() {
  return get_nodes_in_element(elm_type(), basis_order());
}

void ElementTest::test_element(ELEM *element, GRID *p_grid) {
  std::stringstream subsec;
  subsec << "TEST ";

  try {
    // force measure() to evaluate before element->GetMeasure() so we can check
    // if it's implemented (as element->GetMeasure() may throw a TALYException)
    // (C++ does not guarantee argument evaluation order)
    double val = measure();
    record("Measure", subsec.str(), val, element->GetMeasure(p_grid));
  } catch(TestNotImplementedException& e) {
    e.print();
  }

  try {
    record("Nodes Per Surface", subsec.str(), nodes_per_surface(),
          element->GetNodesPerSurface());
  } catch(TestNotImplementedException& e) {
    e.print();
  }

  try {
    record("Surface Count", subsec.str(), surface_count(),
          element->GetSurfaceCount());
  } catch(TestNotImplementedException& e) {
    e.print();
  }

  for (int i = 0; i < element->GetSurfaceCount(); i++) {
    try {
      const int* surf_check = element->GetSurfaceCheckArray();
      int surf_id = surf_check[i * (element->GetNodesPerSurface() + 1)];

      std::stringstream surfss;
      surfss << "Normal for surface " << surf_id;
      record(surfss.str(), subsec.str(), normal(surf_id),
             element->CalculateNormal(p_grid, surf_id));
    } catch(TestNotImplementedException& e) {
      e.print();
    }
  }

  try {
    std::vector< std::pair<ZEROPTV, bool> > test_points = inner_points();
    for (unsigned int i = 0; i < test_points.size(); i++) {
      std::stringstream ptss;
      ptss << "Inner point test #" << (i + 1);
      record(ptss.str(), subsec.str(), test_points[i].second,
             element->IsInnerPoint(p_grid, test_points[i].first));
    }
  } catch (TestNotImplementedException& e) {
    e.print();
  }

  try {
    std::vector<ZEROPTV> localPts;
    if (element->elmType() == kElem2dTriangle || element->elmType() == kElem3dTetrahedral) {
      double start = 0.0;
      double end = 1.0;
      double step = 0.1;
      int elem_nsd = element->nsd();

      // loop over the space from (0, 0, 0) to (1, 1, 1) in uniform steps
      // since some points in this space are obviously outside the unit
      // triangle/tetrahedron, we discard any points in that space that do
      // not lie inside the triangle/tet.

      // Any point inside a unit tetrahedron should satisfy:
      // 0 <= x <= 1,
      // 0 <= y <= 1,
      // 0 <= z <= 1,
      // 0 <= (1 - x - y - z) <= 1.
      // (This is a special case of barycentric coordinates - we don't need
      // to normalize x/y/z since we are dealing with a unit triangle/tet.)
      // The first 3 cases must be true by our loop bounds, so all we
      // need to check is the last case.

      for (double k = (elem_nsd >= 3 ? start : 0.0); k <= (elem_nsd >= 3 ? end : 0.0); k += step) {
        for (double j = (elem_nsd >= 2 ? start : 0.0); j <= (elem_nsd >= 2 ? end : 0.0); j += step) {
          for (double i = (elem_nsd >= 1 ? start : 0.0); i <= (elem_nsd >= 1 ? end : 0.0); i += step) {
            if (1.0 - i - j - k < 0)
              continue;
            localPts.push_back(ZEROPTV(i, j, k));
          }
        }
      }

    } else if (element->elmType() == kElem1d
               || element->elmType() == kElem2dBox
               || element->elmType() == kElem3dHexahedral) {

      double start = -1.0;
      double end = 1.0;
      double step = 0.25;
      int elem_nsd = element->nsd();
      for (double k = (elem_nsd >= 3 ? start : 0.0); k <= (elem_nsd >= 3 ? end : 0.0); k += step) {
        for (double j = (elem_nsd >= 2 ? start : 0.0); j <= (elem_nsd >= 2 ? end : 0.0); j += step) {
          for (double i = (elem_nsd >= 1 ? start : 0.0); i <= (elem_nsd >= 1 ? end : 0.0); i += step) {
            localPts.push_back(ZEROPTV(i, j, k));
          }
        }
      }
    } else {
      throw TestNotImplementedException("GetLocalPtv", name());
    }

    for (const ZEROPTV& localPt : localPts) {
      ZEROPTV ptvg;
      {
        FEMElm fe(p_grid, BASIS_ALL);
        fe.refill(element, BASIS_LINEAR, 0);
        fe.calc_at(localPt);
        ptvg = fe.position();
      }

      std::stringstream sec;
      sec << "GetLocalPtv for " << ptvg;

      try {
        ZEROPTV found_ptvl;
        p_grid->GetLocalPtv(ptvg, found_ptvl, element->elm_id());

        // use L2 norm as error
        double err = 0.0;
        for (int i = 0; i < nsd(); i++) {
          double diff = (found_ptvl(i) - localPt(i));
          err += diff*diff;
        }
        err = sqrt(err);

        record_abs_error(sec.str(), subsec.str(), 0.0, err);
      } catch (TALYException& e) {
        record_abs_error(sec.str(), subsec.str(), 0.0, std::numeric_limits<double>::quiet_NaN());
      }
    }
  } catch (TestNotImplementedException& e) {
    e.print();
  }

  // need to test contains
  // need to test number of nodes
  // need to test node id array (seems dumb, but this could catch off by ones.
}

int ElementTest::report(std::ostream& stream) const {
  int errors = 0;
  for (auto it = checks_.begin(); it != checks_.end(); it++) {
    stream << it->first << "\n";
    for (auto pt = it->second.begin(); pt != it->second.end(); pt++) {
      const bool passed = pt->second.error <= ERROR_THRESHOLD;
      stream << "   " << (passed ? "" : FAIL_COLOR) << pt->first
             << ": " << pt->second.msg << ", err: " << pt->second.error
             << END_COLOR << "\n";

      if (!passed)
        errors++;
    }
  }
  return errors;
}
