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

#include <string.h>

#include <vector>

#include <talyfem/talyfem.h>
#include <stabilizertest.h>
#include <globals.h>

// tests
#include <tests/box1d_linear.h>
#include <tests/box1d_quadratic.h>
#include <tests/box1d_cubic.h>

#include <tests/box2d_linear.h>
#include <tests/box2d_quadratic.h>
#include <tests/box2d_cubic.h>

#include <tests/box3d_linear.h>
#include <tests/box3d_quadratic.h>
#include <tests/box3d_cubic.h>

#include <tests/tri2d_linear.h>
#include <tests/tet3d_linear.h>

// tests list
SUPGTest* tests[] = {
  new Box1DLinearStabilizerTest(),
  new Box1DQuadraticStabilizerTest(),
  new Box1DCubicStabilizerTest(),

  new Box2DLinearStabilizerTest(),
  new Box2DQuadraticStabilizerTest(),
  new Box2DCubicStabilizerTest(),

  new Box3DLinearStabilizerTest(),
  new Box3DQuadraticStabilizerTest(),
  new Box3DCubicStabilizerTest(),

  new Tri2DLinearStabilizerTest(),

  new Tet3DLinearStabilizerTest(),
};

int main(int argc, char **args) {
  PetscInitialize(&argc, &args, NULL, NULL);

  // print comparison messages always, not just for errors?
  PetscBool verbose;
  PetscOptionsHasName(NULL, NULL, "-verbose", &verbose);

  int total_errors = 0;
  int total_missing = 0;
  for (unsigned int i = 0; i < sizeof(tests) / sizeof(tests[0]); i++) {
    SUPGTest* test = tests[i];
    PrintInfo("Testing ", test->name());

    const int mesh_order = basis_get_mesh_order(test->basis_function());
    const int nodes_per_elm = get_nodes_in_element(test->elm_type(),
                                                   mesh_order);

    GRID* pGrid = new GRID(mesh_order);
    pGrid->set_nsd(test->nsd());
    pGrid->set_grid_type(test->grid_type());
    pGrid->set_n_nodes(nodes_per_elm);
    pGrid->set_n_elements(1);
    pGrid->redimArrays(pGrid->n_nodes(), pGrid->n_elements());

    // create nodes for element
    std::vector<LocalNodeID> node_ids;
    for (int n = 0; n < nodes_per_elm; n++) {
      NODE* node = new NODE();
      ZEROPTV node_pos = test->node_position(n);
      node->setCoor(node_pos.x(), node_pos.y(), node_pos.z());
      pGrid->node_array_[n] = node;
      node_ids.push_back(n);
    }

    // create test element
    ELEM* elem = make_elem_of_type(test->elm_type());
    elem->redim(nodes_per_elm, node_ids.data());
    elem->set_elm_id(0);
    pGrid->elm_array_[0] = elem;

    FEMElm femelm(pGrid, BASIS_ALL);

    try {
      if (test->surface_id() == 0) {
        // 0 = element ID
        femelm.refill(0, test->basis_rel_order());
      } else {
        SurfaceIndicator surface(test->surface_id());

        // 0 = element ID
        femelm.refill_surface(0, &surface,
                              test->basis_rel_order());
      }
      assert(femelm.basis_function() == test->basis_function());
    } catch (TALYException& e) {
      std::cout << MISSING_COLOR << "FEMElm::refill: " << e.what()
                << END_COLOR << std::endl;
      total_missing += 1;
      delete pGrid;
      continue;
    }

    while (femelm.next_itg_pt()) {
      test->process_point(femelm, femelm.cur_itg_pt_num());
    }

    std::stringstream ss;
    int errors;
    int missing;
    test->report(ss, &errors, &missing, verbose);

    if (errors > 0 || verbose) {
      std::cout << ss.str();
    }
    if (missing > 0) {
      std::cout << MISSING_COLOR << "  (" << missing << " missing tests)"
                << END_COLOR << std::endl;
    }

    total_errors += errors;
    total_missing += missing;

    delete pGrid;
  }

  std::cout << (total_errors > 0 ? FAIL_COLOR : PASS_COLOR)
            << "\nTotal number of errors: " << total_errors
            << END_COLOR << "\n";

  if (total_missing > 0) {
    std::cout << MISSING_COLOR << total_missing << " tests have not been "
              << "implemented, and were skipped. Run with "
              << "./basis_test -verbose for more details.\n" << END_COLOR;
  }

  PetscFinalize();
  return (total_errors != 0);
}

