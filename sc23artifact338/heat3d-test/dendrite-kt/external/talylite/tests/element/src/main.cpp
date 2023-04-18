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
#include <element_test.h>
#include <tests/elem_1d.h>
#include <tests/elem_2d_box.h>
#include <tests/elem_2d_triangle.h>
#include <tests/elem_3d_hex.h>
#include <tests/elem_3d_tet.h>
#include <globals.h>


// tests list
ElementTest* tests[] = {
  new Elem1DTest(),
  new Elem2DBoxTest(),
  new Elem2DTriangleTest(),
  new Elem3DHexTest(),
  new Elem3DTetTest(),
};

int main(int argc, char **args) {
  PetscInitialize(&argc, &args, NULL, NULL);

  // print comparison messages always, not just for errors?
  PetscBool verbose;
  PetscOptionsHasName(NULL, NULL, "-verbose", &verbose);

  int total_errors = 0;
  for (unsigned int i = 0; i < sizeof(tests) / sizeof(tests[0]); i++) {
    // the entire test msequence is wrapped in a try block in order to
    // catch the case where an entire element is not tested. This will
    // be shown by an exception while gettting the node positions.
    try {
      ElementTest* test = tests[i];
      PrintInfo("Testing ", test->name());

      const int nodes_per_elm = get_nodes_in_element(test->elm_type(),
                                                    test->basis_order());

      GRID* pGrid = new GRID(test->basis_order());
      pGrid->set_nsd(test->nsd());
      pGrid->set_n_nodes(nodes_per_elm);
      pGrid->set_n_elements(1);
      pGrid->redimArrays(pGrid->n_nodes(), pGrid->n_elements());

      // create nodes for element
      std::vector<LocalNodeID> node_ids;
      for (int n = 0; n < nodes_per_elm; n++) {
        NODE* node = new NODE();
        ZEROPTV node_pos = test->node_positions()[n];
        node->setCoor(node_pos.x(), node_pos.y(), node_pos.z());
        pGrid->node_array_[n] = node;
        node_ids.push_back(n);
      }

      // create test element
      ELEM* elem = make_elem_of_type(test->elm_type());
      elem->redim(nodes_per_elm, node_ids.data());
      elem->set_elm_id(0);
      pGrid->elm_array_[0] = elem;

      test->test_element(elem, pGrid);

      std::stringstream ss;
      const int errors = test->report(ss);

      if (errors > 0 || verbose) {
        std::cout << ss.str();
      }

      total_errors += errors;

      delete pGrid;
    } catch (TestNotImplementedException& e) {
      // catch a completely not implemented element
      e.print();
    }
  }

  std::cout << (total_errors > 0 ? FAIL_COLOR : PASS_COLOR)
            << "\nTotal number of errors: " << total_errors
            << END_COLOR << "\n";

  PetscFinalize();
  return (total_errors != 0);
}
