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
#include <talyfem/talyfem.h>

#include <algorithm>  // for std::sort
#include <vector>

#include <pbc_test_input_data.h>


using namespace TALYFEMLIB;  // NOLINT

// Implement Solve so we can create an object.
template<class NodeData>
class CEquationTest : public CEquation<NodeData> {
 public:
  int get_n_periodic_vars() const {
    return this->NumPeriodicVars();
  }
  int get_ith_periodic_var(int index) {
    return this->IthPeriodicVar(index);
  }
  virtual void Solve(double delta_t, double current_time) { }
};

int main(int argc, char **args) {
  PetscInitialize(&argc, &args, NULL, NULL);

  PBCTestInputData inputData;
  GRID* pGrid = NULL;
  {
    GridField < NODEData > data;
    CEquationTest<NODEData> equation;
    inputData.ReadFromFile();

    CreateGrid(pGrid, &inputData);

    // set up periodic boundary object
    PeriodicBounds *pbc = new PeriodicBounds(pGrid, inputData.perBounds,
                                             inputData.nPerBounds);

    data.redimGrid(pGrid);
    data.redimNodeData();

    // find out how many periodic variables and boundaries there are so we can
    // make that info available to preallocate
    int periodic_var_count = inputData.nPerVars;
    int periodic_bound_count = inputData.nPerBounds;

    PeriodicData periodic_data(pbc, inputData.nDoF);

    // convert the list of periodic variables to a vector in order to sort
    // the results
    std::vector<int> sorted_var_list(inputData.perVars,
                                     inputData.perVars + inputData.nPerVars);
    std::sort(sorted_var_list.begin(), sorted_var_list.end());

    for (int i = 0; i < inputData.nPerVars; i++) {
      periodic_data.SetVarIndexPeriodic(sorted_var_list[i]);
    }

    equation.redimSolver(pGrid, inputData.nDoF, false,
                         inputData.basisRelativeOrder, &periodic_data);
    equation.SetPreallocator();
    equation.PresetPeriodicData(periodic_bound_count, periodic_var_count,
                                inputData.nDoF);
    equation.setData(&data);

    // fill values based on input file, then print values
    // list of periodic variables should be sorted when printed
    std::cout << "number of periodic boundaries: "
        << pbc->n_periodic_bounds() << std::endl;
    for (int i = 0; i < 2 * inputData.nsd; i++) {
      std::cout << "is boundary " << i + 1 << " periodic? "
          << pbc->IsBoundaryPeriodic(i) << std::endl;
    }

    std::cout << "number of periodic variables: "
        << equation.get_n_periodic_vars() << std::endl;
    for (int i = 0; i < inputData.nDoF; i++) {
      std::cout << "is variable " << i + 1 << " periodic? "
          << periodic_data.is_var_periodic()[i] << std::endl;
    }

    for (int i = 0; i < equation.get_n_periodic_vars(); i++) {
      std::cout << "periodic variable " << i + 1 << " = "
          << equation.get_ith_periodic_var(i) << std::endl;
    }
    delete pbc;
  }
  DestroyGrid(pGrid);
  PetscFinalize();
  return 0;
}
