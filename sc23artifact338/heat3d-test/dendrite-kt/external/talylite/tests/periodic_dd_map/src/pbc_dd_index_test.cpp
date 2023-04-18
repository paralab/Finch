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

#include <pbc_input_data.h>
#include <pbc_equation.h>

using namespace TALYFEMLIB;  // NOLINT

int main(int argc, char **args) {
  PetscInitialize(&argc, &args, NULL, NULL);

  PBCInputData inputData;
  GRID* pGrid = NULL;
  {
    GridField < NODEData > data;
    PBCEquation equation(&inputData);

    inputData.ReadFromFile();

    CreateGrid(pGrid, &inputData);

    // set up periodic boundary object
    int* pbc_indices = new int[3];
    int pbc_count = 0;
    for (int k = 0; k < inputData.nsd; k++) {
      if (inputData.boundaryType[k] == "periodic") {
        pbc_indices[pbc_count] = 2 * k + 2;
        ++pbc_count;
      }
    }
    PeriodicBounds *pbc = new PeriodicBounds(pGrid, pbc_indices, pbc_count);
    delete [] pbc_indices;

    data.redimGrid(pGrid);
    data.redimNodeData();

    int nOfDofPerNode = inputData.dof;
    PeriodicData periodic_data(pbc, nOfDofPerNode);

    // find out how many periodic variables and boundaries there are so we can
    // make that info available to preallocate
    int periodic_var_count = 0;
    for (int k = 0; k < inputData.dof; k++) {
      if (inputData.isPerVar[k] == 1) {
        periodic_var_count++;
        periodic_data.SetVarIndexPeriodic(k);
      }
    }

    equation.redimSolver(pGrid, nOfDofPerNode, false,
                         inputData.basisRelativeOrder, &periodic_data);
    equation.SetPreallocator();
    equation.PresetPeriodicData(pbc_count, periodic_var_count,
                                nOfDofPerNode);

    equation.setData(&data);
    equation.Solve();
    delete pbc;
  }

  DestroyGrid(pGrid);
  PetscFinalize();
  return 0;
}
