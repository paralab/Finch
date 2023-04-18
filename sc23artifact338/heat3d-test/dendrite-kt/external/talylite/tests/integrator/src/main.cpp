/*
  Copyright 2017 Baskar Ganapathysubramanian

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
#include <string>

#include <cmath>  // for abs

#include <talyfem/talyfem.h>

#include <integrator_test_grid_field.h>
#include <integrator_test_input_data.h>
#include <volume_function.h>


int main(int argc, char **args) {
  int return_value = 0;
  PetscInitialize(&argc, &args, NULL, NULL);  // set up PETSc environment

  IntegratorTestInputData input_data;  // input data object
  try {
    GRID* p_grid = NULL;  // pointer to grid object

    // Read input data from file "config.txt"
    if (!input_data.ReadFromFile() || !input_data.CheckInputData()) {
      throw(std::string("[ERR] Problem with reading input data!"));
    }

    CreateGrid(p_grid, &input_data);

    // Construct gridfield based on Grid
    IntegratorTestGridField data;  // the grid field data object
    data.redimGrid(p_grid);
    data.redimNodeData();

    // Set Initial Conditions
    data.SetValues();

    // Set integration parameter
    AssemblyMethod assembly_method = kAssembleGaussPoints;
    if (input_data.use_elemental_assembly_) {
      assembly_method = kAssembleElements;
    }

    // expected sin = 2^nsd
    double expected_sin = 1.0;
    for (int i = 0; i < p_grid->nsd(); i++) {
      expected_sin *= 2.0;
    }
    double expected_cos = 0.0;
    double expected_measure = input_data.expected_measure_;

    double tolerance = 0.0001;
    PrintInfo("using tolerance: ", tolerance);

    ValueFunction<IntegratorTestNodeData> sin_integrator(
        input_data.basisFunction, input_data.basisRelativeOrder, &data, 0,
        input_data.do_optimize_, assembly_method);
    sin_integrator.set_pGrid(p_grid);
    double sin_value = sin_integrator.Solve();
    double sin_error = std::abs(sin_value - expected_sin);
    PrintInfo("sin field = ", sin_value, " (expected = ", expected_sin,
              " error = ", sin_error, ")");

    ValueFunction<IntegratorTestNodeData> cos_integrator(
        input_data.basisFunction, input_data.basisRelativeOrder, &data, 1,
        input_data.do_optimize_, assembly_method);
    cos_integrator.set_pGrid(p_grid);
    double cos_value = cos_integrator.Solve();
    double cos_error = std::abs(cos_value - expected_cos);
    PrintInfo("cos field = ", cos_value, " (expected = ", expected_cos,
              " error = ", cos_error, ")");

    VolumeFunction volume_integrator(input_data.basisFunction,
                                     input_data.basisRelativeOrder,
                                     input_data.do_optimize_, assembly_method);
    volume_integrator.set_pGrid(p_grid);
    double measure = volume_integrator.Solve();
    double measure_error = std::abs(measure - expected_measure);
    std::string measure_names[4] = {"", "length", "area", "volume"};
    PrintInfo(measure_names[p_grid->nsd()], " = ", measure, " (expected = ",
              expected_measure, " error = ", measure_error, ")");

    // check that values are correct within tolerance
    if (sin_error > tolerance) { return_value = 1; }
    if (cos_error > tolerance) { return_value = 1; }
    if (measure_error > tolerance) { return_value = 1; }

    DestroyGrid(p_grid);  // clean up
  } catch (const std::string& s) {  // a string was thrown due to error
    PetscFinalize();
    std::cerr << s << std::endl;
    return -1;
  } catch (std::bad_alloc e) {  // problem allocating memory
    PetscFinalize();
    std::cerr << "Problem with memory allocation " << e.what();
    return -1;
  // TALYException includes FileIOException
  } catch(const TALYException& e) {
    PetscFinalize();
    return -1;
  } catch (...) {  // some other, unknown exception
    std::cerr << "Unknown error" << std::endl;
    PetscFinalize();
    return -1;
  }
  PetscFinalize();  // clean up PETSc environment
  return return_value;  // program exits normally
}
