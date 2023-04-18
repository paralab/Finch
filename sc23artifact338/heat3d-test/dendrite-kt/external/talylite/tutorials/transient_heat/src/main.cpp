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
#include <math.h>

#include <fstream>
#include <string>
#include <sstream>  // std::stringstream

#include <talyfem/talyfem.h>

#include "HTAnalyticSolution.h"
#include "HTInputData.h"
#include "HTGridField.h"
#include "HTEquation.h"

using namespace TALYFEMLIB;

static char help[] = "Solves a transient transfer problem!";


/**
 * Sets the initial conditions for the system
 *
 * There are two options for the initial conditions and the choice is controlled
 * by the typeOfIC input value. When set to 1, this will call the initial
 * condition function in HTGridField. When set to 0, it will read in the values
 * from a file.
 *
 * @param data the grid field data object
 * @param input_data the input data object
 * @return true if the initial conditions were set correctly
 */
inline bool SetInitialConditions(HTGridField& data, HTInputData& input_data) {
  switch (input_data.typeOfIC) {
    case 1:  // calculate intial conditions from the analytic solution
      data.SetInitialConditions();
      return true;
    case 0:  // read the intial conditions from a file
      try {
        load_gf(&data, &input_data);
      } catch(TALYException& e) {
        e.print();
        PrintError("Failed to load GridField data.");
        return false;
      }

      // to prepare for the first time step, copy u to u_pre in the node data
      data.UpdateDataStructures();
      return true;
    default:  // invalid value given for typeOfIC
      if (GetMPIRank() == 0) {
        std::cerr << "IC not set up " << std::endl;
      }
      return false;
  }
}


int main(int argc, char **args) {
  PetscInitialize(&argc, &args, NULL, help);  // set up PETSc environment

  Repro r(argc, args);
  r.write();

  PrintInfo("sizeof(PetscInt): ", sizeof(PetscInt));
  /*int a = 1300;
  a *= 1300;
  a *= 1300;
  PrintInfo("Attempted overflow");*/

  HTInputData input_data;  // input data object
  try {
    GRID* p_grid = NULL;  // pointer to grid object

    // the following { } are needed to limit the scope of the solver variable
    // and ensure is goes out of scope prior to calling PetscFinalize.
    // Otherwise, one of the destructors associated with the solver may call
    // something PETSc related after the finalize routine has been called.
    // This would lead to a crash.
    {
      // Read input data from file "config.txt"
      if (!input_data.ReadFromFile()) {
        throw(std::string("[ERR] Problem with reading input data!"));
      }

      // Check if input data is complete
      if (!input_data.CheckInputData()) {
        throw(std::string(
            "[ERR] Problem with input data, check the config file!"));
      }

      if (input_data.use_isobox_) {
        PrintInfo("Creating Isobox mesh to avoid using ParMETIS...");

        p_grid = new GRID();
        CMeshPartition part;

        PetscInt scale = (PetscInt) std::round(std::pow(GetMPISize(), 1.0 / 3.0));
        assert(scale * scale * scale == GetMPISize());

        PetscInt lclN = input_data.Nelem[0] / scale;
        assert(input_data.Nelem[0] % scale == 0);

        PrintInfo("  lclN: ", lclN);
        PrintInfo("  L: ", input_data.L[0]);
        part.CreateIsoBox3D(lclN, input_data.L[0]);

        part.TransferToGrid(p_grid);
        part.PartitionFree();
      } else {
        // Based on inputdata create Grid
        CreateGrid(p_grid, &input_data);
      }

      // Create analytic solution
      HTAnalyticSolution sol(input_data.nsd, p_grid);

      // Construct gridfield based on Grid
      HTGridField data(&sol);  // the grid field data object
      data.redimGrid(p_grid);
      data.redimNodeData();

      // if we're not generating the grid automatically,
      // we need to manually set boundary indicators
      if (!(input_data.ifBoxGrid) && !(input_data.ifLoadNodeIndicators)) {
        if (input_data.ifDD)
          PrintWarning("Attempting to guess node boundaries with DD - this will fail if we're on more than one process!");
        else
          PrintInfo("Generating node boundaries for loaded mesh (will fail if not cube-shaped)...");
        data.SetBoundaries();
      }

      // Set Initial Conditions
      if (!SetInitialConditions(data, input_data)) {
        delete p_grid;
        PrintResults("Loading GridField data failed.", false,
                     input_data.should_fail_);
        throw(std::string("Problem with IC, not loaded"));
      }

      // Set Solver parameters
      int n_dof = 1;  // number of degree of freedom per node
      bool do_accelerate = input_data.use_elemental_assembly_;
      AssemblyMethod assembly_method = kAssembleGaussPoints;
      if (input_data.use_elemental_assembly_) {
        assembly_method = kAssembleElements;
      }

      HTEquation heat_equation(&input_data, &sol, do_accelerate,
                               assembly_method);  // the equation solver
      heat_equation.redimSolver(p_grid, n_dof, false,
                                input_data.basisRelativeOrder);  // intialize solver
      heat_equation.setData(&data);
      PrintStatus("Solver initialized.");

      double dt = input_data.dt_;  // time step for time stepping solver
      double t = 0;  // intial time
      // if specified in the input file, start from a time other than 0
      if (input_data.restart_from_t_ > 0) {
        t = input_data.restart_from_t_;
      }

      // write out the initial data to file
      if (!input_data.output_extension_.empty() && input_data.output_ic_
          && input_data.output_extension_ != "none") {
        std::string filename = "data_init" + input_data.output_extension_;
        save_gf(&data, &input_data, filename.c_str(), t);
        PrintStatus("Initial data sent to file.");
      }

      int n_time_steps = input_data.n_time_steps_;  // number of time steps
      // perform solving loop
      for (int i = 0; i < n_time_steps; i++) {
        // turn off matrix recalculation if this is the 2nd step
        // we need it on for the 1st step to caclulate the initial matrix
        if (i == 1 && input_data.use_vector_only_assembly_) {
          heat_equation.set_recalc_matrix(false);
        }
        t += dt;  // update time to new value
        PrintInfo("Solving timestep ", i, "...");
        heat_equation.Solve(dt, t);  // solve the equation
        PrintInfo("   Updating data structures...");
        data.UpdateDataStructures();  // update data for the next cycle
      }
      PrintStatus("Time stepping completed.");

      // calculate the analytic solution at the end to check solver accuracy
      data.SetAnalyticalSolution(t, input_data.nsd);
      PrintStatus("Analytical solution set.");

      // write out the final data to file
      if (!input_data.output_extension_.empty() &&
          input_data.output_extension_ != "none") {
        std::string filename = "data_final" + input_data.output_extension_;
        save_gf(&data, &input_data, filename.c_str(), t);
        PrintStatus("Final data sent to file.");
      }

      // Check if test passed and print result
      PrintStatus("Checking solution...");
      // calculate error in final solution
      const double error = data.CalcMaxError(input_data.ifDD);
      std::stringstream string_stream;  // error string to print
      string_stream << "error = " << error;
      const double max_error = 0.01;  // error tolerance
      // there should always be some error in the calculation. A value of zero
      // for the error nearly always means there's a bug in the code.
      const double min_error = 0.00001;
      const bool test_pass = (error < max_error) && (error > min_error);
      // finally, print the result (pass/fail) of the calculation
      PrintResults(string_stream.str(), test_pass, input_data.should_fail_);
    }
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
    PrintResults(e.what(), false, input_data.should_fail_);
    PetscFinalize();
    return -1;
  }

  PetscFinalize();  // clean up PETSc environment
  return 0;  // program exits normally
}
