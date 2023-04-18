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
/* ------------------------------------------------------------------------

    Solid Fuel Ignition (SFI) problem.  This problem is modeled by
    the partial differential equation

            -Laplacian u - lambda*exp(u) = 0,  0 < x,y < 1,

    with boundary conditions

             u = 0  for  x = 0, x = 1, y = 0, y = 1.
             
    For checking with analytical solution, we solve 1D problem given by
    
            u" - pi^2 * exp(u) = 0,  0 < x < 1,
            u(0) = u(1) = 0
            
            u(x) = -ln(1+cos(pi*(0.5+x)))
    It is solved for -0.4 <= x <= 0.4 range, as u(x=0.5)=inf

  ------------------------------------------------------------------------- */

#include <talyfem/talyfem.h>
#include <math.h>
#include <fstream>
#include <string>

#include <BTInputData.h>
#include <BTNodeData.h>
#include <BTGridField.h>
#include <BTEquation.h>


using namespace TALYFEMLIB;
static char help[] = "Solves a transient transfer problem!";


inline bool SetIC(BTGridField& data, BTInputData& idata) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  switch (idata.typeOfIC) {
    case 1:
      data.SetIC(idata.nsd);
      return true;
      break;
    case 0:
      try {
        load_gf(&data, &idata);
      } catch(const TALYException& e) {
        e.print();
        PrintWarning("Failed to load GridField data!");
        return false;
      }

      data.UpdateDataStructures();
      return true;
      break;
    default:
      if (rank == 0) std::cerr << "IC not set up " << std::endl;
      return false;
  }
  return false;
}

int main(int argc, char **args) {
  PetscInitialize(&argc, &args, (char *)0, help);
  try {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    BTInputData  inputData;
    GRID*  pGrid = NULL;
    {
    BTGridField  data;

    if ((argc == 2) && (strcmp(args[1], "--printParameters") == 0)) {
      if (rank == 0) { inputData.printAll(std::cout); }
      PetscFinalize();
      return -1;
    }

    // Read input data from file "config.txt"
    inputData.ReadFromFile();

    // Loging inputdata to std::clog or file
    if (rank == 0) std::clog << inputData;
    // std::ofstream f_out_logs("InputData.log");
    // if(rank == 0) f_out_logs << inputData;
    // f_out_logs.close();

    //  Check if inputdata is complete
    if (!inputData.CheckInputData()) {
      throw(std::string(
          "[ERR] Problem with input data, check the config file!"));
    }

    // Based on inputdata create Grid
    CreateGrid(pGrid, &inputData);

    // Translating the grid coordinates from [0,Lx] to [-0.4,0.4]
    double xl = -0.4, xr = 0.4;
    double Lx = inputData.L[0];
    for (int A = 0; A < pGrid->n_nodes(); A++) {
      double x = pGrid->GetCoord(A, 0);
      double x_new = xl + x*(xr-xl)/Lx;
      pGrid->node_array_[A]->setCoor(0, x_new);
    }

    // Construct gridfield based on Grid
    data.redimGrid(pGrid);
    data.redimNodeData();

     // Set Initial Conditions
    if (!SetIC(data, inputData)) {
      PrintResults("Failed to set initial conditions", false,
                   inputData.shouldFail);
      delete pGrid;
      throw(std::string("Problem with IC, not loaded"));
    }

    // Set Solver parameters
    int nOfDofPerNode = 1;  // number of degree of freedom per node
    bool do_accelerate = inputData.use_elemental_assembly_;
    AssemblyMethod assembly_method = kAssembleGaussPoints;
    if (inputData.use_elemental_assembly_) {
      assembly_method = kAssembleElements;
    }

    BTEquation heatEq(&inputData, do_accelerate,
                      assembly_method);  // the equation solver

    heatEq.redimSolver(pGrid, nOfDofPerNode, false, inputData.basisRelativeOrder);
    heatEq.setData(&data);
    PrintStatus("solver initialized!", rank);

    heatEq.Solve();
    PrintStatus("Finished Solve!", rank);

    data.SetAnalyticalSolution();
    PrintStatus("analytical solution set!", rank);
    // if (inputData.ifPrintPltFiles) {
    //   data.writeNodeDataToFile(&inputData, "data_final.plt", t);
    // }

    // Check if test passed and print result
    PrintStatus("Checking solution...");
    // calculate error in final solution
    const double error = data.CalcMaxError(inputData.ifDD);
    std::stringstream string_stream;  // error string to print
    string_stream << "error = " << error;
    const double max_error = 2e-2;  // error tolerance
    // there should always be some error in the calculation. A value of zero
    // for the error nearly always means there's a bug in the code.
    const double min_error = 0.00001;
    const bool test_pass = (error < max_error) && (error > min_error);
    // finally, print the result (pass/fail) of the calculation
    PrintResults(string_stream.str(), test_pass, inputData.shouldFail);

    if (inputData.ifPrintPltFiles) {
      save_gf(&data, &inputData, "data.plt", 0.0);
    }
    }
    // clean
    DestroyGrid(pGrid);
  }
  catch (const std::string& s) {
    PetscFinalize();
    std::cerr << s << std::endl;
    return -1;
  }
  catch (std::bad_alloc e) {
    PetscFinalize();
    std::cerr << "Problem with memory allocation " << e.what();
    return -1;
  }
  catch (const TALYException& e) {
    e.print();
    PetscFinalize();
    return -1;
  }
  catch (...) {
    std::cerr << "Unknown error" << std::endl;
    PetscFinalize();
    return -1;
  }

  PetscFinalize();
  return 0;
}
