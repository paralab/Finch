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


using namespace TALYFEMLIB;  // NOLINT


/**
 * This class stores the data values for a single node.
 */
class CDCNodeData : public NODEData {
 public:
  double data_[2];  ///< the data for the object.

  virtual double& value(int index) {
    return data_[index];
  }

  virtual const double& value(int index) const {
    return data_[index];
  }

  static int valueno() {
    return 2;
  }
};

// for debugging, print values (only works with domain decomposition)
void print_values_dd(GRID *p_grid, GridField<CDCNodeData> *p_data) {
  for (int i_rank = 0; i_rank < GetMPISize(); i_rank++) {
    set_RankOfPrinter(i_rank);
    PrintInfo("rank ", i_rank);
    for (int i = 0; i < p_grid->n_nodes(); i++) {
      int owner;
      if (p_grid->node_belong_(i).is_owned) {
        owner = i_rank;
      } else {
        owner = p_grid->node_belong_(i).share_data(0).grid_id();
      }
      PrintInfo("i = ", i, " value: ", p_data->GetNodeData(i).value(1),
              " owner:", owner);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

int main(int argc, char **args) {
  PetscInitialize(&argc, &args, NULL, NULL);  // set up PETSc environment

  InputData input_data;  // input data object
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

      // Based on inputdata create Grid
      CreateGrid(p_grid, &input_data);

      // Construct gridfield based on Grid
      GridField<CDCNodeData> data;  // the grid field data object
      data.redimGrid(p_grid);
      data.redimNodeData();

      // store the process rank in all node data
      double rank = static_cast<double>(p_grid->grid_id());
      for (int i = 0; i < p_grid->n_nodes(); i++) {
        data.GetNodeData(i).value(1) = rank;
      }

//       print_values_dd(p_grid, &data);  // debug

      data.Communicate(1);  // sync across processes

//       print_values_dd(p_grid, &data);  // debug

      // confirm that all values are correct
      // every node not shared should have a value equal to the rank
      // every node that is shared, should have a value equal to the
      // rank of the process that owns the node.
      for (int i = 0; i < p_grid->n_nodes(); i++) {
        if (p_grid->parallel_type_ == kNoDomainDecomp ||
            p_grid->node_belong_(i).is_owned) {
          if (data.GetNodeData(i).value(1) != rank) {
            throw(std::string("incorrect value after exchange"));
          }
        } else {  // not owned byt this process
          int owner = p_grid->node_belong_(i).share_data(0).grid_id();
          if (data.GetNodeData(i).value(1) != static_cast<double>(owner)) {
            throw("incorrect value after exchange");
          }
        }
      }
    }
    DestroyGrid(p_grid);  // clean up
  } catch (const std::string& s) {  // a string was thrown due to error
    std::cerr << s << std::endl;
    PetscFinalize();
    return -1;
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    PetscFinalize();
    return -1;
  } catch (...) {  // some other, unknown exception
    std::cerr << "Unknown error" << std::endl;
    PetscFinalize();
    return -1;
  }

  PetscFinalize();  // clean up PETSc environment
  return 0;  // program exits normally
}
