/*
  Copyright 2014-2017 Baskar Ganapathysubramanian

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
#ifndef TESTS_PERIODICDDMAP_INCLUDE_PBC_EQUATION_H_
#define TESTS_PERIODICDDMAP_INCLUDE_PBC_EQUATION_H_

#include <mpi.h>

#include <sstream>

// since assert might be disabled by the -DNDEBUG flag, this is used instead.
inline void _always_assert(const char* expression, const char* file, int line) {
  fprintf(stderr, "Assertion '%s' failed, file '%s' line '%d'.", expression,
          file, line);
  abort();
}
#define always_assert(EXPRESSION) ((EXPRESSION) ? (void)0 : \
    _always_assert(#EXPRESSION, __FILE__, __LINE__))


class PBCEquation : public CEquation<NODEData> {
 public:
  PBCInputData* input_data_;

  explicit PBCEquation(PBCInputData* input_data)
      : input_data_(input_data) {
  }

  virtual void fillEssBC() {
    for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
      for (int k = 0; k < input_data_->nsd; k++) {
        if (input_data_->boundaryType[k] == "direchlet") {
          // index of boundary nodes (BoNode values) starts at 1, while the
          // index of boundary values starts at zero.
          if (p_grid_->BoNode(nodeID, 2 * k + 1)) {
            specifyValue(nodeID, 0, 0.0);
          }
          if (p_grid_->BoNode(nodeID, 2 * k + 2)) {
            specifyValue(nodeID, 0, 0.0);
          }
        }
      }
    }
  }

  virtual void Solve(double delta_t = 0.0, double current_time = 0.0) {
    int* pbc_indices = new int[3];
    int pbc_count = 0;

    for (int k = 0; k < input_data_->nsd; k++) {
      if (input_data_->boundaryType[k] == "periodic") {
        pbc_indices[pbc_count] = 2 * k + 2;
        ++pbc_count;
      }
    }

    initEssBC();
    fillEssBC();

    if (pbc_count > 0) {
      int per_var_count = 0;
      for (int k = 0; k < input_data_->dof; k++) {
        if (input_data_->isPerVar[k] == 1) {
          per_var_count++;
        }
      }
      if (per_var_count == 0) {
        PrintError("must have at least one periodic variable");
        exit(1);
      }
    }

    ApplyEssBCToSolution();
    Assemble();
    ApplyAllBC();

    CheckMapping();

    delete[] pbc_indices;
  }

  // check the mapping of where nodes are mapped in the index_from mapping
  void CheckMapping() {
    const int total_dof = p_grid_->n_total_nodes() * n_dof_;

    const int kUnknown = -1;
    always_assert(kUnknown < 0);  // this MUST be negative or it will not work.

    // create array to store from values for each unknown on this process. This
    // mapping is done based on the physical id of the data. (-1) means this
    // unknown is not on this process.
    ZEROARRAY<int> from_mapping;
    from_mapping.redim(total_dof);
    from_mapping.fill(kUnknown);
    ZEROARRAY<int> sol_mapping;
    sol_mapping.redim(total_dof);
    sol_mapping.fill(kUnknown);

    // fill the from mapping. the index is based on the physical mapping.
    for (int i_node = 0; i_node < p_grid_->n_nodes(); i_node++) {
      for (int i_dof = 0; i_dof < n_dof_; i_dof++) {
        const int lcl_idx = i_node * n_dof_ + i_dof;
        const int phys_idx = p_grid_->physical_map(i_node) * n_dof_ + i_dof;
        const int sol_idx = p_grid_->solution_map(i_node) * n_dof_ + i_dof;
        from_mapping(phys_idx) = idx_from_(lcl_idx);
        sol_mapping(phys_idx) = sol_idx;
      }
    }

    const int mpi_size = GetMPISize();
    const int mpi_rank = GetMPIRank();

    // this stores all the data needed by rank 0
    ZEROMATRIX<int> from_all(mpi_size, total_dof);
    ZEROMATRIX<int> sol_all(mpi_size, total_dof);

    // copy all data to rank 0
    MPI_Gather(from_mapping.data(), total_dof, MPI_INT, from_all.data_ptr(),
               total_dof, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(sol_mapping.data(), total_dof, MPI_INT, sol_all.data_ptr(),
               total_dof, MPI_INT, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
      // make sure all received data is sane
      for (int i = 0; i < mpi_size; i++) {
        for (int j = 0; j < total_dof; j++) {
          // value must be between 0 and total_dof - 1
          // OR value must be kUnknown
          always_assert(from_all(i, j) == kUnknown ||
                        (from_all(i, j) >= 0 && from_all(i, j) < total_dof));
          always_assert(sol_all(i, j) == kUnknown ||
                        (sol_all(i, j) >= 0 && sol_all(i, j) < total_dof));
        }
      }

      // combine all arrays to get the total mapping for all processes.
      // this amounts to taking the max value for each unknown and putting it
      // into the combined array
      ZEROARRAY<int> combined_from;
      combined_from.redim(total_dof);
      combined_from.fill(kUnknown);
      ZEROARRAY<int> combined_sol;
      combined_sol.redim(total_dof);
      combined_sol.fill(kUnknown);
      for (int j = 0; j < total_dof; j++) {
        int from_value = kUnknown;
        int sol_value = kUnknown;
        for (int i = 0; i < mpi_size; i++) {
          if (from_all(i, j) > from_value) { from_value = from_all(i, j); }
          if (sol_all(i, j) > sol_value) { sol_value = sol_all(i, j); }
        }
        // this should always be greater >= 0 since we have the global value
        always_assert(from_value >= 0);
        always_assert(sol_value >= 0);
        combined_from(j) = from_value;
        combined_sol(j) = sol_value;
      }

      // There should only be one value for each unknown. Make sure every
      // process that has a value assigned to the unknown has the correct value.
      // The value a process has must be either kUnknown or the value from the
      // combined array
      for (int j = 0; j < total_dof; j++) {
        const int expected_from_value = combined_from(j);
        const int expected_sol_value = combined_sol(j);
        for (int i = 0; i < mpi_size; i++) {
          always_assert(from_all(i, j) == kUnknown ||
                        from_all(i, j) == expected_from_value);
          always_assert(sol_all(i, j) == kUnknown ||
                        sol_all(i, j) == expected_sol_value);
        }
//         PrintError( j, " " ,combined_from(j), " ", combined_sol(j));
      }

      // the information below is hard coded for a 16 node (9 element) 2D mesh.
      // if the code isn't being run on that size mesh, these tests are invalid.
      always_assert(p_grid_->n_total_nodes() == 16);
      // this only supports up to 4 degrees of freedom
      always_assert(n_dof_ <= 4);

      // now that we know the values are reasonable, we can check the values
      // for the proper mapping. First the expected mapping of nodes is set up
      // based on the periodic bounds.
      ZEROARRAY<int> expected_from_map;

      expected_from_map.redim(p_grid_->n_total_nodes() * n_dof_);
      // we're not testing the solution mapping, so assume this is correct
      // and copy this into the exepected from mapping. These are the correct
      // values when there is no perioidic bounds.
      for (int j = 0; j < total_dof; j++) {
        expected_from_map(j) = combined_sol(j);
      }

      // apply the periodic mapping by copying values from
      if (input_data_->boundaryType[0] == "periodic" &&
          input_data_->boundaryType[1] == "periodic") {  // periodic in X and Y
        int overwrite[7] = {12, 13, 14, 15, 3, 7, 11};
        int source[7] = {0, 1, 2, 0, 0, 4, 8};
        for (int i_idx = 0; i_idx < 7; i_idx++) {
          for (int i_unk = 0; i_unk < n_dof_; i_unk++) {
            // skip this if the variable is not periodic
            if (!input_data_->isPerVar[i_unk]) { continue; }
            int overwrite_idx = overwrite[i_idx] * n_dof_ + i_unk;
            int source_idx = source[i_idx] * n_dof_ + i_unk;
            expected_from_map(overwrite_idx) = expected_from_map(source_idx);
          }
        }
      } else if (input_data_->boundaryType[0] == "periodic") {  // periodic in X
        int overwrite[4] = {3, 7, 11, 15};
        int source[4] = {0, 4, 8, 12};
        for (int i_idx = 0; i_idx < 4; i_idx++) {
          for (int i_unk = 0; i_unk < n_dof_; i_unk++) {
            // skip this if the variable is not periodic
            if (!input_data_->isPerVar[i_unk]) { continue; }
            int overwrite_idx = overwrite[i_idx] * n_dof_ + i_unk;
            int source_idx = source[i_idx] * n_dof_ + i_unk;
            expected_from_map(overwrite_idx) = expected_from_map(source_idx);
          }
        }
      } else if (input_data_->boundaryType[1] == "periodic") {  // periodic in Y
        int overwrite[4] = {12, 13, 14, 15};
        int source[4] = {0, 1, 2, 3};
        for (int i_idx = 0; i_idx < 4; i_idx++) {
          for (int i_unk = 0; i_unk < n_dof_; i_unk++) {
            // skip this if the variable is not periodic
            if (!input_data_->isPerVar[i_unk]) { continue; }
            int overwrite_idx = overwrite[i_idx] * n_dof_ + i_unk;
            int source_idx = source[i_idx] * n_dof_ + i_unk;
            expected_from_map(overwrite_idx) = expected_from_map(source_idx);
          }
        }
      }

      // first make sure the global values match the expected
      for (int j = 0; j < total_dof; j++) {
        always_assert(expected_from_map(j) == combined_from(j));
      }

      // now make sure that each individual process has the correct mapping for
      // the unknowns that it has. It should not be possible for this to fail
      // if the above tests have passed, but check it just to be safe.
      for (int i_rank = 0; i_rank < mpi_size; i_rank++) {
        for (int j = 0; j < total_dof; j++) {
          // either this should be unknown beacause it's not on this process
          // or it should be equal to the expected value
          always_assert(from_all(i_rank, j) == kUnknown ||
                        expected_from_map(j) == from_all(i_rank, j));
        }
      }

      // if we've made it here, all the asserts above have passed and the
      // mapping is correct, so print a success message
      PrintStatusStream(std::cerr, "", "Mapping passed all checks.");
    }
  }

  void Integrands(const FEMElm& fe, ZeroMatrix<double>& Ae,
                  ZEROARRAY<double>& be) {}
};

#endif  // TESTS_PERIODICDDMAP_INCLUDE_PBC_EQUATION_H_
