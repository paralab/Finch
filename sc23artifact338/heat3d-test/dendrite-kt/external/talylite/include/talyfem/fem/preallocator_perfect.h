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
#pragma once

#include <talyfem/fem/preallocator.h>

#include <algorithm>
#include <set>
#include <vector>

#include <talyfem/grid/femelm.h>

#ifdef PETSC_USE_64BIT_INDICES
#define PETSC_INT_MAX LONG_LONG_MAX
#else
#define PETSC_INT_MAX INT_MAX
#endif

// prevents warnings when asserts are turned off
#define UNUSED(x) (void)(x)

// This might be a useful read if you're trying to understand this code:
// http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatMPIAIJSetPreallocation.html

namespace TALYFEMLIB {

/**
 * Repressents an (i, j) matrix index.
 * Used to give matrix indices an ordering (via the '<' operator) so
 * we can sort them.
 * Conveniently, sizeof(MatIndex) == sizeof(PetscInt * 2), so this class can be
 * transmitted with MPI as two PetscInts instead of a custom data type.
 */
class MatIndex {
 public:
  PetscInt row;  ///< row in the matrix (i)
  PetscInt col;  ///< column in the matrix (j)

  /**
   * Construct a MatIndex at (0, 0).
   */
  MatIndex() : row(0), col(0) {}

  /**
   * Construct a MatIndex at (r, c).
   * @param r row
   * @param c column
   */
  MatIndex(PetscInt r, PetscInt c) : row(r), col(c) {}

  /**
   * Check if lhs comes before rhs.
   * Used for sorting. Sorts by row, then breaks ties by column.
   * @param rhs other index
   * @returns true if lhs comes before rhs
   */
  bool operator<(const MatIndex& rhs) const {
    if (row < rhs.row)
      return true;
    if (row == rhs.row && col < rhs.col)
      return true;
    return false;
  }

  /**
   * Check if two indices are identical.
   * @param rhs other index
   * @returns true if lhs.row == rhs.row and lhs.col == rhs.col
   */
  bool operator==(const MatIndex& rhs) const {
    return (row == rhs.row && col == rhs.col);
  }

  /**
   * Check if two indices are not identical.
   * @param rhs other index
   * @returns true if lhs.row != rhs.row or lhs.col != rhs.col
   */
  bool operator!=(const MatIndex& rhs) const {
    return !(*this == rhs);
  }
};


/**
 * The "perfect" preallocator simulates a timestep of your problem to try and
 * determine exactly which parts of the Ae matrix can be non-zero.
 * In theory, this preallocates the smallest sparse matrix possible.
 *
 * To figure out which matrix indices can ever be non-zero, it calls
 * CEquation->IntegrandsPreallocator(elm, Ae) for every element, which you
 * should override to match your Integrands.
 */
template<class T>
class PreallocatorPerfect : public Preallocator {
 public:
  /**
   * Construct a PreallocatorPerfect.
   * eq must remain valid until calc_preallocation is called.
   * @param eq equation to preallocate for
   */
  explicit PreallocatorPerfect(CEquation<T>* eq) {
    row_start_ = -1;
    row_end_ = -1;
    equation_ = eq;
  }

  virtual void calc_preallocation(Mat& mat, PetscInt blockSize) override {
    const int mpi_rank = GetMPIRank();
    const int mpi_size = GetMPISize();

    GRID* grid = equation_->p_grid_;
    init_mat_ranges(mat);

    // stores nonzero global matrix indices found from this process's elements
    // these indices may be owned by another process
    std::vector< std::vector<MatIndex> > nonzero_indices(mpi_size);

    // initialize counts to zero
    diag_nonzero_count_.resize(local_row_count(), 0);
    offdiag_nonzero_count_.resize(local_row_count(), 0);

    // add boundary condition indices
    // this needs to be done whenever there's a boundary condition
    for (PetscInt row = row_start_; row < row_end_; row++) {
      nonzero_indices[mpi_rank].push_back(MatIndex(row, row));
    }

    // MPITimer calc_nonzeros_timer("Calculate nonzeros");
    // calc_nonzeros_timer.Start();

    // find all nonzeros caused by elements on this process
    FEMElm fe(equation_->p_grid_, 0);
    for (int elmID = 0; elmID < grid->n_elements(); elmID++) {
      if (!equation_->IsMyElement(elmID))
        continue;

      fe.refill(elmID, equation_->rel_order_of_itg());
      std::vector<MatIndex> indices = find_elem_nonzeros(fe);

      for (unsigned int i = 0; i < indices.size(); i++) {
        const MatIndex& idx = indices.at(i);
        nonzero_indices[get_row_owner(idx.row)].push_back(idx);
      }
    }

    for (int i = 0; i < mpi_size; i++) {
      std::sort(nonzero_indices[i].begin(), nonzero_indices[i].end());
      nonzero_indices[i].erase(std::unique(nonzero_indices[i].begin(),
                               nonzero_indices[i].end()),
                               nonzero_indices[i].end());
    }
    // calc_nonzeros_timer.Stop();
    // calc_nonzeros_timer.PrintTotalTimeSeconds();

    // Make sure there's no padding because that would break MPI_Allgather.
    // If this is ever false, to fix it you'll have to learn about custom
    // MPI data types for the MPI_Allgather, or just use an array of PetscInts.
    // This is a compile time check.
    static_assert(sizeof(MatIndex) == sizeof(PetscInt) * 2, "MatIndex size");

    for (int proc = 0; proc < mpi_size; proc++) {
      // gather all nonzeros we have for this process
      // MPITimer find_timer("Find timer");
      // find_timer.Start();
      std::vector<MatIndex>& indices = nonzero_indices[proc];
      // find_timer.Stop();

      // this must be an int to match the MPI_Gatherv sendcount type
      // multiply counts by 2 because we are counting ints, not MatIndexes
      int indices_size = static_cast<int>(indices.size() * 2);

      // gather how many indices each process will send into an array on proc
      std::vector<int> counts(mpi_size);
      MPI_Gather(&indices_size, 1, MPI_INT, &counts[0], 1,
                 MPI_INT, proc, PETSC_COMM_WORLD);

      if (proc == mpi_rank) {
        // calculate offsets into local_indices each process will send to
        std::vector<int> recv_offsets(mpi_size);
        recv_offsets[0] = 0;
        for (int i = 1; i < mpi_size; i++) {
          recv_offsets[i] = recv_offsets[i - 1] + counts[i - 1];
        }

        // divide by 2 because we were counting ints, not MatIndexes
        const int local_indices_total = (recv_offsets[mpi_size - 1] +
                                         counts[mpi_size - 1]) / 2;

        // receive
        std::vector<MatIndex> local_indices(local_indices_total);

        // MPITimer gather_timer("Gather timer");
        // gather_timer.Start();
        MPI_Gatherv(reinterpret_cast<int*>(&indices[0]), indices_size,
                    MPI_TALYFEM_INT, reinterpret_cast<int*>(&local_indices[0]),
                    &counts[0], &recv_offsets[0],
                    MPI_TALYFEM_INT, proc, PETSC_COMM_WORLD);
        // gather_timer.Stop();

        // remove duplicates
        // MPITimer unique_timer("Remove duplicates");
        // unique_timer.Start();
        std::sort(local_indices.begin(), local_indices.end());
        local_indices.erase(std::unique(local_indices.begin(),
                            local_indices.end()), local_indices.end());
        // unique_timer.Stop();
        // const int n_dupes = local_indices_total -
        //                     static_cast<int>(local_indices.size());
        // PetscSynchronizedPrintf(PETSC_COMM_WORLD, " proc %d removed %d "
        //                         "duplicates in %fs\n", mpi_rank, n_dupes,
        //                          unique_timer.GetTotalTimeSeconds());

        // process local indices
        // MPITimer process_timer("process timer");
        // process_timer.Start();
        for (unsigned int i = 0; i < local_indices.size(); i++) {
          const MatIndex& index = local_indices.at(i);
          MatOwnershipType type = get_ownership_type(index);
          // PetscSynchronizedPrintf(PETSC_COMM_WORLD, "  counting (%d, %d) as"
          //                        " %s\n", index.row, index.col,
          //                        type == OWNED_DIAG ? "diag" : "off-diag");

          if (type == OWNED_DIAG) {
            diag_nonzero_count_[index.row - row_start_]++;
          } else if (type == OWNED_OFFDIAG) {
            offdiag_nonzero_count_[index.row - row_start_]++;
          } else {
            std::cerr << "Preallocator warning: process " << proc <<
                         " received unowned index (" << local_indices[i].row
                         << ", " << local_indices[i].col << ")\n";
          }
        }
        // process_timer.Stop();
        // PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Find time: %f\n",
        //                         find_timer.GetTotalTimeSeconds());
        // PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Gather time: %f\n",
        //                         gather_timer.GetTotalTimeSeconds());
        // PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Process time: %f\n",
        //                         process_timer.GetTotalTimeSeconds());
      } else {
        // send indices to the process that owns them
        MPI_Gatherv(reinterpret_cast<int*>(&indices[0]), indices_size,
                    MPI_TALYFEM_INT,
                    NULL, 0, NULL, MPI_TALYFEM_INT, proc, PETSC_COMM_WORLD);
      }
      // PetscSynchronizedFlush(PETSC_COMM_WORLD);
    }

    // done, use get_dnz() and get_onz() for the results
  }

  const PetscInt* get_dnz() const override {
    return &diag_nonzero_count_.at(0);
  }

  const PetscInt* get_onz() const override {
    return &offdiag_nonzero_count_.at(0);
  }

  /**
   * Print out each process's calculated number of diagonal non-zeroes and
   * off-diagonal non-zeroes.
   * This is a collective call (all processors must call this function).
   * This is only useful after calc_preallocation(...) has been called.
   */
  void print_debug() const {
    const int mpi_rank = GetMPIRank();
    const int mpi_size = GetMPISize();

    PetscInt proc_total = 0;
    for (int proc = 0; proc < mpi_size; proc++) {
      if (proc == mpi_rank) {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD, "process %d:\n", proc);
        for (PetscInt i = 0; i < local_row_count(); i++) {
          PetscSynchronizedPrintf(PETSC_COMM_WORLD, "  dnz[%d] = %d,   "
                                  "onz[%d] = %d\n", row_start_ + i,
                                  get_dnz()[i], row_start_ + i, get_onz()[i]);
          proc_total += get_dnz()[i];
          proc_total += get_onz()[i];
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD, "  total for process %d: "
                                "%d\n", mpi_rank, proc_total);
      }
    }
    PetscSynchronizedFlush(PETSC_COMM_WORLD, stdout);

    PetscInt total = 0;
    MPI_Reduce(&proc_total, &total, 1, MPI_TALYFEM_INT,
               MPI_SUM, 0, PETSC_COMM_WORLD);

    PrintInfo("Total allocated by perfect preallocator: ", total);
  }

 protected:
  /**
   * Used to indicate if a matrix index is owned by us, and if it is,
   * if it's on the diagonal or off the diagonal.
   */
  enum MatOwnershipType {
    NOT_OWNED,     ///< not owned by this process
    OWNED_DIAG,    ///< owned by us, on the diagonal
    OWNED_OFFDIAG  ///< owned by us, off the diagonal
  };

  /**
   * Figure out if a given matrix index is owned by us, and if it's on
   * the diagonal or off-diagonal.
   * NOTE: You must call init_mat_ranges(...) first for this to be valid.
   * @param row matrix row (i)
   * @param col matrix column (j)
   * @returns this process's ownership of (row, col)
   */
  MatOwnershipType get_ownership_type(PetscInt row, PetscInt col) const {
    if (row < row_start_ || row >= row_end_)
      return NOT_OWNED;
    if (col < row_start_ || col >= row_end_)
      return OWNED_OFFDIAG;
    return OWNED_DIAG;
  }

  /**
   * Figure out if a given MatIndex is owned by us, and if it's on
   * the diagonal or off-diagonal.
   * @param idx matrix index to get ownership for
   * @returns this process's ownership of idx
   */
  inline MatOwnershipType get_ownership_type(const MatIndex& idx) const {
    return get_ownership_type(idx.row, idx.col);
  }

  /**
   * Return which process owns the given row.
   * NOTE: You must call init_mat_ranges(...) first for this to be valid.
   * @param row row to get owner for
   * @returns which process owns row
   */
  int get_row_owner(PetscInt row) const {
    return row_owners_.at(row);
  }

  /**
   * Get how many rows this process owns.
   * Also happens to be the length of diag_nonzero_count_ and
   * offdiag_nonzero_count_.
   * NOTE: You must call init_mat_ranges(...) first for this to be valid.
   * @returns how many rows this process own
   */
  inline PetscInt local_row_count() const {
    return row_end_ - row_start_;
  }

  /**
   * Initializes information we need to calculate which processes own what
   * in the matrix.
   * @param mat PETSc matrix to initialize with
   */
  void init_mat_ranges(const Mat& mat) {
    // we can't use MatGetOwnershipRange because that only works after the
    // matrix has been preallocated.

    // If CEquation used PETSC_DECIDE to pick which rows went to which process,
    // we would use PetscSplitOwnership().
    // This is the technique recommended in the MatGetOwnershipRange()
    // documentation.
    // But CEquation actually uses grid->n_owned_nodes() * n_dof to decide
    // how many local rows there are.

    PetscErrorCode ierr;
    UNUSED(ierr);  // prevent unused warnings when asserts are off

    /*
    // PETSC_DECIDE method to find owned_rows
    PetscInt total_rows, total_cols;
    ierr = MatGetSize(mat, &total_rows, &total_cols);
    assert(ierr == 0);

    // Calculate the number of rows this process will own after preallocation.
    PetscInt this_owned_rows = PETSC_DECIDE;
    ierr = PetscSplitOwnership(PETSC_COMM_WORLD, &this_owned_rows, &total_rows);
    assert(ierr == 0);
    */

    // grid method to find number of owned rows for this process
    PetscInt this_owned_rows =
        equation_->p_grid_->n_owned_nodes() * equation_->n_dof();

    const int mpi_rank = GetMPIRank();
    const int mpi_size = GetMPISize();

    std::vector<PetscInt> owned_rows(mpi_size);
    ierr = MPI_Allgather(&this_owned_rows, 1, MPI_TALYFEM_INT,
                         &owned_rows[0], 1, MPI_TALYFEM_INT, PETSC_COMM_WORLD);
    assert(ierr == MPI_SUCCESS);

    unsigned int total_rows = 0;
    for (int i = 0; i < mpi_size; i++) {
      total_rows += owned_rows[i];
    }
    row_owners_.resize(total_rows);

    int row = 0;
    for (int proc = 0; proc < mpi_size; proc++) {
      if (proc == mpi_rank)
        row_start_ = row;
      for (int i = 0; i < owned_rows[proc]; i++) {
        row_owners_[row] = proc;
        row++;
      }
      if (proc == mpi_rank)
        row_end_ = row;
    }
  }

  /**
   * Returns the list of non-zero (global) matrix indices for a given element.
   * @param fe element to get indices for
   * @returns non-zero (global) matrix indices for fe
   */
  virtual std::vector<MatIndex> find_elem_nonzeros(FEMElm& fe) const {
    ZeroMatrix<bool> Ae;
    Ae.redim(equation_->n_rows(fe), equation_->n_cols(fe));
    Ae.fill(0);

    // mimic eq->Assemble
    while (fe.next_itg_pt()) {
      equation_->IntegrandsPreallocator(fe, Ae);
    }

    // calculate index mappings for Ae into the big matrix
    ZEROARRAY<PetscInt> rows, cols;
    equation_->CalcAebeIndices(fe, rows, cols);

    // find nonzeros and return their corresponding "global" matrix indices
    std::vector<MatIndex> nz;

    for (int row = 0; row < Ae.nx(); row++) {
      for (int col = 0; col < Ae.ny(); col++) {
        if (Ae(row, col) == true && rows(row) != -1 && cols(col) != -1) {
          nz.push_back(MatIndex(rows(row), cols(col)));
        }
      }
    }

    return nz;
  }

  // see init_mat_ranges for where these are set
  PetscInt row_start_;  ///< our first owned row
  PetscInt row_end_;  ///< our last owned row + 1
  std::vector<int> row_owners_;  ///< row_owners_[i] == process that owns row i

  std::vector<PetscInt> diag_nonzero_count_;  ///< number of digaonal non-zeros
  std::vector<PetscInt> offdiag_nonzero_count_;  ///< off-diagonal non-zeros

  CEquation<T>* equation_;  ///< kept so we can call IntegrandsPreallocator
};

}  // namespace TALYFEMLIB
