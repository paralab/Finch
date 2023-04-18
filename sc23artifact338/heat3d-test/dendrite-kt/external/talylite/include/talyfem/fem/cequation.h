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
#ifndef FEM_CEQUATION_H_
#define FEM_CEQUATION_H_

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petsc.h>
#else

#include <petsc.h>

#endif


#include <time.h>

#include <vector>  // for std::vector
#include <memory>  // for std::shared_ptr

#include <talyfem/common/petsc_logging.h>  // for PetscEventLogger
#include <talyfem/data_structures/zeroarray.h>
#include <talyfem/data_structures/zeromatrix.h>
#include <talyfem/fem/boundary_conditions.h>  // class of member variable
#include <talyfem/fem/periodic_bounds.h>  // class of member variable
#include <talyfem/fem/periodic_data.h>  // used in redim
#include <talyfem/fem/periodic_exchanger.h>  // class of member variable
#include <talyfem/fem/preallocator.h>  // class of member variable
#include <talyfem/fem/preallocator_original.h>  // default preallocator
#include <talyfem/grid/grid_common.h>
#include <talyfem/grid/grid-types.h>
#include <talyfem/grid/grid_types/grid.h>  // for GRID class
#include <talyfem/grid/gridfield.h>  // class of member variable
#include <talyfem/grid/nodeid_types.h>
#include <talyfem/utils/utils.h>  // for PrintError()
#include <talyfem/utils/timers.h>  // for MPITimer


namespace TALYFEMLIB {

/**
 * Method of assembly
 */
    enum AssemblyMethod {
        kAssembleGaussPoints = 0,  // use Gauss point based assembly
        kAssembleElements = 1,  // use element based assembly
    };

/**
 * This is the base class for an equation that is to be solved with the library.
 *
 * Each equation to be solved should derive from this class.
 *
 * This class is responsible for the assembly and solution of the finite
 * element matrices used to the solve the problem. The derived class will need
 * to provide implementations for the following:
 *  fillEssBC - to fill essential boundary conditions (if needed)
 *  Integrands - to assemble the matrix values for each Gauss point
 *  Integrands4side - to assemble values for boundary conditions (if needed)
 *  SetIC - to set the initial conditions
 *  Solve - to call the solver.
 *
 * Several examples of each of these are shown in the tutorials.
 */
    template<class NodeData>
    class CEquation {
    public:
        /**
         * Construct the equation object
         *
         * @param has_uniform_mesh whether all the elements are the same size and
         *                          shape (this triggers some optimization)
         * @param equation_assembly_method which method of assembly to use for the
         *                                 solver (either kAssembleGaussPoints or
         *                                 kAssembleElements)
         */
        explicit CEquation(bool has_uniform_mesh = false,
                           AssemblyMethod equation_assembly_method = kAssembleGaussPoints) {
            has_uniform_mesh_ = has_uniform_mesh;
            assembly_method_ = equation_assembly_method;
            is_serial_ = false;
            p_grid_ = NULL;
            n_total_dof_ = 0;
            p_Ag_ = NULL;
            p_bg_ = NULL;
            p_x_map_ = NULL;
            created_ = false;

            is_periodic_ = false;
            exchanger_ = NULL;

            this->initEssBC();
            this->initPerBC();

            rel_order_ = 0;  // relative order of integration
            basis_flags_ = BASIS_DEFAULT;
            hermite_n_dof_scale_ = 1;
            preallocator_ = NULL;
            matPreallocationDone_ = false;
            preallocator_is_owned_ = false;

            recalc_matrix_ = true;
        }

        virtual ~CEquation() {
            Destroy();
        }

        /**
         * Initialize the equation solver on the given grid.
         *
         * This allocates the memory needed for the solver and sets up the arrays and
         * matrices used by the class.
         *
         * @param p_grid pointer of the grid to use
         * @param n_dof_value number of degrees of freedom per node
         * @param hermite whether or not we should use hermite basis functions
         *                (note this requires a GRID of linear elements)
         * @param rel_order relative order of integration (affects # of gauss points),
         *                  usually 0
         * @param periodic_data periodic data object
         */
        virtual void redimSolver(GRID *p_grid, int n_dof_value,
                                 bool hermite, int rel_order,
                                 PeriodicData *periodic_data = NULL) {
            this->p_grid_ = p_grid;

            if (hermite) {
                hermite_n_dof_scale_ = (int) pow(2, p_grid_->nsd());
            } else {
                hermite_n_dof_scale_ = 1;
            }
            rel_order_ = rel_order;

            this->n_dof_ = n_dof_value;
            n_total_dof_ = p_grid->n_nodes() * n_dof();
            solution_.redim(n_total_dof_);
            solution_.fill(0);
            setGlobalAb(&Ag_, &bg_, p_grid->GetSolutionMap());
            Create();
            if (has_uniform_mesh_) { InitializeAcceleration(); }
            if (periodic_data != NULL) {
                InitializePeriodicData(periodic_data);
            }
        }

        /**
         * Initialize the equation solver on the given grid.
         *
         * This is is an overload left in for compatibility reasons.
         * The function signature matches the old implementation of redimSolver,
         * (when BF was a constant for the entire grid).
         *
         * @param p_grid pointer of the grid to use
         * @param n_dof_value number of degrees of freedom per node
         * @param bf basis function to use (only significant for BASIS_HERMITE)
         * @param rel_order relative order of integration (affects # of gauss points),
         *                  usually 0
         * @param periodic_data periodic data object
         */
        inline void redimSolver(GRID *p_grid, int n_dof_value,
                                kBasisFunction bf, int rel_order,
                                PeriodicData *periodic_data = NULL) {
            redimSolver(p_grid, n_dof_value, bf == BASIS_HERMITE,
                        rel_order, periodic_data);
        }

        /**
         * Precalculate basis arrarys to speed up calculations
         *
         * This is only valid when all of the elements are the same size and shape.
         */
        void InitializeAcceleration() {
            FEMElm fe(p_grid_, basis_flags_);

            if (is_hermite()) {
                fe.refill(0, BASIS_HERMITE, rel_order_);
            } else {
                fe.refill(0, rel_order_);
            }

            int nsd = p_grid_->nsd();
            int n_basis_functions = fe.nbf();
            int n_gauss_points = fe.n_itg_pts();

            // initialize feAccelerate_ with copies of fe
            feAccelerate_.resize(n_gauss_points, fe);
            dNdNArray_.redim(n_gauss_points);
            sumdNdNArray_.redim(n_gauss_points);
            NNArray_.redim(n_gauss_points);
            NdNArray_.redim(n_gauss_points);

            while (fe.next_itg_pt()) {
                const int gauss_point_id = fe.cur_itg_pt_num();
                for (int i = 0; i < n_gauss_points; i++) {
                    if (i >= gauss_point_id) {
                        feAccelerate_.at(i).next_itg_pt();
                    }
                }
                dNdNArray_(gauss_point_id).redim(n_basis_functions, n_basis_functions);
                sumdNdNArray_(gauss_point_id).redim(n_basis_functions, n_basis_functions);
                NNArray_(gauss_point_id).redim(n_basis_functions, n_basis_functions);
                NdNArray_(gauss_point_id).redim(n_basis_functions, n_basis_functions);
                for (int a = 0; a < n_basis_functions; a++) {
                    for (int b = 0; b < n_basis_functions; b++) {
                        dNdNArray_(gauss_point_id)(a, b).redim(nsd, nsd);
                        NdNArray_(gauss_point_id)(a, b).redim(nsd);
                        sumdNdNArray_(gauss_point_id)(a, b) = 0;
                        for (int i = 0; i < nsd; i++) {
                            for (int j = 0; j < nsd; j++) {
                                dNdNArray_(gauss_point_id)(a, b)(i, j) =
                                        fe.dN(a, i) * fe.dN(b, j);
                            }
                            sumdNdNArray_(gauss_point_id)(a, b) += fe.dN(a, i) * fe.dN(b, i);
                            NdNArray_(gauss_point_id)(a, b)(i) = fe.N(a) * fe.dN(b, i);
                        }
                        NNArray_(gauss_point_id)(a, b) = fe.N(a) * fe.N(b);
                    }
                }
            }
            // precalculate simple summation values for elemental assembly
            if (assembly_method_ == kAssembleElements) {
                // set up memory
                ea_NN_.redim(n_basis_functions, n_basis_functions);
                ea_NN_.fill(0.0);
                ea_dNdN_.redim(n_basis_functions, n_basis_functions);
                ea_dNdN_.fill(0.0);
                ea_NNdN_gauss_sum_.redim(n_basis_functions, n_basis_functions);
                for (int a = 0; a < n_basis_functions; a++) {
                    for (int b = 0; b < n_basis_functions; b++) {
                        ea_NNdN_gauss_sum_(a, b).redim(n_basis_functions);
                        ea_NNdN_gauss_sum_(a, b).fill(0.0);
                    }
                }

                for (int g = 0; g < n_gauss_points; g++) {
                    const double detJxW = feAccelerate_.at(g).detJxW();

                    for (int a = 0; a < n_basis_functions; a++) {
                        const double M = NNArray_(g)(a, a) * detJxW;
                        const double K = sumdNdNArray_(g)(a, a) * detJxW;
                        ea_NN_(a, a) += M;
                        ea_dNdN_(a, a) += K;
                        for (int c = 0; c < n_basis_functions; c++) {
                            ea_NNdN_gauss_sum_(a, a)(c) += M * feAccelerate_.at(g).N(c);
                        }
                    }

                    for (int a = 0; a < n_basis_functions - 1; a++) {
                        for (int b = a + 1; b < n_basis_functions; b++) {
                            const double M = NNArray_(g)(a, b) * detJxW;
                            const double K = sumdNdNArray_(g)(a, b) * detJxW;
                            ea_NN_(a, b) += M;
                            ea_dNdN_(a, b) += K;
                            ea_NN_(b, a) += M;
                            ea_dNdN_(b, a) += K;
                            for (int c = 0; c < n_basis_functions; c++) {
                                ea_NNdN_gauss_sum_(a, b)(c) += M * feAccelerate_.at(g).N(c);
                                ea_NNdN_gauss_sum_(b, a)(c) += M * feAccelerate_.at(g).N(c);
                            }
                        }
                    }
                }
            }
        }

        /**
         * Set the grid field data to be used by the solver
         *
         * @param p_data pointer to the grid field data
         */
        virtual void setData(GridField<NodeData> *p_data) {
            this->p_data_ = p_data;
        }

        /**
         * Fill essential boundary conditions
         *
         * This function needs to be implemented for any solver that uses essential
         * boundary conditions. It should use the specifyValue function to set the
         * essential boundary conditions for each node where they apply.
         * Implementations of this function will commonly use the p_grid_->BoNode
         * function to determine which nodes are on a boundary.
         *
         * Prior to calling this function, initEssBC() should be called to remove
         * any previously defined boundaries.
         *
         * The tutorials have several example implementations of this function.
         */
        virtual void fillEssBC() {}

        /**
         * Set the initial conditions for the system.
         *
         * This function needs to be implemented for any solver in order to set the
         * values of the initial conditions for the system.
         */
        virtual void SetIC() {}

        /**
         * Tell the solver information about the periodic structure of the system
         *
         * Since the periodic boundary conditions are not set until after the
         * Preallocate function is called, Preallocate is not aware of the system's
         * periodic structure. Without that information, it cannot accurately
         * determine memory usage. This function should be called prior to
         * redim in order to ensure Preallocate knows about the periodic bounds.
         *
         * @param n_per_bounds number of periodic boundaries
         * @param n_per_vars number of periodic variables per node
         * @param n_total_vars number of total variables per node
         */
        virtual void PresetPeriodicData(int n_per_bounds, int n_per_vars,
                                        int n_total_vars) {
            if (preallocator_ == NULL) {
                throw TALYException() <<
                                      "The preallocator must be set before calling PresetPeriodicData!";
            }
            preallocator_->PresetPeriodicData(n_per_bounds, n_per_vars, n_total_vars);
        }

        /**
         * Sets the preallocator for this solver to the given value
         *
         * This can be used to specify a preallocator other than the default one.
         * If the no preallocator is specified, or if the value is NULL, the default
         * preallocation will be performed.
         *
         * @param prealloc pointer to the preallocator object to use
         */
        void SetPreallocator(Preallocator *prealloc = NULL) {
            if (preallocator_ != NULL && prealloc == NULL) {
                PrintWarning("Ignoring SetPreallocator"
                             " (a preallocator was already set).");
                return;
            }

            if (prealloc == NULL) {
                assert(p_grid_ != NULL);  // if redimSolver wasn't called this can happen
                PreallocatorOriginal *prealloc_orig = new PreallocatorOriginal();
                prealloc_orig->redim(p_grid_, n_dof());
                preallocator_ = prealloc_orig;
                preallocator_is_owned_ = true;
            } else {
                preallocator_ = prealloc;
            }
        }

        /**
         * Allocates memory required for the solver object.
         *
         * Creates the PETSc matrix, but does not preallocate memory for it until
         * UpdateMatPreallocation() is called.
         * @return 0 on success or PETSc error code if there are any errors
         */
        virtual PetscErrorCode Create() {
            if (created_) Destroy();

            PetscErrorCode ierr;

            const PetscInt globalRows = p_grid_->n_total_nodes() * n_dof();
            const PetscInt globalCols = p_grid_->n_total_nodes() * n_dof();

            // create the matrix
            if (p_grid_->parallel_type_ == kNoDomainDecomp) {
                // set matrix in the way suggested by Petsc 3.3+
                ierr = MatCreate(is_serial_ ? PETSC_COMM_SELF : PETSC_COMM_WORLD, &Ag_);
                CHKERRQ(ierr);
                ierr = MatSetSizes(Ag_, PETSC_DECIDE, PETSC_DECIDE,
                                   globalRows, globalCols);
                CHKERRQ(ierr);
                ierr = MatSetFromOptions(Ag_);
                CHKERRQ(ierr);
                ierr = VecCreateMPI(is_serial_ ? PETSC_COMM_SELF : PETSC_COMM_WORLD,
                                    PETSC_DECIDE, globalRows, &bg_);
                CHKERRQ(ierr);
            } else {
                const int myNodesNumber = p_grid_->n_owned_nodes();

                // set matrix in the way suggested by Petsc 3.3+
                ierr = MatCreate(PETSC_COMM_WORLD, &Ag_);
                CHKERRQ(ierr);
                ierr = MatSetSizes(Ag_, myNodesNumber * n_dof(), myNodesNumber * n_dof(),
                                   globalRows, globalCols);
                CHKERRQ(ierr);
                ierr = MatSetFromOptions(Ag_);
                CHKERRQ(ierr);

                // set solution vector in the way suggested by Petsc 3.3+
                ierr = VecCreate(PETSC_COMM_WORLD, &bg_);
                CHKERRQ(ierr);
                ierr = VecSetSizes(bg_, myNodesNumber * n_dof(),
                                   p_grid_->n_total_nodes() * n_dof());
                CHKERRQ(ierr);
                ierr = VecSetBlockSize(bg_, n_dof());
                CHKERRQ(ierr);
                ierr = VecSetFromOptions(bg_);
                CHKERRQ(ierr);
            }

            ierr = VecDuplicate(bg_, &xg_);
            CHKERRQ(ierr);
            if (!is_serial_) {
                KSPCreate(PETSC_COMM_WORLD, &ksp_);
            } else {
                KSPCreate(PETSC_COMM_SELF, &ksp_);
            }

            p_Ag_ = &Ag_;
            p_bg_ = &bg_;
            p_x_map_ = NULL;
            created_ = true;

            has_ess_bc_.redim(n_total_dof_);
            has_ess_bc_.fill(false);
            has_per_bc_.redim(n_total_dof_);
            has_per_bc_.fill(false);
            is_node_periodic_.redim(p_grid_->n_nodes());
            is_node_periodic_.fill(false);
            if (p_grid_->parallel_type_ == kWithDomainDecomp) {
                p_x_map_ = &x_map_;
                x_map_.redim(n_total_dof_);
                idx_from_.redim(n_total_dof_);
                idx_to_.redim(n_total_dof_);
                for (LocalNodeID n = 0; n < p_grid_->n_nodes(); n++) {
                    for (int i = 0; i < n_dof(); i++) {
                        LocalVarIdx lid = (n) * n_dof() + i;
                        GlobalVarIdx gid = (p_grid_->solution_map(n)) * n_dof() + i;
                        idx_from_(lid) = gid;
                        idx_to_(lid) = lid;
                        x_map_(lid) = gid;
                    }
                }
                ISCreateGeneral(PETSC_COMM_SELF, n_total_dof_, idx_from_.data(),
                                PETSC_COPY_VALUES, &from_);
                ISCreateGeneral(PETSC_COMM_SELF, n_total_dof_, idx_to_.data(),
                                PETSC_COPY_VALUES, &to_);
            }
            return 0;
        }

        /**
         * Update the preallocation of the matrix.
         * If this is not done before the matrix is accessed for assembly, everything
         * will slow to a crawl as PETSc constantly re-allocates the entire matrix
         * for every unexpected nonzero and you will get gigabytes of warning logs.
         *
         * This is not automatically called in Create() because the ess/per boundary
         * condition data may not be set yet.
         * @returns 0 on success, or PETSc error code if there was a problem
         */
        PetscErrorCode UpdateMatPreallocation() {
            if (matPreallocationDone_)
                return 0;

            PetscEventLogger ev_prealloc("Preallocation");

            PetscInt bs = -1;
            // PetscOptionsGetInt(NULL, "-mat_block_size", &bs, NULL);

            MPITimer timer("total preallocation time");
            timer.Start();

            PetscErrorCode ierr;
            if (p_grid_->parallel_type_ == kNoDomainDecomp) {
                // this approximation kind of sucks
                int AroundNeighbor = static_cast<int>(pow(3.0f, p_grid_->nsd()));
                AroundNeighbor *= (3 * this->p_grid_->basis_order() - 2);
                if (p_grid_->nsd() == 3 && this->p_grid_->basis_order() > 1) {
                    AroundNeighbor *= this->p_grid_->basis_order();
                }
                if (p_grid_->grid_type() == kGrid2dTriangle) {
                    AroundNeighbor += 1;
                }

                if (bs == -1) {
                    ierr = MatSeqAIJSetPreallocation(Ag_, AroundNeighbor * n_dof(),
                                                     PETSC_NULL);
                    CHKERRQ(ierr);
                    ierr = MatMPIAIJSetPreallocation(Ag_, AroundNeighbor * n_dof(),
                                                     PETSC_NULL, AroundNeighbor * n_dof(),
                                                     PETSC_NULL);
                    CHKERRQ(ierr);
                } else {
                    ierr = MatSeqBAIJSetPreallocation(Ag_, bs,
                                                      AroundNeighbor * n_dof() / bs,
                                                      PETSC_NULL);
                    CHKERRQ(ierr);
                    ierr = MatMPIBAIJSetPreallocation(Ag_, bs,
                                                      AroundNeighbor * n_dof() / bs,
                                                      PETSC_NULL,
                                                      AroundNeighbor * n_dof() / bs,
                                                      PETSC_NULL);
                    CHKERRQ(ierr);
                }
            } else {
                // currently, we only use the preallocator in the DD case
                if (preallocator_ == NULL) {
                    SetPreallocator();
                    if (is_periodic_) {
                        PresetPeriodicData(
                                periodic_data_->periodic_bounds()->n_periodic_bounds(),
                                periodic_data_->NumPeriodicVars(), n_dof_);
                    }
                }
                preallocator_->calc_preallocation(Ag_, bs);

                if (bs == -1) {
                    // use normal sequential AIJ format
                    ierr = MatSeqAIJSetPreallocation(Ag_, 0, preallocator_->get_dnz());
                    CHKERRQ(ierr);
                    ierr = MatMPIAIJSetPreallocation(Ag_, 0, preallocator_->get_dnz(), 0,
                                                     preallocator_->get_onz());
                    CHKERRQ(ierr);
                } else {
                    // use block compressed row format
                    ierr = MatSeqBAIJSetPreallocation(Ag_, bs, 0, preallocator_->get_dnz());
                    CHKERRQ(ierr);
                    ierr = MatMPIBAIJSetPreallocation(Ag_, bs, 0, preallocator_->get_dnz(),
                                                      0, preallocator_->get_onz());
                    CHKERRQ(ierr);
                }
            }

            // treat new nonzeros that would cause a PetscMalloc as errors
            // (as this means there was a problem with preallocation)
            MatSetOption(Ag_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
            // Not sure we want to do this. This will ignore legitamately
            // calculated zero values. That may or may not be a problem
            // MatSetOption(Ag_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);

            matPreallocationDone_ = true;

            timer.Stop();
            timer.PrintTotalTimeSeconds();
            return 0;
        }

        /**
         * Frees all of the solver memory
         *
         * @return 0 on success or PETSc error code if there are any errors
         */
        virtual PetscErrorCode Destroy() {
            if (!created_) { return 0; }
            PetscErrorCode ierr;
            ierr = KSPDestroy(&ksp_);
            CHKERRQ(ierr);
            ierr = MatDestroy(&Ag_);
            CHKERRQ(ierr);
            ierr = VecDestroy(&bg_);
            CHKERRQ(ierr);
            ierr = VecDestroy(&xg_);
            CHKERRQ(ierr);
            if (p_grid_ != NULL && p_grid_->parallel_type_ == kWithDomainDecomp) {
                ierr = ISDestroy(&from_);
                CHKERRQ(ierr);
                ierr = ISDestroy(&to_);
                CHKERRQ(ierr);
            }
            created_ = false;
            matPreallocationDone_ = false;

            this->periodic_vars_.cleanup();

            if (exchanger_ != NULL) {
                delete exchanger_;
            }
            if (preallocator_is_owned_) {
                delete preallocator_;
            }
            preallocator_ = NULL;

            boundary_conditions_.DeleteAllConditions();

            return 0;
        }

        /**
         * Initializes essential boundary conditions by removing any previous values.
         *
         * This should always be called prior to filling the essential boundary
         * conditions. This will NOT remove essential conditions created from
         * periodic bounds. To remove them, call initPerBC().
         */
        void initEssBC() {
            boundary_conditions_.DeleteDirichletConditions();
            has_ess_bc_.fill(false);
        }

        /**
         * Initializes the periodic boundary conditions
         *
         * This should always be called prior to filling periodic boundary
         * conditions.
         */
        void initPerBC() {
            boundary_conditions_.DeletePeriodicConditions();
            has_per_bc_.fill(false);
            is_node_periodic_.fill(false);
        }

        /**
         * Specifies an essential boundary condition for an unknown
         *
         * @param nodeID which node we are adding a boundary condition to (start at 0)
         * @param dir which degree of freedom this will apply to (start at 0)
         * @param val the value of the essential boundary condition
         * @param coeff The value of the diagonal term A(*,*). (Default value is 1)
         *              Some times, A(*,*)=1 may make the system ill-conditioned
         *              because of this essential bc. At that case, you can change
         *              the value of coeff to some other value. Usually you don't
         *              have to specify the value of the diagonal iter. 1 is fine for
         *              most cases.
         */
        void specifyValue(LocalNodeID nodeID, int dir, double val = 0,
                          double coeff = 1) {
            LocalVarIdx index = n_dof() * (nodeID) + dir;
            boundary_conditions_.AddDirichletCondition(index, coeff, val * coeff);

            // NOTE: for domain decomposition, this may lead to all processes that share
            // a particular node trying to simultaneously apply a boundary condition.
            // This probably doesn't matter as long as they are being set to the same
            // value on all processes. Otherwise, this could lead to incorrect results.
            has_ess_bc_.set(index, true);
        }

        /**
         * Specifies a periodic boundary condition for an unknown
         *
         * @param nodeID which node we are adding a boundary condition to (start at 0)
         * @param dir which degree of freedom this will apply to (start at 0)
         * @param coeff this is not used
         * TODO: remove useless coeff param
         */
        void specifyPerValue(LocalNodeID nodeID, int dir, double coeff = 1) {
            LocalVarIdx index = n_dof() * (nodeID) + dir;
            boundary_conditions_.AddPeriodicCondition(index);
            has_per_bc_.set(index, true);  // value in correct place
            is_node_periodic_.set(nodeID, true);
        }

        /**
         * Prints all boundary conditions to stdout
         */
        void printEssBC() const {
            boundary_conditions_.PrintAllConditions();
        }

        /**
         * Uses PETSc's KSP Solver to solve the system.
         *
         * This function is typically called from the Solve function implemented by
         * a derived class.
         *
         * @param solution The vector to store the solution.
         * @param print_timings whether to print timings
         * @return 0 on success or PETSc error code if there are any errors
         *
         * TODO: consider changing timings to use timing code and Print* functions
         */
        PetscErrorCode SolveKSP(ZEROARRAY<double> &solution, bool print_timings = false) {
            PetscErrorCode ierr;
            time_t tm1;
            time(&tm1);

            KSPSetOperators(ksp_, Ag_, Ag_);
            PC pc;
            KSPGetPC(ksp_, &pc);
            PCSetReusePreconditioner(pc, recalc_matrix_ ? PETSC_FALSE : PETSC_TRUE);

            KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);
            VecZeroEntries(xg_);

            KSPSetFromOptions(ksp_);
            KSPSetUp(ksp_);
            KSPSetTolerances(ksp_, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT,
                             PETSC_DEFAULT);

            time_t tm2;
            time(&tm2);
            if (p_grid_->grid_id() == 1 && print_timings) {
                printf("(%3d)", static_cast<int>(difftime(tm2, tm1)));
            }
            ierr = KSPSolve(ksp_, bg_, xg_);
            time_t tm3;
            time(&tm3);
            if (p_grid_->grid_id() == 1 && print_timings) {
                printf("(%3d)", static_cast<int>(difftime(tm3, tm2)));
            }

            if (p_grid_->parallel_type_ == kWithDomainDecomp) {
                VecScatter scatter;
                Vec SolutionVec;
                ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, n_total_dof_,
                                             solution.data(), &SolutionVec);
                CHKERRQ(ierr);
                ierr = VecScatterCreate(xg_, from_, SolutionVec, to_, &scatter);
                CHKERRQ(ierr);
                ierr = VecScatterBegin(scatter, xg_, SolutionVec, INSERT_VALUES,
                                       SCATTER_FORWARD);
                CHKERRQ(ierr);
                ierr = VecScatterEnd(scatter, xg_, SolutionVec, INSERT_VALUES,
                                     SCATTER_FORWARD);
                CHKERRQ(ierr);
                ierr = VecDestroy(&SolutionVec);
                CHKERRQ(ierr);
                ierr = VecScatterDestroy(&scatter);
                CHKERRQ(ierr);
            } else {
                VecScatter scatter;
                Vec solution_vec;
                ierr = VecScatterCreateToAll(xg_, &scatter, &solution_vec);
                CHKERRQ(ierr);
                ierr = VecScatterBegin(scatter, xg_, solution_vec, INSERT_VALUES,
                                       SCATTER_FORWARD);
                CHKERRQ(ierr);
                ierr = VecScatterEnd(scatter, xg_, solution_vec, INSERT_VALUES,
                                     SCATTER_FORWARD);
                CHKERRQ(ierr);
                double *array;
                ierr = VecGetArray(solution_vec, &array);
                CHKERRQ(ierr);
                memcpy(solution.data(), array, sizeof(double) * n_total_dof_);
                ierr = VecRestoreArray(solution_vec, &array);
                CHKERRQ(ierr);
                ierr = VecScatterDestroy(&scatter);
                CHKERRQ(ierr);
                ierr = VecDestroy(&solution_vec);
                CHKERRQ(ierr);
                // map values across periodic bounds
                if (is_periodic_) { CopyPeriodicSolution(); }
            }
            time_t tm4;
            time(&tm4);
            if (p_grid_->grid_id() == 1 && print_timings) {
                printf("(%3d)", static_cast<int>(difftime(tm4, tm3)));
            }

            return 0;
        }

        /**
         * Uses PETSc's KSP Solver to solve the system.
         * This is an overload included to match the old signature that had an
         * unused "iter" argument. You should use the simpler signature.
         *
         * @param solution
         * @param print
         * @param iter ignored
         * @returns 0 if successful
         */
        inline PetscErrorCode SolveKSP(ZEROARRAY<double> &solution, int iter, int print) {
            return SolveKSP(solution, print == 1);
        }

        /**
         * Creates the Ax=b system using the finite element method
         *
         * This function assembles the matrix and vector used to solve the system.
         * In addition, it applies the previously defined essentially boundary
         * conditions to the computation matrix and vector.
         */
        virtual void makeSystem() {
            ApplyEssBCToSolution();
            Assemble();
            ApplyEssBC();
        }

        /**
         * Assembles the linear matrix system
         *
         * This function ultimately calls Integrands and Integrands4side to
         * calculate the matrix and vector elements for the matrix system.
         * Following the construction of the values for each element, the data is
         * sent to PETSc and assembled into the global matrix.
         *
         * @param assemble_surface whether to do the surface assembly. Setting this
         *                         to false when there are no surface conditions will
         *                         result in a slight performance improvement.
         * @throw TALYException on error
         *
         * TODO: there is a petsc call for zeroing bg_ that is better than the method
         *       used here.
         * TODO: is assemble surface actually used? this is also done with the
         *       function AssembleElementSurface called from AssembleElement
         */
        virtual void Assemble(bool assemble_surface = true) {
            PetscEventLogger ev_assemble("Assemble");
            PetscPushErrorHandler(&PetscExceptionErrorHandler, NULL);

            try {
                PetscErrorCode ierr;
                ierr = UpdateMatPreallocation();
                CHKERRV(ierr);

                // zero the stiffness matrix and load
                if (recalc_matrix_) { MatZeroEntries(Ag_); }

                // trying to do VecZeroEntries (bg_);
                PetscInt nlocal;
                double *array;
                VecGetLocalSize(bg_, &nlocal);
                VecGetArray(bg_, &array);
                memset(array, 0, sizeof(double) * nlocal);
                VecRestoreArray(bg_, &array);

                if (has_uniform_mesh_) {
                    AssembleVolumeUniformMesh(assemble_surface, basis_flags_);
                } else {
                    AssembleVolume(assemble_surface, basis_flags_);
                }
                if (assemble_surface) {
                    AssembleSurface();
                }
                if (recalc_matrix_) {
                    ierr = MatAssemblyBegin(*p_Ag_, MAT_FLUSH_ASSEMBLY);
                    CHKERRV(ierr);
                    ierr = MatAssemblyEnd(*p_Ag_, MAT_FLUSH_ASSEMBLY);
                    CHKERRV(ierr);
                }
                ierr = VecAssemblyBegin(*p_bg_);
                CHKERRV(ierr);
                ierr = VecAssemblyEnd(*p_bg_);
                CHKERRV(ierr);
            } catch (TALYException &e) {
                // if an exception was thrown (by the error handler or otherwise),
                // we need to pop the error handler on the way up the stack
                PetscPopErrorHandler();
                throw e;
            }

            PetscPopErrorHandler();
        }

        /**
         * Initialize the periodic bounds data for the equation
         *
         * This is responsible for setting up the periodic boundary conditions
         * for the solver. It also sets up the periodic exchange mapping if requested.
         *
         * @param periodic_data the PeriodicData object to use when initializing
         */
        void InitializePeriodicData(PeriodicData *periodic_data) {
            is_periodic_ = periodic_data->is_periodic();
            if (!is_periodic_) { return; }

            periodic_data_ = periodic_data;

            // if requested, set up periodic exchange object
            if (periodic_data->enable_exchange()) {
                exchanger_ = new PeriodicExchanger(periodic_data_->periodic_bounds(),
                                                   &xg_);
                exchanger_->Initialize();  // initialize object
            }

            // copy periodic variable information
            int n_periodic_vars = static_cast<int>(
                    periodic_data->periodic_var_list().size());
            periodic_vars_.fill_from_array(&periodic_data->periodic_var_list().front(),
                                           n_periodic_vars);

            FillPeriodicBoundaries();

            if (p_grid_->parallel_type_ == kWithDomainDecomp) {
                // set the correct solution partners and values for the domain
                // decomposition 'from_' vector
                const PeriodicSolutionMap *sol_partners =
                        periodic_data_->periodic_bounds()->pbc_sol_partners();
                PeriodicSolutionMap::const_iterator iter;
                for (iter = sol_partners->begin(); iter != sol_partners->end(); ++iter) {
                    LocalNodeID src_node = iter->first;
                    SolutionNodeID dst_node = iter->second;
                    for (int i_var = 0; i_var < this->periodic_vars_.size(); i_var++) {
                        int var = periodic_vars_(i_var);
                        LocalVarIdx from_val = (src_node) * n_dof() + var;
                        GlobalVarIdx to_val = (dst_node) * n_dof() + var;
                        idx_from_(from_val) = to_val;
                    }
                }
                // update the from_ index set with data for periodic bounds
                ISDestroy(&from_);
                ISCreateGeneral(PETSC_COMM_SELF, n_total_dof_, idx_from_.data(),
                                PETSC_COPY_VALUES, &from_);
            }
        }

        /**
         * Sets up the periodic boundaries based on the periodic system data
         */
        void FillPeriodicBoundaries() {
            // set the variables as periodic
            PeriodicPhysicalMap::const_iterator iter;
            const PeriodicPhysicalMap *partners = periodic_data_->periodic_bounds()->
                    pbc_partners();
            for (iter = partners->begin(); iter != partners->end(); ++iter) {
                LocalNodeID nodeID = iter->first;
                for (int k = 0; k < periodic_vars_.size(); k++) {
                    specifyPerValue(nodeID, periodic_vars_(k));
                }
            }
        }

        /**
         * Syncs values across periodic boundaries
         *
         * This is simply a wrapper around the DoPBCExchange method from the
         * PeriodicExchanger class. The send buffer needs to be filled prior to
         * calling this function and after calling this function, the resulting
         * data needs to be taken from the receive buffer.
         *
         * See comment for PeriodicExchanger class for details on how the exchange
         * should be set up and see DoPBCExchangeByID() for an example of
         * implementation.
         */
        void DoPBCExchange() {
            if (exchanger_ == NULL) {
                throw TALYException() << "Attempting to exchange data in CEquation "
                                      << "with uninitialized exchange object.";
            }
            exchanger_->DoPBCExchange();
        }

        /**
         * Syncs values across periodic boundaries
         *
         * This ensures that periodic node partners have the same data.
         * See comment for PeriodicExchanger class for more details.
         *
         * This function is used when the data to be exchanged is stored as one
         * of the items in the node data structure. To exchange values that are
         * not stored in the node data structure, use DoPBCExchange() instead
         *
         * @param id_idx index of data item to exchange
         */
        void DoPBCExchangeByID(int id_idx) {
            if (exchanger_ == NULL) {
                throw TALYException() << "Attempting to exchange data in CEquation "
                                      << "with uninitialized exchange object.";
            }

            // store data from gridfield in the send buffer
            PetscScalar *send_buffer = exchanger_->send_buffer();
            const LocalNodeID *send_indices = exchanger_->send_indices();
            // Store some data to the share array
            for (int i = 0; i < exchanger_->n_send_values(); i++) {
                LocalNodeID node_id = send_indices[i];
                send_buffer[i] = p_data_->GetNodeData(node_id).value(id_idx);
            }

            // do the actual exchange
            exchanger_->DoPBCExchange();

            // store exchanged data in the gridfield
            PetscScalar *recv_buffer = exchanger_->recv_buffer();
            const LocalNodeID *recv_indices = exchanger_->recv_indices();
            for (int i = 0; i < exchanger_->n_recv_values(); i++) {
                LocalNodeID node_id = recv_indices[i];
                p_data_->GetNodeData(node_id).value(id_idx) = recv_buffer[i];
            }
        }

        /**
         * Applies boundary conditions to the system
         *
         * This needs to be called after assembly
         * Note: Currently this is the same as this->ApplyEssBC()
         */
        virtual void ApplyAllBC() {
            this->ApplyEssBC();
        }

        /**
         * Copies solution values from nodes to their periodic partners.
         *
         * This needs to be called after 'xg_' is copied to 'solution_'
         * and before 'solution_' is copied to gridfield.
         * After this call, all periodic boundary nodes should have the correct
         * values.
         *
         * This should not be called when using domain decomposition.
         *
         * @throw TALYException if called with domain decomposition
         */
        virtual void CopyPeriodicSolution() {
            // for domain decomposition, this action is already handled by the
            // from/to copying, so this function should not be called.
            if (p_grid_->parallel_type_ == kWithDomainDecomp) {
                throw TALYException() <<
                                      "Incorrect call to CopyPeriodicSolution with domain decomposition.";
            }

            PeriodicPhysicalMap::const_iterator iter;
            const PeriodicPhysicalMap *partners = periodic_data_->periodic_bounds()->
                    pbc_partners();
            for (iter = partners->begin(); iter != partners->end(); ++iter) {
                LocalNodeID nodeID = iter->first;
                PhysicalNodeID newID = iter->second;
                for (int k = 0; k < this->periodic_vars_.size(); k++) {
                    LocalVarIdx old_loc = this->n_dof() * nodeID + this->periodic_vars_(k);
                    GlobalVarIdx new_loc = this->n_dof() * newID + this->periodic_vars_(k);
                    this->solution_(old_loc) = this->solution_(new_loc);
                }
            }
        }

        /**
         * Returns the number of rows in the Ae matrix for the given element
         *
         * This is used in the assembly process to create the Ae matrix. This would
         * need to be overridden if the Ae matrix is not in the standard shape of NxN
         * where N = n_basis_functions * n_degrees_of_freedom
         *
         * @param fe the finite element whose row count will be returned
         * @return number of rows in the Ae matrix
         */
        virtual int n_rows(const FEMElm &fe) const {
            return fe.nbf() * n_dof_;  // nbf already contains hermite scale
        }

        /**
         * Returns the number of columns in the Ae matrix for the given element
         *
         * This is used in the assembly process to create the Ae matrix. This would
         * need to be overridden if the Ae matrix is not in the standard shape of NxN
         * where N = n__basis_functions * n_degrees_of_freedom
         *
         * @param fe the finite element whose column count will be returned
         * @return number of columns in the Ae matrix
         */
        virtual int n_cols(const FEMElm &fe) const {
            return fe.nbf() * n_dof_;  // nbf already contains hermite scale
        }

        /**
         * Assembles the system from the finite elements.
         *
         * This iterates over all the elements in the system and performs the matrix
         * and vector assembly. For each element, the Integrands function is called
         * for every Gauss point in the element.
         * After assembling all the Gauss points, the surfaces for the system are
         * assembled which includes natural boundary conditions.
         * Finally, the resulting values are sent to the PETSc data structures.
         *
         * @param assemble_surface whether to do the surface assembly. Setting this
         *                         to false when there are no surface conditions will
         *                         result in a slight performance improvement.
         * @param basis_flags What values to have the basis function calculate.
         */
        virtual void AssembleVolume(bool assemble_surface = true,
                                    unsigned int basis_flags = BASIS_DEFAULT) {
            FEMElm fe(p_grid_, basis_flags);
            for (int elmID = 0; elmID < p_grid_->n_elements(); elmID++) {
                if (!IsMyElement(elmID)) continue;

                if (is_hermite()) {
                    fe.refill(elmID, BASIS_HERMITE, rel_order_);
                } else {
                    fe.refill(elmID, rel_order_);
                }

                ZeroMatrix<double> Ae;
                if (recalc_matrix_) {
                    Ae.redim(n_rows(fe), n_cols(fe));
                    Ae.fill(0);
                }
                ZEROARRAY<double> be;
                be.redim(n_rows(fe));
                be.fill(0);
                AssembleElement(elmID, Ae, be, fe, assemble_surface);
            }
        }

        /**
         * Assembles the system from the finite elements assuming unifrom elements.
         *
         * See AssembleVolume for details.
         * This version performs some simple optimizations under the assumption that
         * all elements are the same size and shape.
         *
         * @param assemble_surface whether to do the surface assembly. Setting this
         *                         to false when there are no surface conditions will
         *                         result in a slight performance improvement.
         * @param basis_flags What values to have the basis function calculate.
         */
        virtual void AssembleVolumeUniformMesh(bool assemble_surface = true,
                                               unsigned int basis_flags = BASIS_DEFAULT) {
            FEMElm fe(p_grid_, basis_flags);

            if (is_hermite()) {
                fe.refill(0, BASIS_HERMITE, rel_order_);
            } else {
                fe.refill(0, rel_order_);  // init with element 0
            }

            ZeroMatrix<double> Ae;
            if (recalc_matrix_) { Ae.redim(n_rows(fe), n_cols(fe)); }

            ZEROARRAY<double> be;
            be.redim(n_rows(fe));

            for (int elmID = 0; elmID < p_grid_->n_elements(); elmID++) {
                if (!(IsMyElement(elmID))) continue;
                fe.set_elem_hack(elmID);
                if (recalc_matrix_) { Ae.fill(0); }
                be.fill(0);
                AssembleElement(elmID, Ae, be, fe, assemble_surface);
            }
        }

        /**
         * Assembles an element.
         *
         * This function is responsible for calling the user functions to assemble
         * the element matrix and vector. The values are then put into the PETSc
         * data structures.
         *
         * @param elmID the id of the element we are assembling
         * @param Ae the element matrix to store data in
         * @param be the element vector to store data in
         * @param fe the element to assemble
         * @param assemble_surface whether to assemble the surface of the element
         */
        virtual void AssembleElement(int elmID, ZeroMatrix<double> &Ae,
                                     ZEROARRAY<double> &be, FEMElm &fe,
                                     bool assemble_surface = true) {
            switch (assembly_method_) {
                case kAssembleGaussPoints:
                    if (has_uniform_mesh_) {
                        for (std::size_t g = 0; g < feAccelerate_.size(); g++) {
                            feAccelerate_[g].set_elem_hack(fe.elem());
                            if (recalc_matrix_) {
                                Integrands(feAccelerate_[g], Ae, be);
                            } else {
                                IntegrandsVectorOnly(feAccelerate_[g], be);
                            }
                        }
                    } else {
                        while (fe.next_itg_pt()) {
                            if (recalc_matrix_) {
                                Integrands(fe, Ae, be);
                            } else {
                                IntegrandsVectorOnly(fe, be);
                            }
                        }
                    }
                    break;
                case kAssembleElements:
                    if (recalc_matrix_) {
                        IntegrandsByElement(fe, Ae, be);
                    } else {
                        IntegrandsByElementVectorOnly(fe, be);
                    }
                    break;
            }
            if (assemble_surface) {
                AssembleElementSurface(fe, elmID, Ae, be);
            }
            AssembleAebe(fe, elmID, Ae, be);
        }

        /**
         * Calculates spots for A,b in PETSc arrays.
         *
         * This determines the location for Ae and be data in the global data PETSc
         * data arrays. In the process, it applies both periodic and essential
         * boundary conditions to manipulate the locations.
         *
         * @param fe the element we are storing data from
         * @param rows_out mapping of Ae/be row to PETSc matrix row (output)
         * @param cols_out mapping of Ae column to PETSc matrix column (output)
         * @param data_ptr TODO: this is used in overloaded functins, why???
         */
        virtual void CalcAebeIndices(const FEMElm &fe, ZEROARRAY<PetscInt> &rows_out,
                                     ZEROARRAY<PetscInt> &cols_out,
                                     void *data_ptr = NULL) {
            const int vertexn = n_dof() * fe.nbf();  // should match Ae.nx() and Ae.ny()
            rows_out.redim(vertexn);
            cols_out.redim(vertexn);
            PetscInt *rows_ptr = rows_out.data();
            PetscInt *cols_ptr = cols_out.data();

            if (p_grid_->parallel_type_ == kNoDomainDecomp) {
                for (ElemNodeID i = 0; i < fe.elem()->n_nodes(); i++) {
                    for (int k = 0; k < n_dof(); k++) {
                        LocalVarIdx lid = fe.elem()->node_id_array(i) * n_dof() + k;
                        int idx = i * n_dof() + k;

                        if (has_ess_bc_.get(lid)) {
                            rows_ptr[idx] = -1;
                        } else {
                            rows_ptr[idx] = lid;
                        }
                        cols_ptr[idx] = lid;

                        if (has_per_bc_.get(lid)) {
                            LocalNodeID lclnodeID = fe.elem()->node_id_array(i);
                            SolutionNodeID newNode = periodic_data_->periodic_bounds()->
                                    GetPeriodicSolPartner(lclnodeID);
                            // we don't want this to override an essential condition
                            if (!has_ess_bc_.get(lid)) {
                                rows_ptr[idx] = newNode * n_dof() + k;
                            }
                            cols_ptr[idx] = newNode * n_dof() + k;
                        }
                    }
                }
            } else {
                const ELEM *elem = fe.elem();
                for (ElemNodeID i = 0; i < elem->n_nodes(); i++) {
                    for (int k = 0; k < n_dof(); k++) {
                        LocalVarIdx lid = elem->node_id_array(i) * n_dof() + k;
                        GlobalVarIdx gid = (p_grid_->solution_map(elem->node_id_array(i)))
                                           * n_dof() + k;
                        int idx = i * n_dof() + k;

                        if (has_ess_bc_.get(lid)) {
                            rows_ptr[idx] = -1;
                        } else {
                            rows_ptr[idx] = gid;
                        }
                        cols_ptr[idx] = gid;

                        if (has_per_bc_.get(lid)) {
                            LocalNodeID lclnodeID = elem->node_id_array(i);
                            SolutionNodeID newNode = periodic_data_->periodic_bounds()->
                                    GetPeriodicSolPartner(lclnodeID);
                            // we don't want this to override an essential condition
                            if (!has_ess_bc_.get(lid)) {
                                rows_ptr[idx] = newNode * n_dof() + k;
                            }
                            cols_ptr[idx] = newNode * n_dof() + k;
                        }
                    }
                }
            }
        }

        /**
         * Calculates spots for A,b in PETSc arrays and sends data to PETSc structures.
         *
         * Once the locations are determined, they and the Ae,be data structures
         * are sent to AssembleAebeWithIndex where they are stored in the PETSc
         * structures.
         *
         * @param fe the element we are storing data from
         * @param elmID the id of the element we are storing data from
         * @param Ae the element matrix data to send to PETSc structures
         * @param be the element vector data to send to PETSc structures
         */
        virtual void AssembleAebe(const FEMElm &fe, int elmID,
                                  const ZeroMatrix<double> &Ae,
                                  const ZEROARRAY<double> &be) {
            ZEROARRAY<PetscInt> vertex_arr, vertex_col_arr;
            CalcAebeIndices(fe, vertex_arr, vertex_col_arr);

            AssembleAebeWithIndex(Ae.data(), be.data(), vertex_arr.size(),
                                  vertex_arr.data(), vertex_col_arr.size(),
                                  vertex_col_arr.data());
        }

        /**
         * Passes data in the Ae and be structures to the PETSc data structures
         *
         * @param Ae_ptr element matrix data to send to PETSc array
         * @param be_ptr element vector data to send to PETSc array
         * @param rowNumber number of rows of data to store
         * @param rowIndex global indices of rows for data storage
         * @param colNumber number of columns of data to store
         * @param colIndex global indices of columns for data storage
         * @return 0 on success or PETSc error code if there are any errors
         *
         * NOTE: This function assumes the PetscExceptionErrorHandler is in use.
         *       As such, it uses the CHKERRV macro to check for errors, which
         *       calls the error handler but DOES NOT return on error.
         *       We rely on PetscExceptionErrorHandler to throw an exception which
         *       will unwind the stack.
         *
         * TODO: does VecSetOption need to be called every time??
         */
        void AssembleAebeWithIndex(const double *Ae_ptr,
                                   const double *be_ptr,
                                   PetscInt rowNumber, const PetscInt rowIndex[],
                                   PetscInt colNumber, const PetscInt colIndex[]) {
            // UpdateMatPreallocation() should already be called from Assemble()

            PetscErrorCode ierr;
            if (recalc_matrix_) {
                ierr = MatSetValues(*p_Ag_, rowNumber, rowIndex, colNumber, colIndex,
                                    Ae_ptr, ADD_VALUES);
                CHKERRV(ierr);
            }
            ierr = VecSetOption(*p_bg_, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
            CHKERRV(ierr);
            ierr = VecSetValues(*p_bg_, rowNumber, rowIndex, be_ptr, ADD_VALUES);
            CHKERRV(ierr);
        }

        /**
         * TODO: is this actually used???
         */
        virtual void AssembleSurface() {}

        /**
         * Assembes the element surface.
         *
         * This function loops over the surfaces of the given element and calls
         * Integrads4side to assemble the surfaces. The surface would need to be
         * assembled when there are boundary conditions that apply at the surface.
         *
         * @param fe the element to assemble
         * @param elmID the id of the element we are assembling
         * @param Ae the element matrix to store data in
         * @param be the element vector to store data in
         */
        virtual void AssembleElementSurface(FEMElm &fe, int elmID,
                                            ZeroMatrix<double> &Ae,
                                            ZEROARRAY<double> &be) {
            ELEM::SurfaceList_type::const_iterator it;
            for (it = fe.elem()->surface_indicator_.begin();
                 it != fe.elem()->surface_indicator_.end(); it++) {
                fe.refill_surface(elmID, &*it, rel_order_);

                while (fe.next_itg_pt()) {
                    for (unsigned int i = 0; i < SurfaceIndicator::MAX_SURFACE_INDICATORS;
                         i++) {
                        if (it->has_indicator(i)) {
                            if (recalc_matrix_) {
                                Integrands4side(fe, i, Ae, be);
                            } else {
                                Integrands4sideVectorOnly(fe, i, be);
                            }
                        }
                    }
                }
            }
        }

        /**
         * Fills the Ae and be structures with data from a Gauss point.
         *
         * This is the main function that needs to be implemented for a derived
         * class in order to apply the finite element method to an equation. The
         * purpose is to compute the contributions to the element matrix
         * and element vector for a given element number and integration point.
         *
         * @param fe the element we are assembling a Gauss point from
         * @param Ae the element matrix to put data in
         * @param be the element vector to put data in
         */
        virtual void Integrands(const FEMElm &fe, ZeroMatrix<double> &Ae,
                                ZEROARRAY<double> &be) {}

        /**
         * Fills only the "be" structures with data from a Gauss point.
         *
         * This version of Integrands only needs to be implimented by derived classes
         * that need an assembly that only assembles the b vector.
         *
         * @param fe the element we are assembling a Gauss point from
         * @param be the element vector to put data in
         */
        virtual void IntegrandsVectorOnly(const FEMElm &fe, ZEROARRAY<double> &be) {}

        /**
         * Fills the Ae and be structures with data from an element.
         *
         * This is element based versino of Integrands. the main function that needs
         * to be implemented for a derived class in order to use element based
         * assembly. The purpose is to compute the contributions to the element matrix
         * and element vector for a given element.
         *
         * @param fe the element we are assembling
         * @param Ae the element matrix to put data in
         * @param be the element vector to put data in
         */
        virtual void IntegrandsByElement(const FEMElm &fe, ZeroMatrix<double> &Ae,
                                         ZEROARRAY<double> &be) {}

        /**
         * Fills only the "be" structures with data from an element
         *
         * This version of Integrands only needs to be implimented by derived classes
         * that need an assembly that only assembles the b vector.
         *
         * @param fe the element we are assembling a Gauss point from
         * @param be the element vector to put data in
         */
        virtual void IntegrandsByElementVectorOnly(const FEMElm &fe,
                                                   ZEROARRAY<double> &be) {}

        /**
         * Fills the Ae and be structures with data from a Gauss point on a surface.
         *
         * Similar to "Integrands", this function can be implemented in the derived
         * class in order to contributions to the element matrix and element vector
         * for a Gauss point on a surface. These Gauss points might be present in the
         * weak formulation of the PDE problem to represent natural boundary
         * conditions. If the system does not have these boundary conditions, there
         * is no need implement this function. The default value will do nothing.
         *
         * @param fe the element we are assembling a Gauss point from
         * @param sideInd index of the side to assemble
         * @param Ae the element matrix to put data in
         * @param be the element vector to put data in
         */
        virtual void Integrands4side(const FEMElm &fe, int sideInd,
                                     ZeroMatrix<double> &Ae,
                                     ZEROARRAY<double> &be) {}

        /**
         * Fills be structures with data from a Gauss point on a surface.
         *
         * Similar to "Integrands", this function can be implemented in the derived
         * class in order to contributions to the element vector
         * for a Gauss point on a surface. These Gauss points might be present in the
         * weak formulation of the PDE problem to represent natural boundary
         * conditions. If the system does not have these boundary conditions, there
         * is no need implement this function. The default value will do nothing.
         *
         * @param [in] fe the element we are assembling a Gauss point from
         * @param [in] side_idx boundary type of the side
         * @param [in] id id of the boundary type
         * @param [out] Ae the element matrix to put data in
        */
        virtual void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, const unsigned int side_idx, const unsigned int id, TALYFEMLIB::ZEROARRAY<double> &be) {

            throw TALYFEMLIB::TALYException() << "[ERROR] : Something is wring in your implementation \n You need to implement the Integrands4Side_be" ;
        }

        /**
         * Fills Ae structures with data from a Gauss point on a surface.
         *
         * Similar to "Integrands", this function can be implemented in the derived
         * class in order to contributions to the element matrix
         * for a Gauss point on a surface. These Gauss points might be present in the
         * weak formulation of the PDE problem to represent natural boundary
         * conditions. If the system does not have these boundary conditions, there
         * is no need implement this function. The default value will do nothing.
         *
         * @param [in] fe the element we are assembling a Gauss point from
         * @param [in] side_idx boundary type of the side
         * @param [in] id id of the boundary type
         * @param [out] Ae the element matrix to put data in
        */
        virtual void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, const unsigned int side_idx, const unsigned int id, TALYFEMLIB::ZeroMatrix<double> &Ae) {
            throw TALYFEMLIB::TALYException() << "[ERROR] : Something is wring in your implementation \n You need to implement the Integrands4Side_Ae" ;
        }




        /**
         * Fills only the "be" structures with data from a Gauss point on a surface.
         * This is called instead of Integrands4side(...) when recalc_matrix_ is
         * set to false.
         *
         * @param fe the element we are assembling
         * @param sideInd index of the side to assemble
         * @param be the element vector to put data in
         */
        virtual void Integrands4sideVectorOnly(const FEMElm &fe, int sideInd,
                                               ZEROARRAY<double> &be) {}


        /**
         * Used to calculate nonzeros for the "perfect" preallocator.
         * Ae is initially filled with false.
         * @param fe the element we are assembling
         * @param Ae set each index that will ever be nonzero in Integrands() to true
         */
        virtual void IntegrandsPreallocator(const FEMElm &fe, ZeroMatrix<bool> &Ae) {
            Ae.fill(true);
        }

        /**
         * Applies essential boundary conditions to the matrix system.
         *
         * This function changes the global A matrix and b vector to apply the
         * effects of essential boundary conditions.
         *
         * The essential condtions from periodic bounds are applied prior to the
         * conditions from Dirichlet boundaries. This ensures that a Dirichlet
         * condition will be made periodic if places on a periodic boundary.
         *
         * @return 0 on success and a PETSc error code if there are any errors
         */
        virtual PetscErrorCode ApplyEssBC() {
            PetscErrorCode ierr;

            boundary_conditions_.ApplyPeriodicConditions(*p_Ag_, *p_bg_,
                                                         p_x_map_, recalc_matrix_);
            boundary_conditions_.ApplyDirichletConditions(*p_Ag_, *p_bg_,
                                                          p_x_map_, recalc_matrix_);

            if (recalc_matrix_) {
                ierr = MatAssemblyBegin(*p_Ag_, MAT_FINAL_ASSEMBLY);
                CHKERRQ(ierr);
                ierr = MatAssemblyEnd(*p_Ag_, MAT_FINAL_ASSEMBLY);
                CHKERRQ(ierr);
            }

            ierr = VecAssemblyBegin(*p_bg_);
            CHKERRQ(ierr);
            ierr = VecAssemblyEnd(*p_bg_);
            CHKERRQ(ierr);

            return 0;
        }

        /**
         * Applies the essential boundary conditions to the solution vector
         *
         * This function loops over the boundary conditions and applies them to the
         * solution vector. Applying the boundary conditions consists of adjusting
         * the solution vector by changing entries to match the values given by the
         * essential boundary conditions.
         *
         * @param pSolution the solution for applying essential boundary condition.
         *
         * TODO: will this ever be applied to anything other than solution??
         */
        virtual void ApplyEssBCToSolution(ZEROARRAY<double> *pSolution = NULL) {
            if (pSolution == NULL) {
                boundary_conditions_.ApplyPeriodicConditionsToSolution(solution_);
                boundary_conditions_.ApplyDirichletConditionsToSolution(solution_);
            } else {
                boundary_conditions_.ApplyPeriodicConditionsToSolution(*pSolution);
                boundary_conditions_.ApplyDirichletConditionsToSolution(*pSolution);
            }
        }

        /**
         * Sets pointers to the solver data structures
         *
         * @param pAg pointer to the global A matrix
         * @param pbg pointer to the global b vector
         * @param pXMap pointer to map between local and global unknown numbers
         */
        void setGlobalAb(Mat *pAg, Vec *pbg, ZEROARRAY<GlobalVarIdx> *pXMap) {
            this->p_Ag_ = pAg;
            this->p_bg_ = pbg;
            this->p_x_map_ = pXMap;
        }

        /**
         * Returns 1(true) if the element belongs to this process
         *
         * This is needed for the case of a parallel system without domain
         * decomposition. In that case, all processes have access to all elements
         * regardless of whether the process actually "owns" the element. This
         * function is used to determine whether an element is owned by the current
         * process.
         *
         * @param ElmID the element ID to determine ownership of
         * @return 1 if the element belongs to this process, 0 otherwise
         *         For a serial process, this always returns 1
         */
        virtual int IsMyElement(int ElmID) const {
            if (is_serial_) return true;
            return p_grid_->IsMyElement(ElmID);
        }

        /**
         * Assemble and solve the system
         *
         * This needs to be implemented for a derived class.
         * In the typical case, this function will assemble the matrix, apply boundary
         * conditions, call the PETSc solver, and transfer the resulting data from
         * the PETSc arrays to the library data structures.
         *
         * @param dt time interval for the solver to advance for timestepping systems
         * @param t current computation time for timestepping systems
         */
        virtual void Solve(double dt, double t) = 0;

        /**
         * Returns the norm(NORM_2) of the solution(xg_)
         *
         * @return The L2 norm of the solution
         */
        double CalculateSolutionNorm() const {
            double ret;
            VecNorm(xg_, NORM_2, &ret);
            return ret;
        }

        /**
         * Sets whether or not to recalculate the A matrix during assembly
         *
         * @param recalc_matrix new value for matrix_recalc_
         */
        inline void set_recalc_matrix(bool recalc_matrix) {
            recalc_matrix_ = recalc_matrix;
        }

        /**
         * Returns true if the system is serial
         *
         * @return true if the system is serial
         */
        inline bool is_serial() const {
            return is_serial_;
        }

        /**
         * Returns the assembly method used by the solver.
         *
         * @return the assembly method (kAssembleElements or kAssembleGaussPoints)
         */
        inline AssemblyMethod assembly_method() const {
            return assembly_method_;
        }

        /**
         * Returns the current simulation time
         *
         * @return the current simulation time
         */
        inline double t() const {
            return t_;
        }

        /**
         * Sets the current simulation time
         *
         * @param new_t the new value of the simulation time
         */
        inline void set_t(double new_t) {
            t_ = new_t;
        }

        /**
         * Returns the value of the timestep
         *
         * @return the value of the timestep
         */
        inline double dt() const {
            return dt_;
        }

        /**
         * Sets the the value of the timestep
         *
         * @param new_dt the new value of the timestep
         */
        inline void set_dt(double new_dt) {
            dt_ = new_dt;
        }

        /**
         * Returns reference to "bg_" array
         *
         * @return reference to "bg_" array
         */
        const Vec &bg() const {
            return bg_;
        }

        /**
         * Returns reference to "Ag_" mat
         *
         * @return reference to "Ag_" mat
         */
        const Mat &Ag() const {
            return Ag_;
        }

        /**
         * Returns reference to "from_" array
         *
         * @return reference to "from_" array
         */
        const IS &from() const {
            return from_;
        }

        /**
         * Returns reference to "to_" array
         *
         * @return reference to "to_" array
         */
        const IS &to() const {
            return to_;
        }

        /**
         * Returns total number of degrees of freedom
         *
         * @return the total number of degrees of freedom
         */
        int n_total_dof() const {
            return n_total_dof_;
        }

        /**
         * Returns reference to "idx_from_" array
         *
         * @return reference to "idx_from_" array
         */
        const ZEROARRAY<PetscInt> &idx_from() const {
            return idx_from_;
        }

        /**
         * Returns reference to "idx_to_" array
         *
         * @return reference to "idx_to_" array
         */
        const ZEROARRAY<PetscInt> &idx_to() const {
            return idx_to_;
        }

        /**
         * Returns order of integration of the equation
         */
        int rel_order_of_itg() const {
            return rel_order_;
        }

        /**
         * Returns number of degrees of freedom in equation
         */
        int n_dof() const {
            return n_dof_ * hermite_n_dof_scale_;
        }

        /**
         * @returns if using hermite basis functions (ndof multiplier)
         */
        bool is_hermite() const {
            return (hermite_n_dof_scale_ != 1);
        }

        /**
         * @returns flags automatically passed to FEMElm
         */
        inline unsigned int basis_flags() const {
            return basis_flags_;
        }

        /**
         * @param flags set of basis flags to use
         */
        inline void set_basis_flags(unsigned int flags) {
            basis_flags_ = flags;
        }

        /**
         * @param flag basis flag(s) to add
         */
        inline void add_basis_flag(BasisFlags flag) {
            basis_flags_ |= flag;
        }

        /**
         * @param flag basis flag(s) to remove
         */
        inline void remove_basis_flag(BasisFlags flag) {
            basis_flags_ = basis_flags_ & (~flag);
        }

        GRID *p_grid_;  ///< pointer to grid data
        GridField<NodeData> *p_data_;  ///< pointer to grid field data

        ZEROARRAY<double> solution_;  ///< solution
        std::vector<FEMElm> feAccelerate_;  ///< array of finite elements used when
        ///< acceleration is turned on
        ZEROARRAY<ZeroMatrix<ZeroMatrix<double> > > dNdNArray_;
        ///< array of dN*dN values used when using uniform mesh acceleration.
        ///< This is an array of length n_gauss_points with each entry being a matrix
        ///< of matrices that store the produce of dN values at each gauss point.
        ///< dNdNArray_[g](a, b)(i, j) = dn(a, i) * dn(b, j)
        ZEROARRAY<ZeroMatrix<ZEROARRAY<double> > > NdNArray_;
        ///< array of cached N*dN values when using uniform mesh acceleration
        ZEROARRAY<ZeroMatrix<double> > sumdNdNArray_;
        ///< sum of dN*dN values for each gauss point using uniform mesh acceleration.
        ///< This is an array of length n_gauss_points with each entry being a matrix
        ///< with terms equal to the sum ofver all i,j of dn(a, i) * dn(b, j)
        ZEROARRAY<ZeroMatrix<double> > NNArray_;
        ///< array of N*N values used when using uniform mesh acceleration.
        ///< This is an array of length n_gauss_points with each entry being a matrix
        ///< of N*N values used in assembly.
        ///< NNArray_[g](a, b) = N(a) * N(b)

        // data for elemental assembly
        ZeroMatrix<double> ea_NN_;  ///< weighted sum of fe.N*fe.N for use during
        ///< elemental assembly
        ZeroMatrix<double> ea_dNdN_;  ///< weighted sum of fe.dN*fe.dN for use during
        ///< elemental assembly
        ZeroMatrix<ZEROARRAY<double> > ea_NNdN_gauss_sum_;  ///< N*N*dN for assembly

    protected:
        /**
         * Returns number of periodic variables (used in test code)
         */
        int NumPeriodicVars() const { return periodic_vars_.size(); }

        /**
         * Returns the ith periodic variable (used in test code)
         *
         * @param index index of periodic value desired
         * @return ith periodic variable
         */
        int IthPeriodicVar(int index) const { return periodic_vars_(index); }

        double dt_;  ///< the size of time interval
        double t_;  ///< the current time

        int n_total_dof_;  ///< The total number of degree freedom.

        Mat Ag_;  ///< stiffness matrix
        Vec bg_;  ///< load
        Vec xg_;  ///< solution in linear system

        int n_dof_;  ///< The number of degrees of freedom per node.
        int hermite_n_dof_scale_;  ///< ndof scaling factor when using Hermite basis

        ZEROARRAY<bool> has_ess_bc_;  ///< An array for quickly judge whether there is
        ///< an ess bc. (Mainly for fast applying essbc)
        ZEROARRAY<bool> has_per_bc_;  ///< An array to quickly judge whether there is
        ///< a periodic bc. (Mainly for applying perbc)
        ZEROARRAY<bool> is_node_periodic_;  ///< An array specifying whether a given
        ///< local node is periodic.

        bool has_uniform_mesh_;  ///< whether the underlying mesh is uniform

        AssemblyMethod assembly_method_;  ///< method of assemblying the matrix
        PeriodicExchanger *exchanger_;  ///< exchange communication object
        PeriodicData *periodic_data_;  ///< periodic data object

        IS from_;  ///< index set with indices for PETSc scatter of values when using
        ///< domain docomposition. These are the indices where the data
        ///< is sent from.
        IS to_;  ///< index set with indices for PETSc scatter of values when using
        ///< domain docomposition. These are the indices where the data
        ///< is sent to.

        int created_;  ///< Whether memory for linear system has been allocated

        Mat *p_Ag_;  ///< pointer to global Ag matrix
        Vec *p_bg_;  ///< pointer to global Bg vector

        KSP ksp_;  ///< KSP system

        BoundaryConditions boundary_conditions_;  ///< boundary condition data

        ZEROARRAY<GlobalVarIdx> *p_x_map_;  ///< &x_map_
        ZEROARRAY<GlobalVarIdx> x_map_;  ///< local matrix to global matrix map

        ZEROARRAY<PetscInt> idx_from_;  ///< indices for exchange of values when using
        ///< domain docomposition. These are the
        ///< indices the data is sent from.
        ZEROARRAY<PetscInt> idx_to_;  ///< indices for exchange of values when using
        ///< domain docomposition. These are the indices
        ///< the data is sent to.

        bool is_serial_;  ///< whether this is a serial calculation

        bool is_periodic_;  ///< whether the system has any periodic boundaries

        ZEROARRAY<int> periodic_vars_;  ///< indicating which variables are periodic

        int rel_order_;  ///< relative order of integration
        unsigned int basis_flags_;  ///< flags to pass into FEMElm constructor

        Preallocator *preallocator_;  ///< preallocator to use for the matrix
        bool matPreallocationDone_;  ///< whether the preallocation has been done
        bool preallocator_is_owned_;  ///< if we own (need to delete) preallocator

        bool recalc_matrix_;  ///< whether to recalculate the 'A' matrix
    };

}  // namespace TALYFEMLIB

#endif  // FEM_CEQUATION_H_
