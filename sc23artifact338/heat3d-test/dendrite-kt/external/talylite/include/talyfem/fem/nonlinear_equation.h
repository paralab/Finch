#pragma once

#include <talyfem/fem/cequation.h>

namespace TALYFEMLIB {

enum SurfaceIntegrationBehavior {
  SKIP_SURFACE_INTEGRATION,  ///< disable surface assembly (Integrands4side) (faster if not needed)
  ENABLE_SURFACE_INTEGRATION,  ///< enable surface assembly (Integrands4side) (default)
};

/**
 * A base class for nonlinear equations that automatically manages a PETSc SNES object.
 * This class provides an implementation for CEquation::Solve().
 * Users should fill in Integrands(), Integrands4side(), and fillEssBC(), as you would for a typical CEquation class.
 * In Integrands, you should calculate the residual in be and the Jacobian in Ae.
 *
 * IMPORTANT NOTE: Essential boundary conditions work differently for nonlinear equations. specifyValue() will set the
 * corresponding row in the Jacobian to be 1 and the value in the residual to the given value.
 * In other words, specifyValue is actually specifying the value of the residual, *not* the solution!
 * Calling specifyValue() for any value other than 0 is usually a mistake.
 *
 * This class makes use of TalyFEM's PetscEventLogger to time different parts of the code. Run your program with
 * -log_view to see a summary of event times. See the PETSc documentation for different output formats which may
 * contain more information (for example, there's a Python file output which is pretty useful).
 */
template <typename NodeData>
class NonlinearEquation : public CEquation<NodeData> {
 public:
  /// a = dof_map_[i].first, b = dof_map_[i].second means degree of freedom a maps to NodeData index b
  typedef std::vector< std::pair<int, int> > DofMap;

  enum DivergedBehavior {
    ON_DIVERGED_IGNORE,  ///< don't do anything if SNESSolve diverges
    ON_DIVERGED_PRINT_WARNING,  ///< print a warning using PrintWarning if SNESSolve diverges
    ON_DIVERGED_ABORT,  ///< call MPI_Abort if SNESSolve diverges
  };

  explicit NonlinearEquation(SurfaceIntegrationBehavior surface_behavior = ENABLE_SURFACE_INTEGRATION,
      bool uniform_mesh = false, AssemblyMethod asm_method = kAssembleGaussPoints)
    : CEquation<NodeData>(uniform_mesh, asm_method), snes_(NULL), surface_behavior_(surface_behavior),
      diverged_behavior_(ON_DIVERGED_PRINT_WARNING) {
  }

  virtual ~NonlinearEquation() {
    cleanupNonlinear();
  }


  // overridden to automatically initialize our SNES object
  using CEquation<NodeData>::redimSolver;
  virtual void redimSolver(GRID* p_grid, int n_dof_value,
                           bool hermite, int rel_order,
                           PeriodicData *periodic_data = NULL) override {
    // call base class version of redimSolver first
    CEquation<NodeData>::redimSolver(p_grid, n_dof_value, hermite, rel_order, periodic_data);

    this->UpdateMatPreallocation();
    initNonlinear();
  }

  virtual void Solve(double dt, double t) override {
    assert (this->snes_ != NULL);
    assert (this->xg_ != NULL);

    // set t_ and dt_ so we can use it in Integrands during assembly
    this->set_t(t);
    this->set_dt(dt);

    // this is the only time fillEssBC() is called, since BC are constant for the whole timestep
    // note that fillEssBC may change p_data_, and thus influence our initial conditions
    this->fillEssBC();

    // copy current state of the GridField into our xg as our initial guess
    copyGridFieldToVec(this->p_data_, this->xg_, dof_map_);

    // do nonlinear solve
    PetscErrorCode ierr = SNESSolve(snes_, NULL, this->xg_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    // copy final solution into the gridfield
    copyVecToGridField(this->xg_, this->p_data_, dof_map_);

    // if we diverged, print a warning
    SNESConvergedReason converged_reason;
    SNESGetConvergedReason(snes_, &converged_reason);
    if (converged_reason < 0) {
      // diverged
      if (diverged_behavior_ == ON_DIVERGED_PRINT_WARNING) {
        PrintWarning("Non-linear solve diverged.");
      } else if (diverged_behavior_ == ON_DIVERGED_ABORT) {
        MPI_Abort(PETSC_COMM_WORLD, converged_reason);
      }
    }
  }

  /**
   * Returns NULL before redimSolver() has been called.
   * @return a reference to the SNES object (useful for manually setting tolerances, etc...)
   */
  inline SNES snes() {
    return snes_;
  }

  /**
   * Change behavior when SNES diverges. Default behavior is to print a warning (ON_DIVERGED_PRINT_WARNING).
   * Set to ON_DIVERGED_IGNORE to do nothing.
   * @param behavior new behavior
   */
  void setDivergedBehavior(DivergedBehavior behavior) {
    diverged_behavior_ = behavior;
  }

  /**
   * Change whether or not to do surface integration during assembly.
   * The initial value is set by the constructor, but you can use this function to change it later.
   * @param behavior new behavior
   */
  void setSurfaceIntegration(SurfaceIntegrationBehavior behavior) {
    surface_behavior_ = behavior;
  }

 protected:
  // NOTE: The usage of "this" in this class is due to weirdness when accessing the base class from a templated class
  // that inherits a templated class (yo dawg...)

  /**
   * Register a degree of freedom <-> NodeData value mapping.
   * This is required for degrees of freedom be automatically copied to and from the GridField.
   * @param dof_idx degree of freedom index (in the range [0..n_dof])
   * @param node_data_idx NodeData::value(i) index (in the range [0..NodeData::valueno()])
   */
  void addDof(int dof_idx, int node_data_idx) {
    dof_map_.push_back({dof_idx, node_data_idx});
  }

  /**
   * This function is called by FormJacobianCallback to evaluate the Jacobian based on some "guess" (potential solution).
   * @param[out] PETSc matrix to fill with the jacobian
   * @param guess "guess solution" to evaluate the jacobian from
   */
  virtual void formJacobian(Mat jacobian, Vec guess) {
    // TODO change talyfem so we can assemble directly into jacobian...
    if (jacobian != this->Ag_)
      throw TALYException() << "Jacobian does not match Ag, unfortunate past design decisions";

    copyVecToGridField(guess, this->p_data_, dof_map_);
    // this->recalc_matrix_ = true;
    this->Assemble(surface_behavior_ == ENABLE_SURFACE_INTEGRATION);
    this->ApplyEssBC();
  }

  /**
   * This function is called by FormFunctionCallback to evaluate the residual based on some "guess" (potential solution).
   * @param[out] residual PETSc vector to fill with the residual
   * @param guess "guess solution" to evaluate the residual from
   */
  virtual void formFunction(Vec residual, Vec guess) {
    copyVecToGridField(guess, this->p_data_, dof_map_);
    // this->recalc_matrix_ = false;
    this->Assemble(surface_behavior_ == ENABLE_SURFACE_INTEGRATION);
    this->ApplyEssBC();

    // TODO change TalyFEM so we can assemble directly into bg...
    VecCopy(this->bg_, residual);
  }

  /**
   * Copies vec[node_id * vec_ndof + dof_map[i].first] into gf.GetNode(node_id).value(dof_map[i].second)
   * @param[in] vec vector to copy from
   * @param[inout] gf gridfield to copy into
   * @param dof_map Which degrees of freedom to copy where. first value is vector dof, second is node data index.
   *                May not cover all degrees of freedom in vec.
   */
  void copyVecToGridField(Vec vec, GridField<NodeData>* gf, const DofMap& dof_map) const {
    PetscEventLogger ev("copyVecToGridField");
    PetscErrorCode ierr;
    MPI_Comm comm = PETSC_COMM_WORLD;

    const PetscInt vec_ndof = this->n_dof();
    assert (vec_ndof > 0);
    const GRID* grid = gf->p_grid_;

    // tmp_data will hold the new data for the gridfield
    std::vector<PetscScalar> tmp_data(this->n_total_dof());

    // create a PETSc vector backed by tmp_data
    Vec tmp_vec;
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, vec_ndof, tmp_data.size(), tmp_data.data(), &tmp_vec);
    CHKERRABORT(comm, ierr);

    // set up a scatter context for copying vec into tmp_vec
    VecScatter scatter;
    if (grid->parallel_type_ == kWithDomainDecomp) {
      ierr = VecScatterCreate(vec, this->from(), tmp_vec, this->to(), &scatter);
      CHKERRABORT(comm, ierr);
    } else {
      // copy all of vec (which each process only has part of) into tmp_vec (which each process has all of)
      ierr = VecScatterCreateToAll(vec, &scatter, NULL);
      CHKERRABORT(comm, ierr);
    }

    // do the scatter
    ierr = VecScatterBegin(scatter, vec, tmp_vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRABORT(comm, ierr);
    ierr = VecScatterEnd(scatter, vec, tmp_vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRABORT(comm, ierr);

    // copy tmp_vec into the gridfield
    for (LocalNodeID node_id = 0; node_id < gf->p_grid_->n_nodes(); node_id++) {
      for (unsigned int i = 0; i < dof_map.size(); i++) {
        int dof_idx = dof_map[i].first;
        int node_data_idx = dof_map[i].second;
        gf->GetNodeData(node_id).value(node_data_idx) = tmp_data[node_id * vec_ndof + dof_idx];
      }
    }

    // cleanup
    ierr = VecScatterDestroy(&scatter); CHKERRABORT(comm, ierr);
    ierr = VecDestroy(&tmp_vec); CHKERRABORT(comm, ierr);
  }

  /**
   * Copies gf.GetNode(node_id).value(dof_map[i].second) into vec[node_id * vec_ndof + dof_map[i].first].
   * Opposite of copyVecToGridField.
   * @param[in] gf gridfield to copy from
   * @param[inout] vec vector to copy into
   * @param dof_map Which degrees of freedom to copy where. first value is vector dof, second is node data index.
   *                May not cover all degrees of freedom in vec.
   */
  void copyGridFieldToVec(const GridField<NodeData>* gf, Vec vec, const DofMap& dof_map) const {
    PetscEventLogger ev("copyGridFieldToVec");

    PetscErrorCode ierr;
    MPI_Comm comm = PETSC_COMM_WORLD;

    const PetscInt vec_ndof = this->n_dof();
    assert (vec_ndof > 0);
    const GRID* grid = gf->p_grid_;

    // copy GridField data into tmp_data accoridng to dof_map
    std::vector<PetscScalar> tmp_data(this->n_total_dof());
    for (LocalNodeID node_id = 0; node_id < gf->p_grid_->n_nodes(); node_id++) {
      for (unsigned int i = 0; i < dof_map.size(); i++) {
        int dof_idx = dof_map[i].first;
        int node_data_idx = dof_map[i].second;
        tmp_data[node_id * vec_ndof + dof_idx] = gf->GetNodeData(node_id).value(node_data_idx);
      }
    }

    if (grid->parallel_type_ == kWithDomainDecomp) {
      // copy tmp_data into vec by using CEquation's to()/from() index sets and a PETSc VecScatter
      // this takes care of all the nasty details of distributed arrays for us

      // make a PETSc vector (backed by our std::vector) so we can scatter
      // (note that PETSc does not take ownership of tmp_data, it is just using it)
      Vec tmp_vec;
      ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, vec_ndof, tmp_data.size(), tmp_data.data(), &tmp_vec);
      CHKERRABORT(comm, ierr);

      // set up the scatter context
      VecScatter scatter;
      ierr = VecScatterCreate(tmp_vec, this->to(), vec, this->from(), &scatter); CHKERRABORT(comm, ierr);

      // do the scatter (i.e. copy tmp_data into vec)
      ierr = VecScatterBegin(scatter, tmp_vec, vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRABORT(comm, ierr);
      ierr = VecScatterEnd(scatter, tmp_vec, vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRABORT(comm, ierr);

      // cleanup
      ierr = VecScatterDestroy(&scatter); CHKERRABORT(comm, ierr);
      ierr = VecDestroy(&tmp_vec); CHKERRABORT(comm, ierr);
    } else {
      // every process has all of the new data in tmp_data, but only owns a portion of vec
      // each process writes tmp_data to the part of vec that it owns, which avoids communication

      // get a pointer to this process's portion of vec as vec_data (should not cause any copying)
      PetscScalar* vec_data;
      VecGetArray(vec, &vec_data);

      // figure out what portion of vec vec_data is pointing to
      PetscInt start, end;
      VecGetOwnershipRange(vec, &start, &end);

      // copy tmp_data into vec_data (for our local portion only)
      for (unsigned int i = start; i < end; i++) {
        vec_data[i - start] = tmp_data[i];
      }

      // cleanup
      VecRestoreArray(vec, &vec_data);
    }
  }

 private:
  /**
   * Create and initialize snes_ and checks dof_map_.
   * This function is called automatically by NonlinearEquation::redimSolver.
   * This function may be called multiple times.
   * Can only be called after the matrix and vector (Ag/bg) have been created (i.e. after redimSolver/Create).
   */
  void initNonlinear() {
    // make sure we don't leak SNES objects if this function is called twice
    if (snes_ != NULL)
      cleanupNonlinear();

    assert (this->Ag() != NULL);
    assert (this->bg() != NULL);

    PetscErrorCode ierr;
    MPI_Comm comm = PETSC_COMM_WORLD;

    // create and set up the SNES object
    ierr = SNESCreate(comm, &snes_); CHKERRABORT(comm, ierr);
    ierr = SNESSetFunction(snes_, this->bg_, &FormFunctionCallback, this); CHKERRABORT(comm, ierr);
    ierr = SNESSetJacobian(snes_, this->Ag_, this->Ag_, &FormJacobianCallback, this); CHKERRABORT(comm, ierr);

    // apply command line arguments to the SNES object (this also updates the KSP and PC objects owned by snes_)
    ierr = SNESSetFromOptions(snes_); CHKERRABORT(comm, ierr);

    // print warnings if dof_map_ isn't set up correctly
    validateDofMap();
  }

  /**
   * Frees the snes_ object allocated by initNonlinear.
   * This function is automatically called by the NonlinearEquation destructor.
   */
  void cleanupNonlinear() {
    SNESDestroy(&snes_);
  }

  /**
   * This function is called automatically during SNESSolve by PETSc to build the Jacobian for a particular "guess."
   * It's hooked up by calling SNESSetJacobian() during initNonlinear.
   * @param snes snes object, assumed to match snes_
   * @param guess "guess solution" to evaluate jacobian for
   * @param[out] jac jacobian we need to fill in
   * @param[out] jac_for_preconditioner matrix used for constructing the preconditioner. In this class, it is always
   *                               the same as jac (see SNESSetJacobian() in initNonlinear), so we ignore it.
   * @param ctx user data set in SNESSetJacobian - in this class, it is our NonlinearEquation object
   * @return 0 on success, as per PETSc convention
   */
  static PetscErrorCode FormJacobianCallback(SNES snes, Vec guess, Mat jac, Mat jac_for_preconditioner, void* ctx) {
    // convert C-style void* userdata to NonlinearEquation type
    auto nl = static_cast<NonlinearEquation<NodeData>*>(ctx);

    // enforce assumption that this function will only ever be called by our NonlinearEquation-managed SNES object
    assert (snes == nl->snes_);

    // enforce assumption that jac_for_preconditioner is unused
    assert (jac == jac_for_preconditioner);

    // call C++-style virtual method to do the actual work
    nl->formJacobian(jac, guess);

    return 0;
  }

  /**
   *This function is called automatically during SNESSolve by PETSc to evaluate the residual for a particular "guess."
   * It's hooked up by calling SNESSetFunction() during initNonlinear.
   * @param snes snes object, assumed to match snes_
   * @param guess "guess solution" to evaluate residual for
   * @param[out] residual vector to put the residual in
   * @param ctx user data set in SNESSetFunction - in this class, it is our NonlinearEquation object
   * @return 0 on success, as per PETSc convention
   */
  static PetscErrorCode FormFunctionCallback(SNES snes, Vec guess, Vec residual, void* ctx) {
    // convert C-style void* userdata to NonlinearEquation type
    auto nl = static_cast<NonlinearEquation<NodeData>*>(ctx);

    // enforce assumption that this function will only ever be called by our NonlinearEquation-managed SNES object
    assert (snes == nl->snes_);

    // call C++-style virtual method to do the actual work
    nl->formFunction(residual, guess);

    return 0;
  }

  /**
   * Called during initialization to validate dof_map_. This isn't done inside addDof() because n_dof is not set
   * until redimSolver() has been called, but the recommended usage of addDof() is from the constructor.
   */
  void validateDofMap() {
    if (dof_map_.empty()) {
      PrintWarning("NonlinearEquation - DOF to NodeData map is empty. "
                   "Try adding NonlinearEquation::addDof() in your constructor?");
    }

    for (unsigned int i = 0; i < dof_map_.size(); i++) {
      if (dof_map_[i].first < 0 || dof_map_[i].first >= this->n_dof()) {
        PrintError((i+1), "th call to NonlinearEquation::addDof() - dof index ", i, " is out of bounds!");
      }
      if (dof_map_[i].second < 0 || dof_map_[i].second >= NodeData::valueno()) {
        PrintWarning((i+1), "th call to NonlinearEquation::addDof() - NodeData value index ", i, " is out of bounds!");
      }
    }
  }

  SNES snes_;  ///< nonlinear solver
  DofMap dof_map_;  ///< mapping of dofs to node data indices for copyVecToGridField/copyGridFieldToVec
  DivergedBehavior diverged_behavior_;  ///< what to do if SNESSolve diverges
  SurfaceIntegrationBehavior surface_behavior_;  ///< if we should assemble surfaces using Integrands4side
};

}
