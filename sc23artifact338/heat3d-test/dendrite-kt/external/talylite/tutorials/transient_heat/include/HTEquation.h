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
#ifndef INCLUDE_HTEQUATION_H_
#define INCLUDE_HTEQUATION_H_

#include "HTAnalyticSolution.h"
#include "HTNodeData.h"
#include "HTInputData.h"
#include <talyfem/fem/preallocator_perfect.h>

/**
 * This class implements a transient heat problem.
 *
 * The class can handle 1, 2, or 3 dimensions. Boundary conditions are
 * Direchlet (fixed value on the boundaries) with the boundary values set to be
 * equal to the anayltic value at the boundaries.
 *
 * The system is solved using a time stepping method with a Backward Euler
 * integration scheme.
 */
class HTEquation : public CEquation<HTNodeData> {
 public:
  // indices to timing arrays. These are locations in the timers_ array that
  // correspond to specific stages of the code that we wish to time.
  static const int kTimerSolve = 0;  ///< index to timer for solve process
  static const int kTimerAssemble = 1;  ///< index to timer for assemble
  static const int kTimerKSPSolve = 2;  ///< index to timer for KSPsolve
  static const int kTimerUpdate = 3;  ///< index to timer for update process

  /**
   * Constucts the solver by setting up timers.
   */
  explicit HTEquation(HTInputData* input_data,
                      HTAnalyticSolution *analytic_solution,
                      bool has_uniform_mesh = false,
                      AssemblyMethod assemble_method = kAssembleGaussPoints)
      : CEquation<HTNodeData>(has_uniform_mesh, assemble_method),
        input_data_(input_data),
        analytic_solution_(analytic_solution) {
    timers_[kTimerSolve].set_label("Solve");
    timers_[kTimerAssemble].set_label("Assemble");
    timers_[kTimerKSPSolve].set_label("KSPSolve");
    timers_[kTimerUpdate].set_label("Update");

    // SetPreallocator(new PreallocatorPerfect<HTNodeData>(this));
  }

  /**
   * Destroys the object
   *
   * This also prints the average times required for each of the processes that
   * we times during execution.
   */
  virtual ~HTEquation() {
    timers_[kTimerSolve].PrintGlobalAverageSeconds();
    timers_[kTimerAssemble].PrintGlobalAverageSeconds();
    timers_[kTimerKSPSolve].PrintGlobalAverageSeconds();
    timers_[kTimerUpdate].PrintGlobalAverageSeconds();
  }

  /**
   * Sets up the essential (Direchlet) boundary conditions
   *
   * The values of the nodes on each boundary are set to analytic values from
   * the known solution.
   */
  virtual void fillEssBC() {
    initEssBC();  // initialize the boundary conditions before filling them
    // Loop over all the nodes. If the node is a boundary node, calculate the
    // analytical value at that point and set that value as the boundary value.
    for (int node_id = 0; node_id < p_grid_->n_nodes(); node_id++) {
      if (p_grid_->BoNode(node_id)) {  // this is a boundary node
        double u_analytical = analytic_solution_->ValueAt(node_id, t_ + dt_);
        specifyValue(node_id, 0, u_analytical);  // set the boundary value
      }
    }
  }

  /**
   * Solves the system for a given time value.
   *
   * @param dt the time step between solve times
   * @param t the current time value
   */
  virtual void Solve(double delta_t, double current_time) {
    timers_[kTimerSolve].Start();  // we want to time the entire solve process
    this->t_ = current_time;  // store these in the object data for later use
    this->dt_ = delta_t;

    PrintInfo("  Filling essential boundary conditions...");
    fillEssBC();  // Set the boundary conditions
    PrintInfo("  Applying Ess BC to solution...");
    ApplyEssBCToSolution();  // apply boundary conditions to the solution vector
    timers_[kTimerAssemble].Start();  // we're timing just the assembly
    PrintInfo("  Assembling...");
    Assemble(false);  // assemble the system
    timers_[kTimerAssemble].Stop();
    PrintInfo("  Applying Ess BC...");
    ApplyEssBC();  // apply boundary conditions to the assembled system

    timers_[kTimerKSPSolve].Start();  // time the KSPSolve step
    PrintInfo("  Solving with KSP...");
    SolveKSP(solution_, 1, 0);  // run the KSP solver (uses PETSc)
    timers_[kTimerKSPSolve].Stop();

    timers_[kTimerUpdate].Start();  // time the process of saving the solution

    // The result of the system solve is in the solution_ vector. We want to
    // store this data in our node data arrays where we keep track of the
    // current and previous heat values. This function call will put the data
    // from solution_ into location 0 of each NodeData object.
    PrintInfo("  Copying data from solution vector into node data...");
    p_data_->NodeDataFromArray(solution_.data(), U, n_dof());

    // Here we copy the hermite derivatives into each node's HTNodeData.
    // The hermite basis function adds extra degrees of freedom to our equation
    // so it can include derivatives in the solution. n_dof() includes these
    // extra degrees of freedom (+1 for 1D, +3 for 2D, and +7 for 3D).
    // (The exact formula is n_dof() = given_n_dof * 2^nsd).
    //
    // The first value is the u value, which we already copied above.
    // The rest of the values are the derivatives. Since this is a
    // multi-timestep problem, we need to keep these derivatives around
    // so we can use them to calculate valueFEM in Integrands().
    //
    // This loop works because HTNodeData is arranged such that
    // DU = derivative 1, DU + 1 = derivative 2, etc.
    for (int i = 1; i < n_dof(); i++) {
      p_data_->NodeDataFromArray(solution_.data() + i, DU + i - 1, n_dof());
    }

    PrintInfo("  Done!");

    timers_[kTimerUpdate].Stop();
    timers_[kTimerSolve].Stop();
  }

  /**
   * Fills the Ae and be structures with data for a single Gauss point.
   *
   * This uses a Backward Euler scheme to discretize the heat equation.
   *
   * @param fe the element we are assembling a Gauss point from
   * @param Ae the element matrix to put data in
   * @param be the element vector to put data in
   */
  virtual void Integrands(const FEMElm& fe, ZeroMatrix<double>& Ae,
                          ZEROARRAY<double>& be) {
    const int n_dimensions = fe.nsd();  // # of dimensions: 1D, 2D, or 3D
    const int n_basis_functions = fe.nbf();  // # of basis functions
    const double detJxW = fe.detJxW();  // (determinant of J) cross W
    double k_val = input_data_->K_;  // thermal diffusivity in heat equation

    // ValueFEM works for hermite because HTNodeData is ordered so that the
    // DU_PRE values come right after the U_PRE value.
    // (U_PRE, DU_PRE_1, DU_PRE_2, etc.)
    const double u_pre_curr = p_data_->valueFEM(fe, U_PRE);

    // in order to assemble the gauss point, we loop over each pair of basis
    // functions and calculate the individual contributions to the 'Ae' matrix
    // and 'be' vector.
    for (int a = 0; a < n_basis_functions; a++) {
      for (int b = 0; b < n_basis_functions; b++) {
        double M = fe.N(a) * fe.N(b) * detJxW;
        double N = 0;
        for (int k = 0; k < n_dimensions; k++) {
          N += k_val * fe.dN(a, k) * fe.dN(b, k) * detJxW;
        }
        // Add term to the A element matrix
        Ae(a, b) += M / dt_ + N;
      }

      // Add term to the b element vector
      be(a) += fe.N(a) / dt_ * u_pre_curr * detJxW;
    }
  }

  /**
   * Fills only the be structures with data for a single Gauss point.
   *
   * This uses a Backward Euler scheme to discretize the heat equation.
   *
   * Warning: will not work for hermite basis functions.
   *
   * @param fe the element we are assembling a Gauss point from
   * @param be the element vector to put data in
   */
  virtual void IntegrandsVectorOnly(const FEMElm& fe, ZEROARRAY<double>& be) {
    const int n_basis_functions = fe.nbf();  // # of basis functions
    const double detJxW = fe.detJxW();  // (determinant of J) cross W

    // in order to assemble the gauss point, we loop over each pair of basis
    // functions and calculate the individual contributions to the 'be' vector.
    for (int a = 0; a < n_basis_functions; a++) {
      for (int b = 0; b < n_basis_functions; b++) {
        double M = fe.N(a) * fe.N(b) * detJxW;
        // get process id that corresponds to the id of the node in the element
        int node_id = fe.elem()->ElemToLocalNodeID(b);
        // Add term to the b element vector
        be(a) += M / dt_ * p_data_->GetNodeData(node_id).u_pre;
      }
    }
  }

  /**
   * Fills the Ae and be structures with data for the entire element.
   *
   * This uses a Backward Euler scheme to discretize the heat equation.
   *
   * Warning: will not work for hermite basis functions.
   *
   * @param fe the element we are assembling a Gauss point from
   * @param Ae the element matrix to put data in
   * @param be the element vector to put data in
   */
  virtual void IntegrandsByElement(const FEMElm& fe, ZeroMatrix<double>& Ae,
                                   ZEROARRAY<double>& be) {
    const int n_basis_funcs = fe.nbf();
    const double k_val = input_data_->K_;  // thermal diffusivity
    double *q_prev = new double[n_basis_funcs];

    for (int i = 0; i < n_basis_funcs; i++) {
      const int idx = fe.elem()->ElemToLocalNodeID(i);
      q_prev[i] = p_data_->GetNodeData(idx).u_pre;
    }

    // fill be values
    const double dt_inv = 1.0 / dt();
    for (int a = 0; a < n_basis_funcs; a++) {
      for (int b = 0; b < n_basis_funcs; b++) {
        be(a) += ea_NN_(a, b) * dt_inv * q_prev[b];
      }
    }

    // fill Ae values
    for (int a = 0; a < n_basis_funcs; a++) {
      const double M_plus_K = ea_NN_(a, a) * dt_inv + k_val * ea_dNdN_(a, a);
      Ae(a, a) += M_plus_K;
    }
    for (int a = 0; a < n_basis_funcs-1; a++) {
      for (int b = a+1; b < n_basis_funcs; b++) {
        const double M_plus_K = ea_NN_(a, b) * dt_inv + k_val * ea_dNdN_(a, b);
        Ae(a, b) += M_plus_K;
        Ae(b, a) += M_plus_K;
      }
    }
    delete [] q_prev;
  }

  /**
   * Fills only the be structures with data for the entire element.
   *
   * This uses a Backward Euler scheme to discretize the heat equation.
   *
   * Warning: will not work for hermite basis functions.
   *
   * @param fe the element we are assembling a Gauss point from
   * @param be the element vector to put data in
   */
  virtual void IntegrandsByElementVectorOnly(const FEMElm& fe,
                                             ZEROARRAY<double>& be) {
    const int n_basis_funcs = fe.nbf();
    double *q_prev_over_dt = new double[n_basis_funcs];

    const double dt_inv = 1.0 / dt();
    for (int i = 0; i < n_basis_funcs; i++) {
      const int idx = fe.elem()->ElemToLocalNodeID(i);
      q_prev_over_dt[i] = p_data_->GetNodeData(idx).u_pre * dt_inv;
    }

    // fill be values
    for (int a = 0; a < n_basis_funcs; a++) {
      for (int b = 0; b < n_basis_funcs; b++) {
        be(a) += ea_NN_(a, b) * q_prev_over_dt[b];
      }
    }
    delete [] q_prev_over_dt;
  }

  virtual void IntegrandsPreallocator(const FEMElm& fe, ZeroMatrix<bool>& Ae) {
    const int n_basis_functions = fe.nbf();  // # of basis functions
    for (int a = 0; a < n_basis_functions; a++) {
      for (int b = 0; b < n_basis_functions; b++) {
        Ae(a, b) = true;
      }
    }
  }

 private:
  MPITimer timers_[4];  ///< for timing several parts of the code
  HTInputData* input_data_;  ///< pointer to input data
  HTAnalyticSolution *analytic_solution_;
};

#endif  // INCLUDE_HTEQUATION_H_
