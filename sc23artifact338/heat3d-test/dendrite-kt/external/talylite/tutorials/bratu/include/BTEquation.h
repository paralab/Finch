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
#ifndef BT_EQUATION_HPP
#define BT_EQUATION_HPP

#include "BTNodeData.h"


MPITimer timers[4];  ///< for timing several parts of the code

// indices to timing array
const int kTimerSolve = 0;
const int kTimerAssemble = 1;
const int kTimerSNESSolve = 2;
const int kTimerUpdate = 3;

// need to track the first assemble and last update time
// see destructor for details.
double first_assemble;
double last_update;

class BTEquation : public CEquation<BTNodeData> {
 public:
  BTInputData* idata_;  ///< pointer to inputdata

  SNES snes;  ///< nonlinear solver context
  int counter;  ///< number of function calls
  double xl, xr, Lx;

  BTEquation(BTInputData* idata, bool has_uniform_mesh = false,
             AssemblyMethod assembly_method = kAssembleGaussPoints);
  virtual ~BTEquation();
  virtual void fillEssBC();
  void Solve(double delta_t, double current_time);
  virtual void Integrands(const FEMElm& fe, ZeroMatrix<double>& Ae,
                          ZEROARRAY<double>& be);
  virtual void IntegrandsByElement(const FEMElm& fe, ZeroMatrix<double>& Ae,
                                   ZEROARRAY<double>& be);
};

PetscErrorCode UpdateGridField(Vec _xg, BTEquation *ceqn, int flag = 0);
PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ceqn);
PetscErrorCode FormJacobian(SNES snes, Vec _xg, Mat jac, Mat B, void *ceqn_);

BTEquation::BTEquation(BTInputData* idata, bool has_uniform_mesh,
                       AssemblyMethod equation_assembly_method)
    : CEquation<BTNodeData>(has_uniform_mesh, equation_assembly_method),
      idata_(idata),
      xl(-0.4),
      xr(0.4),
      Lx(1) {
  counter = 0;
  timers[kTimerSolve].set_label("Solve");
  timers[kTimerAssemble].set_label("Assemble");
  timers[kTimerSNESSolve].set_label("SNESSolve");
  timers[kTimerUpdate].set_label("Update");
  first_assemble = 0.0;
  last_update = 0.0;
}

BTEquation::~BTEquation() {
  timers[kTimerSolve].PrintGlobalTotalSeconds();
  timers[kTimerAssemble].PrintGlobalTotalSeconds();
  // SNES solve time includes assemble and update times. To get an accurate
  // total, these times need to be removed from the SNESSolve total time.
  // The first assemble value and last update value are not part of the
  // SNESSolve total and should not be subtracted from it.
  double assemble_mod = (timers[kTimerAssemble].GetTotalTimeSeconds()
      - first_assemble) * -1.0;
  double update_mod = (timers[kTimerUpdate].GetTotalTimeSeconds() - last_update)
      * -1.0;
  timers[kTimerSNESSolve].AddToTotalTime(assemble_mod);
  timers[kTimerSNESSolve].AddToTotalTime(update_mod);
  timers[kTimerSNESSolve].PrintGlobalTotalSeconds();
  timers[kTimerUpdate].PrintGlobalTotalSeconds();
}

void BTEquation::fillEssBC() {
  // Setting boundary conditions. for 1D problem
  this->initEssBC();
  for (int nodeID = 0; nodeID < this->p_grid_->n_nodes(); nodeID++) {
    if (this->p_grid_->BoNode(nodeID, 1)) {
      specifyValue(nodeID, 0, 0);
      this->p_data_->GetNodeData(nodeID).u[0] = -0.668371;
    }
    if (this->p_grid_->BoNode(nodeID, 2)) {
      specifyValue(nodeID, 0, 0);
      this->p_data_->GetNodeData(nodeID).u[0] = 3.017089;
    }
  }
}

void BTEquation::Solve(double delta_t = 0.0, double current_time = 0.0) {
  timers[kTimerSolve].Start();
  PetscErrorCode ierr;

  // Assemble the Jacobian matrix and Function vector
  this->fillEssBC();
  timers[kTimerAssemble].Start();
  this->Assemble();
  timers[kTimerAssemble].Stop();
  first_assemble = timers[kTimerAssemble].GetLastTimeSeconds();
  this->ApplyEssBC();

  // Create nonlinear solver context
  ierr = SNESCreate(PETSC_COMM_WORLD, &snes);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // Set initial guess while checking if mesh
  // partitioning is used i.e. parallel_type_==kWithDomainDecomp
  if (this->p_grid_->parallel_type_ == kWithDomainDecomp) {
    PetscScalar *array;
    ierr = PetscMalloc(this->n_total_dof() * sizeof(PetscScalar), &array);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    Vec initVec;
    VecScatter scatter;
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, this->n_total_dof(),
                                 array, &initVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    // setting initial guess
    for (int nodeID = 0; nodeID < this->p_grid_->n_nodes(); nodeID++) {
      for (int i = 0; i < n_dof(); i++) {
        array[nodeID * n_dof() + i] = this->p_data_->GetNodeData(nodeID).u[i];
      }
    }
    ierr = VecScatterCreate(initVec, this->to(), this->xg_, this->from(), &scatter);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = VecScatterBegin(scatter, initVec, this->xg_, INSERT_VALUES,
                           SCATTER_FORWARD);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = VecScatterEnd(scatter, initVec, this->xg_, INSERT_VALUES,
                         SCATTER_FORWARD);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    ierr = VecDestroy(&initVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    ierr = VecScatterDestroy(&scatter);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = PetscFree(array);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  } else {
    double *val = new double[this->n_dof()];
    PetscInt *index = new PetscInt[this->n_dof()];
    for (int nodeID = 0; nodeID < this->p_grid_->n_nodes(); nodeID++) {
      // setting initial guess
      for (int i = 0; i < n_dof(); i++) {
        val[i] = this->p_data_->GetNodeData(nodeID).u[i];
        index[i] = nodeID * n_dof() + i;
      }
      ierr = VecSetValues(this->xg_, this->n_dof(), index, val, INSERT_VALUES);
      CHKERRABORT(PETSC_COMM_WORLD, ierr);
    }
    delete[] val;
    delete[] index;
    ierr = VecAssemblyBegin(this->xg_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = VecAssemblyEnd(this->xg_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  }

  // MatView(this->Ag,PETSC_VIEWER_STDOUT_SELF);

  // Set function evaluation routine and vector.
  ierr = SNESSetFunction(snes, this->bg_, FormFunction, this);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // Set jacobian matrix
  ierr = SNESSetJacobian(snes, this->Ag_, this->Ag_, FormJacobian, this);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // Set nonlinear solver tolerances
  double atol = 1e-12, rtol = 1e-12, stol = 1e-12;
  int maxit = 30, maxf = 1000;
  counter = 0;  // Resetting function counter
  SNESSetTolerances(snes, atol, rtol, stol, maxit, maxf);
  // SNESMonitorSet(snes,SNESMonitorDefault,PETSC_NULL,PETSC_NULL);
  SNESSetFromOptions(snes);
  timers[kTimerSNESSolve].Start();
  ierr = SNESSolve(snes, PETSC_NULL, this->xg_);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  timers[kTimerSNESSolve].Stop();

  // Get number of iterations used
  PetscInt its;
  SNESGetIterationNumber(snes, &its);
  PetscPrintf(PETSC_COMM_WORLD, "number of iterations = %d\n", its);

  // Store solution to p_data_
  timers[kTimerUpdate].Start();
  UpdateGridField(this->xg_, this);
  timers[kTimerUpdate].Stop();
  last_update = timers[kTimerUpdate].GetLastTimeSeconds();

  ierr = SNESDestroy(&snes);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // Calculate residual
  // ~ double resX=0; VecNorm(this->bg_,NORM_2,&resX);
  // ~ return resX;
  timers[kTimerSolve].Stop();
}

void BTEquation::Integrands(const FEMElm& fe, ZeroMatrix<double>& Ae,
                            ZEROARRAY<double>& be) {
  const int spatial_dims = fe.nsd();
  const int nbf = fe.nbf();
  const double detJxW = fe.detJxW();
  const double lambda = -M_PI * M_PI;

  ZEROPTV du;

  double u = this->p_data_->valueFEM(fe, 0);
  const double lambda_exp_u = lambda * exp(u);

  for (int i = 0; i < spatial_dims; i++) {
    du(i) = this->p_data_->valueDerivativeFEM(fe, 0, i);
  }

  double N;
  for (int a = 0; a < nbf; a++) {
    for (int b = 0; b < nbf; b++) {
      N = 0;
      for (int k = 0; k < spatial_dims; k++) {
        N += fe.dN(a, k) * fe.dN(b, k);
      }
      Ae(a, b) += (N - lambda_exp_u * fe.N(a) * fe.N(b)) * detJxW;
    }

    double du_dw = 0;
    for (int i = 0; i < spatial_dims; i++) {
      du_dw += du(i) * fe.dN(a, i);
    }
    be(a) += -(-du_dw + lambda_exp_u * fe.N(a)) * detJxW;
  }
}

void BTEquation::IntegrandsByElement(const FEMElm& fe, ZeroMatrix<double>& Ae,
                                    ZEROARRAY<double>& be) {
  
  const int nbf = fe.nbf();
  const int spatial_dims = fe.nsd();
  const double lambda = -M_PI * M_PI;

  ZEROPTV du;
  // pulling this part out results in the majority of the the performance
  // improvement. The other parts would be harder to pull out because of the
  // dependence on exp(u).
  for (int a = 0; a < nbf; a++) {
    for (int b = 0; b < nbf; b++) {
      const double N = ea_dNdN_(a,b);
      Ae(a, b) += N;
    }
  }

  for (int g = 0; g < nbf; g++) {
    feAccelerate_[g].set_elem_hack(fe.elem());
    const double detJxW = feAccelerate_[g].detJxW();

    const double u = this->p_data_->valueFEM(feAccelerate_[g], 0);
    const double lambda_exp_u = lambda * exp(u);

    for (int i = 0; i < spatial_dims; i++) {
      du(i) = this->p_data_->valueDerivativeFEM(feAccelerate_[g], 0, i);
    }

    for (int a = 0; a < nbf; a++) {
      const double Na_detJxW = feAccelerate_[g].N(a) * detJxW;
      for (int b = 0; b < nbf; b++) {
        const double M = feAccelerate_[g].N(b) * Na_detJxW;
        Ae(a, b) -= lambda_exp_u * M;
      }

      double du_dw = 0.0;
      for (int i = 0; i < spatial_dims; i++) {
        du_dw += du(i) * feAccelerate_[g].dN(a, i);
      }
      be(a) += -(-du_dw * detJxW + lambda_exp_u * Na_detJxW);
    }
  }
}

PetscErrorCode UpdateGridField(Vec _xg, BTEquation *ceqn, int flag) {
  PetscErrorCode ierr;

  // collect solution from 'xg_' and store into 'solution'
  if (ceqn->p_grid_->parallel_type_ == kWithDomainDecomp) {
    VecScatter scatter;
    Vec SolutionVec;
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, ceqn->n_total_dof(),
                                 ceqn->solution_.data(), &SolutionVec);
    CHKERRQ(ierr);
    ierr = VecScatterCreate(_xg, ceqn->from(), SolutionVec, ceqn->to(), &scatter);
    CHKERRQ(ierr);
    ierr = VecScatterBegin(scatter, _xg, SolutionVec, INSERT_VALUES,
                           SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterEnd(scatter, _xg, SolutionVec, INSERT_VALUES,
                         SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecDestroy(&SolutionVec);
    CHKERRQ(ierr);

    ierr = VecScatterDestroy(&scatter);
    CHKERRQ(ierr);
  } else {
    VecScatter scatter;
    Vec solution_vec;
    ierr = VecScatterCreateToAll(_xg, &scatter, &solution_vec);
    CHKERRQ(ierr);
    ierr = VecScatterBegin(scatter, _xg, solution_vec, INSERT_VALUES,
                           SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterEnd(scatter, _xg, solution_vec, INSERT_VALUES,
                         SCATTER_FORWARD);
    CHKERRQ(ierr);
    double* array;
    ierr = VecGetArray(solution_vec, &array);
    CHKERRQ(ierr);
    memcpy(ceqn->solution_.data(), array,
           sizeof(double) * (ceqn->n_total_dof()));
    ierr = VecRestoreArray(solution_vec, &array);
    CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter);
    CHKERRQ(ierr);

    ierr = VecDestroy(&solution_vec);
    CHKERRQ(ierr);
  }

  // Pass ceqn->solution_ to the p_data_ structure
  for (int i = 0; i < ceqn->p_grid_->n_nodes(); i++) {
    for (int j = 0; j < ceqn->n_dof(); j++) {
      ceqn->p_data_->GetNodeData(i).u[j] = ceqn->solution_(i * ceqn->n_dof() + j);
    }
  }

  return ierr;
}

PetscErrorCode FormFunction(SNES snes, Vec _xg, Vec _f, void *ceqn_) {
  PetscErrorCode ier;
  BTEquation *ceqn = (BTEquation*) ceqn_;
  ceqn->counter++;

  // Computes Ag and bg_
  timers[kTimerUpdate].Start();
  UpdateGridField(_xg, ceqn);
  timers[kTimerUpdate].Stop();

  ceqn->fillEssBC();
  timers[kTimerAssemble].Start();
  ceqn->Assemble(false);
  timers[kTimerAssemble].Stop();
  ceqn->ApplyEssBC();

  ier = VecCopy(ceqn->bg(), _f);
  CHKERRQ(ier);
  return 0;
}

PetscErrorCode FormJacobian(SNES snes, Vec _xg, Mat jac, Mat B,
                            void *ceqn_) {
  //*flag = SAME_NONZERO_PATTERN;
  return 0;
}

#endif
