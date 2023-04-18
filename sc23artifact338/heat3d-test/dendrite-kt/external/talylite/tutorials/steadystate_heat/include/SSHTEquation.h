#pragma once

#include "SSHTNodeData.h"

#include <talyfem/talyfem.h>

#include "SSHTAnalyticSol.h"

class SSHTEquation : public CEquation<SSHTNodeData> {
 public:
  SSHTEquation(SSHTAnalyticSolution* sol) : analytic_sol_(sol) {}

  void setBC(const std::map<int, SSHTInputData::BoundaryCondition> &bc) {
    bc_ = bc;
  }

  virtual void Solve(double delta_t, double current_time) override {
    fillEssBC();

    // Assemble Ae matrix and be vector
    makeSystem();

    // Calculate the solution vector by using the Krylov subspace Solve
    SolveKSP(solution_, 1, 0);

    // The result of the system solve is in the solution_ vector. We want to
    // store this data in our NodeData arrays so we can save and load the data
    // later. This helper function puts the data from solution_  into location 0
    // of each NodeData object.
    //p_data_->NodeDataFromArray(solution_.data(), 0);
    p_data_->NodeDataFromArray(solution_.data(), TEMPERATURE_IDX, n_dof());
  }

  // Specify the Dirichlet boundary
  virtual void fillEssBC() override {
    initEssBC();

    // if the node is on *any* boundary do not specify a bundary parameter when calling BoNode
    // to check a specific boundary, use BoNode(node_id, boudary id)
    /* boundary ids
        1 == -ve x
        2 == +ve x
        3 == -ve y
        4 == +ve y */

    // At least one dirichlet boundary condition should be set, or there are
    // multiple valid solutions for this problem. This variable is used to
    // warn the user if that is the case.
    bool found_dirichlet = false;

    // loop over all the nodes...
    for (LocalNodeID node_id = 0; node_id < p_grid_->n_nodes(); node_id++) {
      // loop over all the boundary conditions
      // it->first = the boundary ID (input to BoNode)
      // it->second = the boundary type (DIRICHLET/NEUMANN)
      for (auto it = bc_.begin(); it != bc_.end(); it++) {
        // does this node lie in the boundary?
        if (p_grid_->BoNode(node_id, it->first)) {
          // yes it does, is this a dirichlet boundary?
          if (it->second == SSHTInputData::DIRICHLET) {
            // yes it is, fix the value as the analytic solution at that point
            const ZEROPTV& pt = p_grid_->GetNode(node_id)->location();
            double val = analytic_sol_->calc_u_at(pt);
            specifyValue(node_id, TEMPERATURE_IDX, val);
            found_dirichlet = true;
          }
          // if it->second == NEUMANN, we handle this in Integrands4side
        }
      }
    }

    // check if any process set found_dirichlet to true
    bool any_dirichlet = false;
    MPI_Allreduce(&found_dirichlet, &any_dirichlet, 1,
                  MPI_CXX_BOOL, MPI_LOR, PETSC_COMM_WORLD);
    if (!any_dirichlet) {
      PrintWarning("No dirichlet boundaries were set. Solution is non-unique.");
      PrintWarning("If you are loading from a file, make sure ifLoadNodeIndicators is set ",
                   "and your mesh has the appropriate node indicators (1-6).");
    }
  }

  virtual void Integrands(const FEMElm &fe, ZeroMatrix<double> &Ae,
                          ZEROARRAY<double> &be) override {
    const int n_dimensions = fe.nsd();      // # of dimensions: 1D, 2D, or 3D
    int n_basis_functions = fe.nbf();  // # of basis functions
    const double detJxW = fe.detJxW();      // (determinant of J) cross W
    const ZEROPTV p = fe.position();

    double force = analytic_sol_->calc_d2u_at(p);

    for (int a = 0; a < n_basis_functions; a++) {
      for (int b = 0; b < n_basis_functions; b++) {
        double N = 0;
        for (int k = 0; k < n_dimensions; k++) {
          N += fe.dN(a, k) * fe.dN(b, k) * detJxW;
        }
        Ae(a, b) -= N;
      }
      be(a) += fe.N(a) * (force) * detJxW;
    }
  }

  // Defining Neumann boundary condition of the form \surf_int w*(grad (u))dA
  // which is a surface integral, where A is the surface area
  virtual void Integrands4side(const FEMElm &fe, int sideInd, ZeroMatrix<double> &Ae, ZEROARRAY<double> &be) override {
    auto it = bc_.find(sideInd);
    if (it != bc_.end() && it->second == SSHTInputData::NEUMANN) {
      const ZEROPTV& p = fe.position();
      const ZEROPTV& normal = fe.surface()->normal();
      double flux = analytic_sol_->calc_grad_u_at(p).innerProduct(normal);

      for (int i = 0; i < fe.nbf(); i++) {
        be(i) -= fe.N(i) * fe.detJxW() * flux;
      }
    }
  }

 private:
  // map side index to boundary condition type
  std::map<int, SSHTInputData::BoundaryCondition> bc_;
  SSHTAnalyticSolution* analytic_sol_;
};
