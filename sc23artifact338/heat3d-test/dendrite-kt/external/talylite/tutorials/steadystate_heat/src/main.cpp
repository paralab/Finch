#include <talyfem/talyfem.h>

#include <SSHTInputData.h>
#include <SSHTAnalyticSol.h>
#include <SSHTEquation.h>
#include <SSHTNodeData.h>

int main(int argc, char **args)
{
  PetscInitialize(&argc, &args, NULL, NULL); // set up PETSc environment

  // This try/catch statement has two uses:
  // 1. It reports errors instead of just aborting when an exception is thrown.
  // 2. It makes sure MyEquation goes out of scope before PetscFinalize()
  //    is called. Otherwise, the destructors associated with the solver may
  //    call something PETSc-related, which causes a crash.
  try
  {
    SSHTInputData input_data;

    if (!input_data.ReadFromFile()) {
      throw TALYException() << "Error reading input data, check config.txt!";
    }

    if (!input_data.CheckInputData()) {
      throw TALYException() << "Invalid input data, check config.txt!";
    }

    // Create a GRID depending on options set in config.txt.
    GRID *p_grid = NULL;             // Declaring pointer of type GRID to store grid data
    CreateGrid(p_grid, &input_data); // fills in p_grid with a GRID based on input_data

    /*Matrix3 mat = Matrix3::axisAngle(ZEROPTV(0, 0, 1), -M_PI / 4);
    for (int i = 0; i < p_grid->n_nodes(); i++) {
      ZEROPTV& loc = p_grid->GetNode(i)->location();
      loc = mat * loc;
    }*/

    SSHTAnalyticSolution analytic_sol(input_data.analytic_sol_type);

    // Create our equation object.
    SSHTEquation equation(&analytic_sol);
    equation.add_basis_flag(BASIS_DIMENSION_REDUCTION);  // support 1D/2D in 3D elements
    equation.SetPreallocator(new PreallocatorPerfect<SSHTNodeData>(&equation));

    // This allocates the matrix and vector (stored as part of the equation object).
    int n_dof = 1;                                                                                // number of degrees of freedom
    equation.redimSolver(p_grid, n_dof, false, input_data.basisRelativeOrder); // resize the equation to match the grid structure

    // Create our GridField.
    //The GridField holds the data for each node in our mesh.
    GridField<SSHTNodeData> grid_field;
    grid_field.redimGrid(p_grid);  // tell the gridfield where to find the grid (for number of nodes)
    grid_field.redimNodeData();    // initialize the node data
    equation.setData(&grid_field); // tell the equation where to get node data

    equation.setBC(input_data.boundary_conditions);

    // Only for transient problems //
    // Solve a single timestep of our equation by calling MyEquation::Solve().
    double t = 0.0;
    double dt = 0.0;
    //    

    //You do still need to give something to solve method (even if you are solving steady state problems) as it has been written in CEquation.h
    equation.Solve(t, dt);

    /*//----------------------------To fill in the analytic solution--------------------------//
                                  Define the analytic solution equation
                                  For laplacian{T} = 0 all non-dimension

                                  //--------------For 1D case ---------------//
                                Both walls dirichlet T_left = 1 ; T_right  = 0;
                                                T_exact = 1 - x;
    //--------------------------------------------------------------------------------------*/

    for (int i = 0; i < p_grid->n_nodes(); i++)
    {
      SSHTNodeData& data = grid_field.GetNodeData(i);  // MyNodeData& type is the reference to nodedata 
                                                     // from GridField for getting data at each node
      NODE* node = p_grid->GetNode(i);  // Grab the node associated with index i 
      ZEROPTV p = node->location();  // Grab the location associated with the node
      data.T_exact = analytic_sol.calc_u_at(p);
    }

    // ----------------------------------------------------------------------------------------//

    {
      std::stringstream ss;
      ss << "data_final" << input_data.output_extension;
      save_gf(&grid_field, &input_data, ss.str().c_str(), 0.0);
    }

    double l2_err = analytic_sol.calc_l2_error(&input_data, p_grid, &grid_field);

    std::stringstream ss;
    ss << "L2 error = " << l2_err;
    PrintResults(ss.str(), l2_err < 1e-2, false);

    // Free the GRID object.
    DestroyGrid(p_grid);
  }
  catch (TALYException &e)
  {
    PrintError(e.what());
  }

  PetscFinalize(); // clean up PETSc environment
  return 0;        // program exits normally
}
