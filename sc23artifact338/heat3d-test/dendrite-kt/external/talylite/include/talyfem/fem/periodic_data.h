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
#ifndef FEM_PERIODIC_DATA_H_
#define FEM_PERIODIC_DATA_H_

#include <algorithm>  // for std::count
#include <vector>  // for std::vector


namespace TALYFEMLIB {

class PeriodicBounds;

/**
 * Periodic data details to pass to CEquation classes.
 *
 * This class is responsible for describing the periodic state of an equation
 * to be solved. The class stores:
 * - the underlying periodic grid structure
 * - which variables are periodic
 * - whether to set up a periodic exchange across the boundaries
 *
 * This class requires an underlying PeriodicBounds object to be created.
 * The same PeriodicBounds object can be used for multiple equation solvers.
 * Each solver will need its own PeriodicData object. By using a different
 * object for each equation to be solved, it is possible to have each equation
 * use a different set of periodic conditions on the same underlying mesh.
 *
 * See the PeriodicBounds class for details on how to create a PeriodicBounds
 * object. Once that is created, the PeriodicData object is created and set up
 * with the list of which variables are periodic.
 *
 * Once this object is created, it is passed to the equation solver.
 * This sets up periodic bounds and should now automatically apply periodic
 * bounds. By default, the CEquation solver will apply the boundary conditions
 * along with the essential boundary conditions. It will also automatically
 * remap the values of periodic partners as needed after the solve is complete.
 *
 * Note: once this object has been passed to the solver, changes to this object
 * will have no effect on the solver.
 *
 * Example usage (1 periodic degree of freedom, no exchange):
 * // using previously defined PeriodicBounds object
 * const int n_total_dof = 1;  // one degree of freedom
 * PeriodicData periodic_data(periodic_bounds, int n_total_dof);
 * // The first variable (index 0) is periodic
 * periodic_data.SetVarIndexPeriodic(0);
 * // pass to CEquation object, along with previously defined p_grid object
 * equation.redim(p_grid, int n_total_dof, &periodic_data);
 *
 * Example usage (3 degree of freedom, 1st and 3rd periodic, with exchange):
 * // using previously define PeriodicBounds object
 * const int n_total_dof = 3;  // three degrees of freedom
 * PeriodicData periodic_data(periodic_bounds, n_total_dof);
 * // The first variable (index 0) and third variable (index 2) are periodic
 * periodic_data.SetVarIndexPeriodic(0);
 * periodic_data.SetVarIndexPeriodic(2);
 * periodic_data.set_enable_exchange(true);  // turn on exhange support
 * // pass to CEquation object, along with previously defined p_grid object
 * equation.redim(p_grid, n_total_dof, &periodic_data);
 */
class PeriodicData {
 public:
  /**
   * Creates the periodic data object with the given information
   *
   * This sets up the the data needed to make a solver periodic.
   * It copies the relevant data from the PeriodicBounds object to this object.
   * It also sets up tracking of which variables are periodic. By default none
   * of the variables are periodic.
   *
   * In order to actually make the system periodic, at least one solver
   * variable must be set to periodic using SetVarIndexPeriodic().
   *
   * @param periodic_bounds pointer to underlying periodic data
   * @param n_degrees_of_freedom number of degrees of freedom in the system
   */
  PeriodicData(PeriodicBounds *periodic_bounds, int n_degrees_of_freedom);

  /**
   * Marks given variable index as periodic and updates object accordingly.
   *
   * This function actually makes a given variable periodic. By default all
   * variables are non-periodic until marked as periodic by this function.
   * Trying to set a variable as periodic twice will lead to an exception.
   *
   * @param index index to mark as periodic
   * @throw TALYException if trying to set a variable that is already periodic
   */
  void SetVarIndexPeriodic(int index);

  /**
   * Marks given variable index as nonperiodic and updates object accordingly.
   *
   * Trying to mark a nonperiodic variable as nonperiodic again will lead to an
   * exception. When the object is first created, all variables are nonperiodic.
   *
   * @param index index to mark as nonperiodic
   * @throw TALYException if trying to mark a variable that is already
   *                      nonperiodic
   */
  void SetVarIndexNonPeriodic(int index);

  /**
   * Toggles whether to request setting up a periodic exchange object.
   *
   * @param enable_exchange_value whether to set up periodic exchange object
   */
  inline void set_enable_exchange(bool enable_exchange_value) {
    enable_exchange_ = enable_exchange_value;
  }

  /**
   * Returns a bool specifying whether exchange is enabled.
   */
  inline bool enable_exchange() const { return enable_exchange_; }

  /**
   * @returns whether the data is periodic.
   *
   * The data is considered periodic if at least one variable has been set
   * to be periodic
   */
  inline bool is_periodic() const { return is_periodic_; }

  /**
   * Returns vector of periodic variables
   *
   * @return vector of periodic variables
   */
  inline const std::vector<int> periodic_var_list() const {
    return periodic_var_list_;
  }

  /**
   * Returns number of periodic variables
   *
   * @return number of periodic variables
   */
  inline int NumPeriodicVars() const {
    return static_cast<int>(periodic_var_list_.size());
  }

  /**
   * Returns vector specifying whether each variable is periodic
   *
   * @return vector specifying whether each variable is periodic
   */
  inline const std::vector<bool> is_var_periodic() const {
    return is_var_periodic_;
  }

  /**
   * Returns a pointer to the underlying periodic bounds object
   */
  inline PeriodicBounds *periodic_bounds() const {
    return periodic_bounds_;
  }

 private:
  PeriodicBounds *periodic_bounds_;  ///< periodic mapping data
  std::vector<int> periodic_var_list_;  ///< list of variables that are periodic
  std::vector<bool> is_var_periodic_;  ///< whether a given variable is periodic
  bool is_periodic_;  ///< whether this data is actually periodic
  bool enable_exchange_;  ///< whether to allow exhange of periodic values
};

}  // namespace TALYFEMLIB

#endif  // FEM_PERIODIC_DATA_H_
