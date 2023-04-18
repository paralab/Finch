/*
  Copyright 2017 Baskar Ganapathysubramanian

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

#include <vector>

#include <talyfem/fem/cequation.h>
#include <talyfem/basis/basis.h>
#include <talyfem/basis/elemnodes.h>
#include <talyfem/grid/femelm.h>


namespace TALYFEMLIB {

/**
 * Contains FunctionIntegrator, the base class for a simple function integrator.
 *
 * The FunctionIntegrator class is used to integrate a function over the entire
 * system domain. This class needs to be extended with functions to calculate
 * the desired integral on a gauss point or elemental based.
 *
 * The class must define either CalcGaussPointIntegral(FEMElm& fe) or
 * CalcGaussPointIntegral(FEMElm& fe) to provide the integral function.
 *
 * An example class to calculate the volume of a system using the Gauss
 * point method would be:
 *
 * class VolumeFunction : public FunctionIntegrator {
 *  public:
 *   VolumeFunction(kBasisFunction basis_function, int rel_order,
 *                  bool do_accelerate=false,
 *                  IntegrationMethod method=FunctionIntegrator::kByGaussPoint)
 *       : FunctionIntegrator(basis_function, rel_order, do_accelerate, method)
 *   { }
 *
 *   virtual double VolumeFunction::CalcGaussPointIntegral(FEMElm& fe) {
 *     double value = 0;
 *     for (int a = 1; a <= fe.elem()->n_nodes(); a++) {
 *       value += fe.N(a);
 *     }
 *     return value * fe.detJxW();
 *   }
 * }
 *
 * Usage of this example class:
 * VolumeFunction volume_integrator(basis_function, rel_order);
 * volume_integrator.set_pGrid(p_grid);
 * double volume = volume_integrator.Solve();
 *
 */


/**
 * Base class for a simple function integrator.
 *
 * This does not support boundary conditions.
 * See file comment above for example usage.
 */
class FunctionIntegrator {
 public:
  /**
   * Constructs the object
   *
   * @param[in] rel_order relative order to use
   * @param[in] do_accelerate whether to use acceleration
   * @param[in] method method of integration. Either kAssembleGaussPoints or
   *            kAssembleElements
   */
  FunctionIntegrator(int rel_order, bool do_accelerate, AssemblyMethod method);

  /**
   * Constructs the object, specifying volume
   *
   * @param[in] rel_order relative order to use
   * @param[in] volume the system volume
   * @param[in] do_accelerate whether to use acceleration
   * @param[in] method method of integration. Either kAssembleGaussPoints or
   *            kAssembleElements
   */
  FunctionIntegrator(int rel_order, double volume, bool do_accelerate, AssemblyMethod method);

  virtual ~FunctionIntegrator() { }

  /**
   * Calculates the integral averaged over the volume
   *
   * @return the volume averaged integral
   */
  double CalcVolumeAverage();

  /**
   * Calculates the integral at a Gauss point.
   *
   * This needs to be dervied by a subclass.
   *
   * @param[in] fe the finite element with the appropriate Gauss point selected
   * @return the value of the function at the gauss point
   */
  virtual double CalcGaussPointIntegral(FEMElm& fe);

  /**
   * Calculates the value of integral for an element.
   *
   * This needs to be dervied by a subclass.
   *
   * @param[in] fe the finite element to integrate over
   * @return the value of the function over the element
   */
  virtual double CalcElementIntegral(FEMElm& fe);

  /**
   * Calculates the value of integral over the entire system
   *
   * This function calculates the total value by taking a sum of the values
   * from all processes.
   * Note: this function is intended to be overloaded to perform any set up
   * needed prior to doing the integration.
   *
   * @return the global value of the function
   */
  virtual double Solve();

  /**
   * Sets the grid pointer to the system grid.
   *
   * If the object is using acceleration, this function also sets up the
   * acceleration parameters for later use.
   *
   * @param[in] pGrid pointer to grid
   */
  void set_pGrid(GRID* pGrid);

  /**
   * Sets the volume of the system.
   *
   * @param[in] volume the total volume of the system
   */
  inline void set_volume(double volume) {
    volume_ = volume;
    volume_set_ = true;
  }

 protected:
  /**
   * Calculate the value of integral over the entire system.
   * This function is called by Solve. Depending on the initalization, this
   * will either integrate over Gauss points or elementsIt calculates the total
   * value by taking a sum of the values from all processes.
   *
   * @return the global value of the function
   */
  double Integrate();

  GRID* pGrid_;  ///< pointer to system grid
  int rel_order_;  ///< relative order of integration

 private:
  /**
   * Set up the data structures for integration when acceleration is in use.
   * This is called by set_pGrid.
   */
  void FillAccelerateElements();

  /**
   * Calculate system integral using elemental integration.
   * This function is called by Integrate.
   *
   * @return the global value of the function
   */
  double IntegrateByElement();

  /**
   * Calculate system integral using Gauss point integration.
   * This function is called by Integrate.
   *
   * @return the global value of the function
   */
  double IntegrateByGaussPoint();

  double volume_;  ///< volume of the system
  bool volume_set_;  ///< whether the volume has been set
  bool do_accelerate_;  ///< whether to use acceleration
  AssemblyMethod integrate_method_;  ///< method of integration
  std::vector<FEMElm> fe_accelerate_;  ///< finite elements for acceleration
};

}  // namespace TALYFEMLIB
