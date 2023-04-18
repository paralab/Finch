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

#include <talyfem/talyfem.h>


/**
 * Class to calculate the volume of the system.
 *
 * This is a simple example usage of the FunctionIntegrator class.
 */
class VolumeFunction : public FunctionIntegrator {
 public:
  /**
   * Construct the volume integrator object.
   *
   * @param[in] basis_function basis function to use
   * @param[in] rel_order relative order to use
   * @param[in] do_accelerate whether to use acceleration
   * @param[in] method method of integration (see IntegrationMethod for options)
   */
  VolumeFunction(kBasisFunction basis_function, int rel_order,
                 bool do_accelerate = false,
                 AssemblyMethod method = kAssembleGaussPoints);

  /**
   * Calculate the volume value at the Gauss point.
   *
   * @param[in] fe the finite element with the appropriate Gauss point selected
   * @return the value of the volume function at the gauss point
   */
  virtual double CalcGaussPointIntegral(FEMElm& fe) override;
};
