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
#include <volume_function.h>


VolumeFunction::VolumeFunction(kBasisFunction basis_function, int rel_order,
                               bool do_accelerate, AssemblyMethod method)
    : FunctionIntegrator(rel_order, do_accelerate, method) { }

double VolumeFunction::CalcGaussPointIntegral(FEMElm& fe) {
  double value = 0.0;
  for (int a = 0, max_a = fe.nbf(); a < max_a; a++) {
    value += fe.N(a);
  }
  return value * fe.detJxW();
}
