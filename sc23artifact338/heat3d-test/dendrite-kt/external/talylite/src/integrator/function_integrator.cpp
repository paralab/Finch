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
#include <talyfem/integrator/function_integrator.h>

#include <talyfem/common/comm.h>  // for DoGlobalSum


namespace TALYFEMLIB {

FunctionIntegrator::FunctionIntegrator(int rel_order, bool do_accelerate, AssemblyMethod method)
    : pGrid_(NULL),
      rel_order_(rel_order),
      volume_(-1.0),
      volume_set_(false),
      do_accelerate_(do_accelerate),
      integrate_method_(method) { }

FunctionIntegrator::FunctionIntegrator(int rel_order, double volume, bool do_accelerate, AssemblyMethod method)
    : pGrid_(NULL),
      rel_order_(rel_order),
      volume_(volume),
      volume_set_(true),
      do_accelerate_(do_accelerate),
      integrate_method_(method)  { }

double FunctionIntegrator::CalcVolumeAverage() {
  if (!volume_set_) {
    PrintError("Volume is not set properly in CalcVolumeAverage");
    exit(1);
  }
  return Integrate() / volume_;
}

void FunctionIntegrator::set_pGrid(GRID* pGrid) {
  pGrid_ = pGrid;
  if (do_accelerate_) {  // set up acceleration, if needed
    FillAccelerateElements();
  }
}

double FunctionIntegrator::Integrate() {
  switch (integrate_method_) {
    case kAssembleGaussPoints:
      return IntegrateByGaussPoint();
    case kAssembleElements:
      return IntegrateByElement();
  }
  PrintError("Unknown integration method");
  exit(1);
}

double FunctionIntegrator::IntegrateByGaussPoint() {
  FEMElm fe(pGrid_);  // finite element for calculating integral
  double value = 0.0;
  if (do_accelerate_) {
    for (int elmID = 0; elmID < pGrid_->n_elements(); elmID++) {
      if (!(pGrid_->IsMyElement(elmID))) continue;
      for (std::size_t g = 0, max = fe_accelerate_.size(); g != max; ++g) {
        fe_accelerate_[g].set_elem_hack(elmID);
        value += CalcGaussPointIntegral(fe_accelerate_[g]);
      }
    }
  } else {  // no acceleration
    for (int elmID = 0; elmID < pGrid_->n_elements(); elmID++) {
      if (!(pGrid_->IsMyElement(elmID))) continue;
      fe.refill(elmID, rel_order_);
      while (fe.next_itg_pt()) {
        value += CalcGaussPointIntegral(fe);
      }
    }
  }
  DoGlobalSum(value, 1);
  return value;
}

double FunctionIntegrator::IntegrateByElement() {
  FEMElm fe(pGrid_);  // finite element for calculating integral
  double value = 0.0;
  for (int elmID = 0; elmID < pGrid_->n_elements(); elmID++) {
    if (!(pGrid_->IsMyElement(elmID))) continue;
    fe.refill(elmID, rel_order_);
    value += CalcElementIntegral(fe);
  }
  DoGlobalSum(value, 1);
  return value;
}

double FunctionIntegrator::Solve() {
  return Integrate();
}

double FunctionIntegrator::CalcGaussPointIntegral(FEMElm& fe) {
  PrintError("CalcGaussPointIntegral not implimented");
  exit(1);
}

// default implementation is the same as the gauss point integration.
double FunctionIntegrator::CalcElementIntegral(FEMElm& fe) {
  double value = 0.0;
  while (fe.next_itg_pt()) {
    value += CalcGaussPointIntegral(fe);
  }
  return value;
}

void FunctionIntegrator::FillAccelerateElements() {
  FEMElm fe(pGrid_);

  int elmID = 0;
  fe.refill(elmID, rel_order_);
  int nbf = fe.nbf();
  int gpno = nbf;

  fe_accelerate_.resize(gpno, fe);

  for (int g = 0; g < gpno; g++) {
    fe_accelerate_[g].refill(elmID, rel_order_);
  }
  while (fe.next_itg_pt()) {
    const int gpID = fe.cur_itg_pt_num();
    for (int g = 0; g < gpno; g++) {
      if (g >= gpID) {
        fe_accelerate_[g].next_itg_pt();
      }
    }
  }
}

}  // namespace TALYFEMLIB
