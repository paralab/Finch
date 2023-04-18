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
#include <talyfem/fem/periodic_data.h>

#include <algorithm>  // for std::remove

#include <talyfem/common/exceptions.h>  // for throwing exception


namespace TALYFEMLIB {

PeriodicData::PeriodicData(PeriodicBounds *periodic_bounds_obj,
                           int n_degrees_of_freedom)
  : periodic_bounds_(periodic_bounds_obj),
    periodic_var_list_(0),
    is_var_periodic_(n_degrees_of_freedom, false),
    is_periodic_(false),
    enable_exchange_(false) { }

void PeriodicData::SetVarIndexPeriodic(int index) {
  // don't allow double setting
  if (is_var_periodic_[index]) {
    throw TALYException() << "Periodic value is already set.";
  }
  is_periodic_ = true;  // now that something is periodic, this is true
  is_var_periodic_[index] = true;  // mark this as periodic
  periodic_var_list_.push_back(index);  // add to list of periodic vars
}

void PeriodicData::SetVarIndexNonPeriodic(int index) {
  // don't allow removing non existing values
  if (!is_var_periodic_[index]) {
    throw TALYException() << "Periodic value is not set.";
  }
  is_var_periodic_[index] = false;  // mark this as nonperiodic
  // remove this value from the list of periodic values
  periodic_var_list_.erase(
    std::remove(periodic_var_list_.begin(), periodic_var_list_.end(), index),
    periodic_var_list_.end());

  // check if there are any remaining periodic variables.
  // If not, this data is no longer periodic
  if (periodic_var_list_.empty()) {
    is_periodic_ = false;
  }
}

}  // namespace TALYFEMLIB
