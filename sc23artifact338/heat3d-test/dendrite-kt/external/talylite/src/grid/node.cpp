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
#include <talyfem/grid/node.h>

#include <talyfem/utils/utils.h>  // for popcount


namespace TALYFEMLIB {

NODE::NODE() {
  indicators_ = 0;
  elem_id_ = 0;
}

void NODE::setCoor(const double x_val, const double y_val , const double z_val , const double t_val){
  this->x() = x_val;
  this->y() = y_val;
  this->z() = z_val;
#ifdef ENABLE_4D
  this->t() = t_val;
#endif
}

int NODE::getIndicatorNo() const {
  return popcount<NodeIndicator>(indicators_);
}

}  // namespace TALYFEMLIB
