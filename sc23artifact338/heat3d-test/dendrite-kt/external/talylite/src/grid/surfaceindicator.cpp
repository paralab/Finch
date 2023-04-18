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
#include <talyfem/grid/surfaceindicator.h>

#include <talyfem/utils/utils.h>  // for popcount


namespace TALYFEMLIB {

SurfaceIndicator::SurfaceIndicator(int surfaceID) {
  surface_id_ = surfaceID;
  indicators_ = 0;
  normal_ = ZEROPTV(0, 0, 0);
}

int SurfaceIndicator::num_indicators() const {
  return popcount<IndicatorType>(indicators_);
}

bool SurfaceIndicator::has_indicator(int id) const {
  return (indicators_ & ((IndicatorType)1 << id)) != 0;
}

void SurfaceIndicator::add_indicator(int id) {
  indicators_ = indicators_ | ((IndicatorType)1 << id);
}

void SurfaceIndicator::set_indicators(IndicatorType indicators) {
  indicators_ = indicators;
}

}  // namespace TALYFEMLIB
