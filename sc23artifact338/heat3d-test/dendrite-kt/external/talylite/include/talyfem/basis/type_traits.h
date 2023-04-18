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
#pragma once

#include <cstdint>
#include <type_traits>
#include <talyfem/utils/macros.h>  // for DEFINE_HAS_SIGNATURE

namespace TALYFEMLIB {

class ZEROPTV;

/// Check if a type has a calc_d2Nde function.
DEFINE_HAS_SIGNATURE(has_calc_d2Nde, T::calc_d2Nde, void (*)(const ZEROPTV&,
                     double (&)[T::nbf][T::nsd * (T::nsd + 1) / 2]));

/// Check if a type has a calc_dXde function.
DEFINE_HAS_SIGNATURE(has_calc_dXde, T::calc_dXde, void (*)(const ZEROPTV&,
                     const double (&)[T::nbf][T::nsd],
                     const ElemNodes&, double (&)[T::nsd][T::nsd]));

/// Check if a type has a calc_position function.
DEFINE_HAS_SIGNATURE(has_calc_position, T::calc_position,
                     void (*)(const ZEROPTV&, const double (&)[T::nbf],
                     const ElemNodes&, ZEROPTV*));

}  // namespace TALYFEMLIB
