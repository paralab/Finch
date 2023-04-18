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
#include <talyfem/basis/itg_pts/itg_pts.h>

// constexpr allows us to declare arrays inline, but it does not give them
// external linkage. We use this .cpp file to make the compiler store the
// itg_pts/weights arrays in this .o file.

// If you are trying to add more integration points and the library
// compiles fine but user code gives a lot of "undefined reference to
// ItgPts<x, y, z>::itg_pts" and "ItgPts<x, y, z>::weights", you probably
// forgot to update this file to include your new integration points.

#define DEF_BOX4D_SURF(order, surf_id) \
constexpr decltype(BoxItgPts<order, 4, surf_id>::itg_pts) BoxItgPts<order, 4, surf_id>::itg_pts; \
constexpr decltype(BoxItgPts<order, 4, surf_id>::weights) BoxItgPts<order, 4, surf_id>::weights;

namespace TALYFEMLIB {

// surfaces
#ifdef ENABLE_4D
DEF_BOX4D_SURF(2, -1)
DEF_BOX4D_SURF(2, +1)
DEF_BOX4D_SURF(2, -2)
DEF_BOX4D_SURF(2, +2)
DEF_BOX4D_SURF(2, -3)
DEF_BOX4D_SURF(2, +3)
DEF_BOX4D_SURF(2, -4)
DEF_BOX4D_SURF(2, +4)

DEF_BOX4D_SURF(3, -1)
DEF_BOX4D_SURF(3, +1)
DEF_BOX4D_SURF(3, -2)
DEF_BOX4D_SURF(3, +2)
DEF_BOX4D_SURF(3, -3)
DEF_BOX4D_SURF(3, +3)
DEF_BOX4D_SURF(3, -4)
DEF_BOX4D_SURF(3, +4)

DEF_BOX4D_SURF(4, -1)
DEF_BOX4D_SURF(4, +1)
DEF_BOX4D_SURF(4, -2)
DEF_BOX4D_SURF(4, +2)
DEF_BOX4D_SURF(4, -3)
DEF_BOX4D_SURF(4, +3)
DEF_BOX4D_SURF(4, -4)
DEF_BOX4D_SURF(4, +4)

#endif

}  // namespace TALYFEMLIB
