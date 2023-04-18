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

#define INSTANTIATE_BOX(order, nsd, surf_id) \
constexpr decltype(BoxItgPts<order, nsd, surf_id>::itg_pts) BoxItgPts<order, nsd, surf_id>::itg_pts; \
constexpr decltype(BoxItgPts<order, nsd, surf_id>::weights) BoxItgPts<order, nsd, surf_id>::weights;

#define INSTANTIATE_TRI(order, surf_id) \
constexpr decltype(TriItgPts<order, surf_id>::itg_pts) TriItgPts<order, surf_id>::itg_pts; \
constexpr decltype(TriItgPts<order, surf_id>::weights) TriItgPts<order, surf_id>::weights;

#define INSTANTIATE_TET(order, surf_id) \
constexpr decltype(TetItgPts<order, surf_id>::itg_pts) TetItgPts<order, surf_id>::itg_pts; \
constexpr decltype(TetItgPts<order, surf_id>::weights) TetItgPts<order, surf_id>::weights;

namespace TALYFEMLIB {

// volumes
    INSTANTIATE_BOX(1, 1, 0)
    INSTANTIATE_BOX(2, 1, 0)
    INSTANTIATE_BOX(3, 1, 0)
    INSTANTIATE_BOX(4, 1, 0)
    INSTANTIATE_BOX(5, 1, 0)
    INSTANTIATE_BOX(6, 1, 0)
    INSTANTIATE_BOX(7, 1, 0)
    INSTANTIATE_BOX(8, 1, 0)
    INSTANTIATE_BOX(9, 1, 0)
#ifdef CLENSHAW_CURTIS
    INSTANTIATE_BOX(10, 1, 0)
#endif
//
    INSTANTIATE_TRI(2, 0)
    INSTANTIATE_TRI(3, 0)

    INSTANTIATE_TET(2, 0)
    INSTANTIATE_TET(3, 0)
    INSTANTIATE_TET(4, 0)

}  // namespace TALYFEMLIB
