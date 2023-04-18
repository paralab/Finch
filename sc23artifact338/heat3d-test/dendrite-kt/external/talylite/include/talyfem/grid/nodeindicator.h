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
#ifndef GRID_NODEINDICATOR_H_
#define GRID_NODEINDICATOR_H_

#include <stdint.h>  // for uint32_t


namespace TALYFEMLIB {

/**
 * Used to store a set of indicators in a bitmask. The nth bit of a
 * NodeIndicator contains the state of indicator number n.
 */
typedef uint32_t NodeIndicator;

#define MAX_NODE_INDICATORS sizeof(NodeIndicator)*8

// Format string for the printf functions.
// If NodeIndicator is ever switched to 64-bit, change this to %ld.
#define NODE_INDICATOR_FORMAT "%d"

/** 
 * Get a boundary indicator's bitmask representation by number.
 * TODO - this should be marked constexpr once we switch to C++11
 * @param id The number of the boundary indicator to get. 
 * @returns the node indicator mask for id
 */
inline NodeIndicator INDICATOR_NUM(uint32_t id) {
  // make sure we're not asking for a higher node indicator than we can store
  // with this setting
  // assert(id >= 0 && id < MAX_NODE_INDICATORS);
  return 1 << id;
}

}  // namespace TALYFEMLIB

#endif  // GRID_NODEINDICATOR_H_
