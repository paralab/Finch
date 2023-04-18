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
#ifndef GRID_SURFACEINDICATOR_H_
#define GRID_SURFACEINDICATOR_H_

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

// The capacity of IndicatorType determines how many surface indicators
// are allowed per surface. Using a 32-bit datatype gives you 32 possible
// indicators, and using a 64-bit datatype gives you 64 possible indicators.

// How to switch between 32-bit and 64-bit:
// 1. Change SURFACE_INDICATOR_FORMAT to the appropriate printf format string.
//    You should probably use the PRIu[32/64] macros for this.
//    (they are defined in inttypes.h by the standard, and are portable).
//    This is used by FileIO/surface_io.cpp to read surface files.
// 2. Change MPI_SURFACE_INDICATOR to the MPI_Datatype for transmitting
//    a surface indicator.
//    MPI_UNSIGNED for 32-bit, MPI_LONG_LONG_INT for 64-bit.
//    This is used in CMeshPartition::LoadElmSurfaceData for sending indicators.
// 3. Change the SurfaceIndicator::IndicatorType typedef.
//    uint32_t for 32-bit or uint64_t for 64-bit.
//    (Do NOT use unsigned long long, as it is not guaranteed to be 64 bits
//     on all architectures.)

#define SURFACE_INDICATOR_FORMAT "%" PRIu64
#define MPI_SURFACE_INDICATOR MPI_UINT64_T

#include <talyfem/grid/zeroptv.h>

namespace TALYFEMLIB {

/**
 * Holds information about a surface for an element.
 */
class SurfaceIndicator {
 public:
  typedef uint64_t IndicatorType;  ///< indicator storage (stored as bitflags)

  /**
   * Maximum surface indicator number possible.
   * Should match IndicatorType.
   */
  static const unsigned int MAX_SURFACE_INDICATORS = sizeof(IndicatorType) * 8;

  /**
   * Create a new SurfaceIndicator.
   * @param surfaceID surface ID
   */
  explicit SurfaceIndicator(int surfaceID = -1);

  /**
   * @returns surface ID for this surface
   */
  inline int surface_id() const {
    return surface_id_;
  }

  /**
   * @returns number of indicators set
   */
  int num_indicators() const;

  /**
   * @param id indicator ID to check
   * @returns if indicator id is set
   */
  bool has_indicator(int id) const;

  /**
   * Set indicator id as "set."
   * @param id indicator ID to set as true
   */
  void add_indicator(int id);

  /**
   * Set the raw indicator flags.
   * You probably don't want to do this unless you're writing a loader.
   * @param indicators raw indicator flags to set
   */
  void set_indicators(IndicatorType indicators);

  /**
   * Set the normal for this surface.
   * @param normal_vector normal to save for this surface
   */
  inline void set_normal(const ZEROPTV& normal_vector) {
    normal_ = normal_vector;
  }

  /**
   * @returns the precalculated normal for this surface (set with set_normal())
   */
  inline const ZEROPTV& normal() const { return normal_; }

 private:
  int surface_id_;  ///< surface ID
  IndicatorType indicators_;  ///< surface indicators (bitflags)
  ZEROPTV normal_;  ///< cached normal calculation for this surface
};

}  // namespace TALYFEMLIB

#endif  // GRID_SURFACEINDICATOR_H_
