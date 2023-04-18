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

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petsc.h>
#else
#include <petsc.h>
#endif

#include <stdio.h>
#include <vector>
#include <talyfem/grid/nodeid_types.h>
#include <talyfem/grid/elem.h>

namespace TALYFEMLIB {

/**
 * Each element keeps a list of SurfaceIndicator objects.
 * Each SurfaceIndicator has:
 *  (1) a surface ID
 *  (2) a list of indicators

 * Surface IDs must correspond to ELEM::GetSurfaceCheckArray(),
 * which defines which element-local nodes correspond to which surface ID.
 * (Certain surface IDs are expected by the element basis functions, so it
 * doesn't make much sense for them to be user-defined.)

 * This file format lets you specify an arbitrary list of indicators
 * for the library-defined surfaces of an element.

 * Files look like this:
 * [ number of elements ]
 * [ type of element ] (matches Tecplot "ET=..." value)

 * (line length = number of possible surfaces)
 * (line count = number of elements)
 * [ elm0_surf_0_indicators ] [ elm0_surf_1_indicators ] ...
 * [ elm1_surf_0_indicators ] [ elm1_surf_1_indicators ] ...
 * ...
 * ============================================================================

 * IMPORTANT NOTE:
 * Indicators are stored as a decimal representation of bitflags.
 * That is, the nth bit corresponds to the nth indicator being set.
 * i.e. indicators = 4 in decimal -> ...000100 in binary -> indicator 2 is set
 * (count from the right, starting at 0).

 * You will need to use bitwise operations to convert integer indicator numbers
 * into their bitflag representation when generating a surfaces file.
 * You can do this with:
 * unsigned int indicators = (1 << indicator_num);
 * Or set multiple indicators with a bitwise or:
 * unsigned int indicators = (1 << indicator_1) | (1 << indicator_2);

 * ----------------------------------------------------------------------------

 * Here's an example of how you would make one of these files for a
 * mesh with 2D triangle elements.

 * First, check that element type's GetSurfaceCheckArray() method for
 * the list of surfaces.
 * The surfaces for a 2D triangle are defined in
 * Grid/elem_types/elem2dtriangle.cpp's GetSurfaceCheckArray() method:
 *   [ surface_id, node_id0, node_id1 ]
 *   static const int T3n2DCheckArray[] = {
 *      -1, 0, 2 ,
 *      +1, 1, 2 ,
 *      -2, 0, 1
 *   };

 * which means:

 * Surface | Surface ID | element-local node IDs
 * ---------------------------------------------
 *   0     |     -1     |  0, 2
 *   1     |      1     |  1, 2
 *   2     |     -2     |  0, 1

 * For a 2D triangle mesh with the following global node IDs:
 *  0--2
 *  | /|    ELEMENT 0: (0, 1, 2)
 *  |/ |    ELEMENT 1: (1, 3, 2)
 *  1--3

 * Let's make a surfaces file that outlines the rectangle this mesh makes,
 * using surface indicators that correspond to direction.
 * This is an arbitrary convention for this example. You can use whatever
 * makes sense to you for your problem.
 * 1 --> X-
 * 2 --> X+
 * 4 --> Y-
 * 8 --> Y+
 */

/*
2            # this file has 2 elements
TRIANGLE     # of the 2D triangle type

8 0 1        # element 0 has surface 0 with indicators == 8
             #           and surface 2 with indicators == 1
             # -> (0, 2) = Y+ and (0, 1) = X-

0 2 4        # element 1 has surface 1 with indicators == 2
             #           and surface 2 with indicators == 4
             # -> (3, 2) = X+ and (1, 3) = Y-
*/

/**
 * Represents a header for the surface file format.
 * Contains the number of elements and their type.
 * From the element type, we can infer the number of surfaces per element,
 * number of nodes per surface, and surface IDs of each surface.
 */
class SurfHeader {
 public:
  /**
   * Construct a SurfHeader with no data set.
   */
  SurfHeader() : num_elements_(0), elem_(NULL) {}
  virtual ~SurfHeader();

  /**
   * Set the data in this SurfHeader.
   * @param num_elements number of elements
   * @param elem_type type of element
   */
  void set(PetscInt num_elements, ElemType elem_type);

  /**
   * @returns number of elements in the file after this header
   */
  inline PetscInt num_elements() const {
    return num_elements_;
  }

  /**
   * @returns type of element in this block
   */
  inline ElemType elem_type() const {
    return elem_->elmType();
  }

  /**
   * @returns number of surfaces for this element type
   */
  inline unsigned int num_surfaces() const {
    return elem_->GetSurfaceCount();
  }

  /**
   * @returns nodes per surface for this element type
   */
  inline unsigned int nodes_per_surface() const {
    return elem_->GetNodesPerSurface();
  }

  /**
   * @param surface surface number
   * @returns the surface ID for the given surface
   */
  inline int surface_id(int surface) const {
    return elem_->GetSurfaceCheckArray()[surface * (nodes_per_surface() + 1)];
  }

 private:
  PetscInt num_elements_;  ///< number of elements after this header
  ELEM* elem_;  ///< fake element that we use to infer surface data
};


/**
 * Reads a surfaces file.
 */
class SurfReader {
 public:
  /**
   * Construct a SurfReader.
   */
  SurfReader();
  virtual ~SurfReader();

  /**
   * Start reading a surfaces file.
   * @param file file to read from
   */
  void open(const char* file);

  /**
   * Close an opened surfaces file.
   * Only valid after calling open().
   * Automatically called on destruction.
   */
  void close();

  /**
   * Read the header from the file.
   * Should be the first thing you do after opening a file.
   */
  void read_header();

  /**
   * Read the next element's surface data from the file.
   * Writes the surface data for ALL POSSIBLE SURFACES
   * into out. So, out should be at least num_surfaces() large.
   * @param[out] out surface indicator output array
                     should be at least num_surfaces() large
   */
  void read_elm(SurfaceIndicator::IndicatorType* out);

  /**
   * Get the number of elements in this file.
   * Only valid after calling read_header().
   * @returns the number of elements in this file
   */
  inline PetscInt num_elements() const {
    return header_.num_elements();
  }

  /**
   * Get the number of surfaces for each element in this file.
   * Use this to know how big the read_elm(...) output argument should be.
   * Only valid after calling read_header().
   * @returns number of surfaces
   */
  inline unsigned int num_surfaces() const {
    return header_.num_surfaces();
  }

  /**
   * Get the number of nodes per surface for each element in this file.
   * Only valid after calling read_header().
   * @returns number of nodes per surface
   */
  inline unsigned int nodes_per_surface() const {
    return header_.nodes_per_surface();
  }

  /**
   * Used to map the indices of read_elm to surface IDs.
   * @param surface idx of the read_elm output array to get the surface ID for
   * @returns surface ID for the surfaceth entry in the read_elm() output array
   */
  inline int surface_id(int surface) const {
    return header_.surface_id(surface);
  }

 private:
  FILE* file_;  ///< file handle for the currently open()'d file
  SurfHeader header_;  ///< header, set by read_header()
};

}  // namespace TALYFEMLIB
