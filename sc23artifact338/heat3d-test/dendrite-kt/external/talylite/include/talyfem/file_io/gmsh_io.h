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

#include <stdint.h>

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petsc.h>
#else
#include <petsc.h>  // for PetscInt
#endif

#include <string>
#include <vector>
#include <map>
#include <fstream>

#include <talyfem/grid/elem_common.h>  // for TALYFEMLIB::ElemType
#include <talyfem/grid/surfaceindicator.h>  // for TALYFEM::SurfaceIndicator

namespace gmsh {

typedef PetscInt NodeID;
typedef uint32_t ElmID;
typedef int32_t Tag;

/**
 * Gmsh element type.
 * Mostly taken from the Gmsh documentation, a couple higher-order 2D types
 * were found by experimentation.
 */
enum ElementType : int {
  NONE = 0,

  LINE_ORDER1 = 1,
  TRIANGLE_ORDER1 = 2,
  BOX_ORDER1 = 3,
  TETRAHEDRON_ORDER1 = 4,
  HEXAHEDRON_ORDER1 = 5,
  PRISM_ORDER1 = 6,
  PYRAMID_ORDER1 = 7,
  LINE_ORDER2 = 8,
  TRIANGLE_ORDER2 = 9,
  BOX_ORDER2 = 10,
  TETRAHEDRON_ORDER2 = 11,
  HEXAHEDRON_ORDER2 = 12,
  PRISM_ORDER2 = 13,
  PYRAMID_ORDER2 = 14,
  POINT = 15,
  BOX_ORDER2_NOCENTER = 16,
  TRIANGLE_9NODE = 20,  // from experimentation
  TRIANGLE_10NODE = 21,  // from experimentation
  LINE_ORDER3 = 26,
  LINE_ORDER4 = 27,
  TETRAHEDRON_ORDER3 = 29,
  TETRAHEDRON_ORDER4 = 30,
  TETRAHEDRON_ORDER5 = 31,
  BOX_ORDER3 = 36,  // found from what gmsh exported - NOT confirmed in docs!
  HEXAHEDRON_ORDER3 = 92,
  HEXAHEDRON_ORDER4 = 93
};

/**
 * Convert a Gmsh element type to a TalyFEM element type.
 * Will implicitly lose order information.
 * @param t element type
 * @returns TalyFEM element type corresponding to t
 */
TALYFEMLIB::ElemType gmsh_elm_to_taly_elm(gmsh::ElementType t);

/**
 * Get number of nodes in an element of the given type.
 * Takes element type's order into account.
 * @param type type of element
 * @returns number of nodes in the given element type
 */
int nodes_per_element(gmsh::ElementType type);

/**
 * Get the order of an element type.
 * @param type element type to get order of
 * @returns order of type
 */
int get_order(gmsh::ElementType type);

/**
 * Converts an element type to its order-1 version.
 * Ex: HEXAHEDRON_ORDER3 -> HEXAHEDRON_ORDER1.
 * Returns type of type is already order 1.
 * @param type element type to convert
 * @returns order 1 version of type
 */
gmsh::ElementType to_order_1(gmsh::ElementType type);

/**
 * Returns the Gmsh element type that "surfaces" for the given element
 * type should be (e.g. line for box, box for hex).
 * @param type type of element to get surface for
 * @returns element type of surface for type
 */
ElementType surface_type_of_element(ElementType type);

/**
 * Reorder Gmsh connectivity to match TalyFEM's basis shape function order.
 * Called automatically when reading an element.
 * @param type element type
 * @param[in,out] connectivity connectivity list to re-order (in-place)
 */
void reorder_connectivity_for_taly(gmsh::ElementType type,
                                   std::vector<gmsh::NodeID>& connectivity);

/**
 * Holds information about a $Section in a Gmsh file.
 */
struct Section {
  std::string name;  ///< name of the section (ex "$Nodes" would be "Nodes")
  std::streamoff start;  ///< byte offset into the file for the start
      ///< of the section - note that for Nodes and Elements, this is the start
      ///< of the *data*, past the "number of nodes"/"number of elms" line
  std::streamoff length;  ///< total number of bytes this section is from start
      ///< (for Nodes/Elements, this does not include the "number of X" line)
};

/**
 * Represents a Gmsh node.
 * Contains all the information found in a node line in a Gmsh file.
 */
struct Node {
  NodeID number;  ///< 1-indexed node ID
  double coords[3];  ///< X/Y/Z coordinates
};

/**
 * Represents a Gmsh element.
 * Contains all the information found in an element line in a Gmsh file.
 */
struct Element {
  ElmID number;  ///< 1-indexed Gmsh element number
  ElementType type;  ///< element type ID (cast to enum)
  std::vector<Tag> tags;  ///< list of all element tags
      /// tags[0] is usually "the physical entity to which the element belongs",
      /// tags[1] is usually "the number of the elementary geometrical entity
      /// to which the element belongs" according to the Gmsh docs
  std::vector<NodeID> connectivity;  ///< 1-indexed list of node IDs
};


/**
 * Convert a Gmsh element's tags field to a TalyFEM SurfaceIndicator bitmask.
 * @param elem element to use
 * @returns TalyFEM surface indicators corresponding to elem.tags
 */
TALYFEMLIB::SurfaceIndicator::IndicatorType tags_to_surface_indicators(
    const gmsh::Element& elem);

/**
 * Reads a Gmsh file.
 * Currently only supports ASCII files.
 */
class Reader {
 public:
  Reader();
  virtual ~Reader();

  /**
   * Open a Gmsh file.
   * @param path path to the Gmsh file (typically *.msh)
   */
  void open(const std::string& path);

  /**
   * Close the open Gmsh file.
   * Does nothing if no file is open.
   */
  void close();

  /**
   * The number of nodes in the Nodes section.
   * This value is cached when the file is first opened.
   * @returns the number of nodes
   */
  inline NodeID n_nodes() const { return n_nodes_; }

  /**
   * Read a node from the file.
   * For the first call, the file will automatically skip to start of the
   * Nodes section. Subsequent read_node() calls will continue from where
   * the last read_node() call finished.
   * Calling read_node() more than n_nodes() times consecutively is undefined.
   * @returns the next Node line in the file
   */
  Node read_node();

  /**
   * The number of elements in the Elements section.
   * This value is cached when the file is first opened.
   * @returns the number of elements
   */
  inline ElmID n_elements() const { return n_elements_; }

  /**
   * Read an element from the file.
   * For the first call, the file will automatically skip to the start of
   * the Elements section. Subsequent read_element() calls will continue from
   * where the lsat read_elements() call finished.
   * Calling read_element() more than n_elements() times consecutively is
   * undefined.
   * NOTE: Element connectivity ordering may not exactly what is in the file.
   *       Connectivity is passed through reorder_connectivity_for_taly(...)
   *       before being returned.
   * @returns the next element in the file
   */
  Element read_element();

  /**
   * Get a list of all sections in the file (including non-standard sections).
   * @returns list of all sections
   */
  inline const std::vector<Section>& sections() const { return sections_; }

  /**
   * Get a map that contains counts for all elements in the file.
   * @returns map of element counts
   */
  inline const std::map<ElementType, NodeID>& elm_type_counts() const {
    return elm_type_counts_;
  }

  /**
   * Get the "primary element type" for the file.
   * The "primary element type" is defined as the element type with the highest
   * number of elements in the file. This is not Gmsh standard.
   * @returns the "primary element type" in this file
   */
  inline ElementType primary_elm_type() const { return primary_elm_type_; }

  /**
   * The "primary surface type" is surface_type_of_element(primary_elm_type).
   * Ex: TRIANGLE_ORDER1 -> LINE_ORDER1, BOX_ORDER2 -> LINE_ORDER2,
   *     HEXAHEDRON_ORDER1 -> BOX_ORDER1.
   * This is not Gmsh standard.
   * @returns the "primary surface type" for the file
   */
  inline ElementType primary_surf_type() const { return primary_surf_type_; }

  /**
   * @returns the nsd of primary_elm_type()
   */
  inline int nsd() const { return nsd_; }

 protected:
  /**
   * Finds all $Sections in the Gmsh file and populates sections_.
   */
  void build_section_list();

  /**
   * Validates sections_ and caches some values.
   */
  void check_section_list();

  /**
   * Moves file_'s cursor to sec.
   * @param sec section to move to
   */
  void seek_section(Section* sec);

  std::ifstream file_;  ///< the open Gmsh file

  // this isn't a map because the Gmsh manual explicitly allows
  // multiple sections with the same name
  std::vector<Section> sections_;  ///< list of all sections
  Section* cur_sec_;  ///< current section
  Section* nodes_sec_;  ///< $Nodes section info (not including count line)
  Section* elements_sec_;  ///< $Elements section info (not including count)

  // set in build_section_list()
  ///< number of each type of element found in the file
  std::map<ElementType, NodeID> elm_type_counts_;

  // cached values, set in check_section_list()
  NodeID n_nodes_;  ///< number of nodes in the Nodes section
  ElmID n_elements_;  ///< number of elements in the Elements section
  int nsd_;  ///< the nsd of primary_elm_type_
  ElementType primary_elm_type_;  ///< most common element in the file
  ElementType primary_surf_type_;  ///< surface_type_of_element(primary_elm)
};

}  // namespace gmsh

