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

#include <vector>
#include <string>

#include <talyfem/grid/grid_common.h>  // for ElemType
#include <talyfem/grid/nodeid_types.h>

/**
 * This file defines an interface for reading and writing Tecplot files.
 * An implementation is in TecplotIO_ascii.h/.cpp.
 */

namespace TALYFEMLIB {

/**
 * Represents the value of the "F=" line of a Tecplot header.
 */
enum ElemFormatType {
  kFiniteElementPoint  ///< only supported value (points)
};


/**
 * Used to hold all the data in a Tecplot header.
 */
struct TecplotHeader {
  std::string title;  ///< The title of the graph.
  std::vector<std::string> variables;  ///< All variables stored in the file.

  /**
   * @returns guess at dimensionality of mesh based on variable names
   */
  inline int coord_nsd() const {
    int nsd = 0;
    if (variables.size() >= 1 && (variables[0] == "x" || variables[0] == "X"))
      nsd++;
    if (variables.size() >= 2 && (variables[1] == "y" || variables[1] == "Y"))
      nsd++;
    if (variables.size() >= 3 && (variables[2] == "z" || variables[2] == "Z"))
      nsd++;
    return nsd;
  }
};


/**
 * Used to hold all the data in a Tecplot zone.
 */
struct TecplotZone {
  PetscInt num_nodes;  ///< The number of nodes in the file.
  PetscInt num_elements;  ///< The number of elements in the file.
  ElemFormatType format;  ///< The value of the "F=" line.
  ElemType elem_type;  ///< The type of element in the file.

  /**
   * @returns the number of nodes in each element based on elem_type
   */
  inline unsigned int get_nodes_in_element() const {
    switch (elem_type) {
      case kElem3dHexahedral:
        return 8;
      case kElem2dBox:
        return 4;
      case kElem2dTriangle:
        return 3;
      case kElem3dTetrahedral:
        return 4;
      case kElem1d:
        return 2;
    }

    return 0;
  }
};


/**
 * Used to read data from a Tecplot file.
 *
 * For an example of how to use this,
 * see Tests/TecplotIO/src/tecplot_io_test.cpp.
 */
class TecplotReader {
 public:
  virtual ~TecplotReader() { }

  /**
   * Access the header data after we've read it.
   */
  inline const TecplotHeader& header() const {
    return _header;
  }

  /**
   * Access the current zone's data after we've read it.
   */
  inline const TecplotZone& zone() const {
    return _zone;
  }

  /**
   * Opens the file to read from.
   * @param path path to file to open
   * @throw FileIOException If the file cannot be opened.
   */
  virtual void open(const char* path) = 0;

  /**
   * Reads the header at the top of Tecplot files.
   * That is, the "TITLE" and "VARIABLES" lines in ASCII files.
   * Should be called once, after opening the file.
   * @throw FileIOException if there was a problem reading the header.
   */
  virtual void read_header() = 0;

  /**
   * Reads a "ZONE" line in a Tecplot file.
   * Note that the Tecplot spec says this doesn't have to strictly be on one
   * line, but this implementation requires it.
   * @throw FileIOException if there was a problem reading the ZONE line.
   */
  virtual void read_zone() = 0;

  /**
   * Reads a line of node values (one for every variable in header).
   *
   * Can be called zone().num_nodes times after reading the header.
   *
   * @param node_data_out Array to write the read data to.
   *                      Should be at least header().variables.size() large.
   * @throw FileIOException if there was an error reading the data.
                            (e.g. malformed data, unexpected EOF, etc.)
   */
  virtual void read_node(double* node_data_out) = 0;

  /**
   * Reads a line of element values (based on header's elem_type value).
   *
   * Can be called zone().num_elements times after reading all node data.
   * @param elem_data_out Array to write the read data to. Should be at least
   *                      header().get_nodes_in_element()large.
   * @throw FileIOException if there was an error reading the data.
   */
  virtual void read_elem(PhysicalNodeID* elem_data_out) = 0;

  /**
   * Close the file we're reading from, freeing any open file resources.
   * Does nothing if the file is not open.
   * open() can be called again after this function is called.
   * IMPLEMENTATIONS SHOULD CALL THIS FUNCTION IN THEIR DESTRUCTOR!
   */
  virtual void close() = 0;

 protected:
  TecplotHeader _header;  ///< header information, valid after read_header()
  TecplotZone _zone;  ///< zone information, valid after read_zone()
};


/**
 * Used to write data to a Tecplot file.
 *
 * For an example of how to use this,
 * see Tests/TecplotIO/src/tecplot_io_test.cpp.
 */
class TecplotWriter {
 public:
  virtual ~TecplotWriter() { }

  /**
   * Identical to what was passed to write_header().
   * @returns header data, only valid after calling write_header()
   */
  inline const TecplotHeader& header() const {
    return _header;
  }

  /**
   * Identical to what was passed in to write_zone().
   * @returns zone data, only valid after calling write_zone()
   */
  inline const TecplotZone& zone() const {
    return _zone;
  }

  /**
   * Open the file to write to.
   * @param path path to output file (including extension)
   * @param append whether or not to append to the file
   * @throw FileIOException if there was an error.
   */
  virtual void open(const char* path, bool append) = 0;

  /**
   * Write the data in header() to the file.
   * Should be called once after opening the file.
   * @param header header data to write
   * @throw FileIOException if there was a problem writing data.
   */
  virtual void write_header(const TecplotHeader& header) = 0;

  /**
   * Write the zone data to the file.
   * Should be called after writing the header, but before writing node data.
   * @param zone zone data to write
   * @throw FileIOException if there was a problem writing data.
   */
  virtual void write_zone(const TecplotZone& zone) = 0;

  /**
   * Write node data to a file.
   *
   * Should be called after write_zone().
   * Should be called zone().num_nodes times.
   * @param node_data node data to write, header().variables.size() large
   * @throw FileIOException if there was a problem writing data
   */
  virtual void write_node(const double* node_data) = 0;

  /**
   * Write element connectivity data to a file.
   *
   * Should be called after writing all node data for the current zone.
   * Should be called zone().num_elements times.
   *
   * @param elem_data Element data to write. Should be
   *                  zone().get_nodes_in_element() large.
   * @throw FileIOException if there was a problem writing data
   */
  virtual void write_elem(const PhysicalNodeID* elem_data) = 0;

  /**
   * Close any file resources that may be open.
   * Does nothing if the file is not open.
   * IMPLEMENTATIONS SHOULD CALL THIS FUNCTION IN THEIR DESTRUCTOR!
   */
  virtual void close() = 0;

 protected:
  TecplotHeader _header;  ///< header data, only valid after write_header()
  TecplotZone _zone;  ///< current zone data, only valid after write_zone()
};

}  // namespace TALYFEMLIB
