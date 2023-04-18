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
#ifndef GRID_NODEDATA_H_
#define GRID_NODEDATA_H_

#include <cstdio>  // for sprintf


namespace TALYFEMLIB {

/**
 * Stores the data for a node in the gridfield
 *
 * This is intended as a base class for node data in user code
 */
class NODEData {
 public:
  /**
   * Returns value of given value for the node.
   *
   * @param index index of data item to obtain
   * @return reference to desired data item
   */
  virtual double& value(int index) {
    return *(reinterpret_cast<double*>(this));
  }

  /**
   * Returns value of given value for the node.
   *
   * @param index index of data item to obtain
   * @return reference to desired data item
   */
  virtual const double& value(int index) const {
    return *(reinterpret_cast<const double*>(this));
  }

  /**
   * Returns name of data stored in the given index
   *
   * @param index index of data to return name of
   * @return name of data at given index
   */
  static const char* name(int index) {
    static char str[256];
    snprintf(str, sizeof(str), "var%d", index);
    return str;
  }

  /**
   * Returns number of data items stored in this object
   *
   * @return number of data items in this object
   */
  static int valueno() {
    return 0;
  }

  /**
   * update data structures
   * TODO: this is vaguely defined... when is this used and why?
   * TODO: is there a point to having this in the base class?
   */
  virtual void UpdateDataStructures() { }

  /**
   * TODO: is there a point to having this in the base class?
   */
  virtual int copy_node_data(char* data, bool bO2D) {
    return 0;
  }

  /**
   * Save the node data
   * TODO: saving is done by the fileIO routines which use a different format
   * TODO: for each file type
   * TODO: is there a point to having this in the base class?
   */
  virtual void save(FILE* fp) { }

  /**
   * Load the node data
   * TODO: loading is done by the fileIO routines which use a different format
   * TODO: for each file type
   * TODO: is there a point to having this in the base class?
   */
  virtual void load(FILE* fp) { }

  virtual ~NODEData() { }
};

}  // namespace TALYFEMLIB

#endif  // GRID_NODEDATA_H_
