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

#ifdef ENABLE_HDF5

#include <hdf5.h>

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petsc.h>
#else
#include <petsc.h>
#endif

#include <vector>
#include <string>

#include <talyfem/common/exceptions.h>
#include <talyfem/grid/grid_common.h>  // for ElemType

namespace TALYFEMLIB {

// Our HDF5 file is laid out like this:
// GROUP [output_id]  ("root group")
//   DATASET node_coordinates
//   GROUP node_data
//     DATASET node_attrib0
//     DATASET node_attrib1
//     ...
//     DATASET node_attribN-1
//   DATASET elem_connectivity

/**
 * Simple struct that holds an attribute's name and opened H5 dataset.
 */
struct H5Attrib {
  std::string name;
  hid_t dset;

  H5Attrib() : dset(0) {}
  H5Attrib(const std::string& n, hid_t ds) : name(n), dset(ds) {}
};

/**
 * This class is used to open or create a "root group" in a file.
 * It is intended for use primarily by the H5Writer/H5Reader classes.
 *
 * Actual HDF5 reading/writing is not performed by this class; this class
 * is only a wrapper for opening/creating a group, providing convenient
 * getters for the various HDF5 objects inside the root group.
 * YOU PROBABLY WANT TO USE H5READER/H5WRITER INSTEAD, which completely wraps
 * the HDF5 calls.
 *
 * It can either open an existing group in a file (with open()) OR
 * create a new group (with create_group()).
 *
 * Once a group has been opened or created, the relevant datasets in the file
 * can be accessed with various getters (root(), coord_dset(), elem_dset(),
 * attribs()[i]). HDF5 functions can be used to read/write to these datasets.
 */
class H5OutputGroup {
 public:
  H5OutputGroup();
  virtual ~H5OutputGroup();

  /**
   * Opens a group in an open HDF5 file.
   * This H5OutputGroup does NOT take responsibility for closing file_id.
   * @param file_id Open HDF5 file descriptor to use.
   * @param group_name Group to read. If NULL or empty string, will open the
   *                   first available group in the file.
   *                   Fails if no group exists in the file.
   */
  void open(hid_t file_id, const char* group_name);

  /**
   * Creates a new group in the HDF5 file opened with file_id.
   * This H5OutputGroup does NOT take responsibility for closing file_id.
   * WILL FAIL IF group_name ALREADY EXISTS IN TEH ROOT OF THE FILE.
   * @param file_id Open HDF5 file descriptor to use.
   * @param group_name Name of the group to create.
   * @param grid_type Type of the grid for the group we're creating (static).
   * @param attribute_names Vector of attribute names present in this group.
   * @param total_node_count Total number of nodes in this group.
   * @param total_elem_count Total number of elements in this group.
   */
  void create_group(hid_t file_id, const char* group_name,
                  int nsd, ElemType elem_type,
                  const std::vector<std::string>& attribute_names,
                  unsigned int total_node_count,
                  unsigned int total_elem_count);

  /**
   * Closes HDF5 resources this group has opened via open()/create_group().
   * WILL NOT CLOSE file_ids passed in from open()/create_group().
   * Called automatically by destructor.
   */
  void close();

  // simple getters
  /** HDF5 group for the "root group." */
  inline hid_t root() { return _root; }
  /** HDF5 dataset for node coordinates. */
  inline hid_t coord_dset() { return _coord_dset; }
  /** HDF5 dataset for element connectivity. */
  inline hid_t elem_dset() { return _elem_dset; }
  /**
   * Vector of HDF5 attributes present in the "root group."
   * Access with attribs()[i].dset.
   */
  inline const std::vector<H5Attrib>& attribs() { return _attribs; }

  /** Type of grid. NSD is automatically calculated from this. */
  inline ElemType elem_type() const { return _elem_type; }
  inline unsigned int nsd() const { return _nsd; }
  inline unsigned int nodes_per_elem() const { return _nodes_per_elem; }

  // getters based on hdf5 file data
  // (these must be read from file every time, so the efficiency is not good)
  PetscInt node_count() const;  ///< Total number of nodes.
  PetscInt elem_count() const;  ///< Total number of elements.

 private:
  // helper function
  void create_datasets(unsigned int nsd,
    const std::vector<std::string>& attrib_names, unsigned int total_nodes,
    unsigned int nodes_per_elem, unsigned int total_elems);

  hid_t _root;
  hid_t _coord_dset;
  hid_t _elem_dset;

  // guaranteed to be in the same order as the initial attribute_names() vector
  // supplied in create_group, even across file/group re-opens
  std::vector<H5Attrib> _attribs;

  // these values also exist on the file, but are accessed frequently when
  // writing (and are determined from grid_type, so may require conversion)
  // so, we cache them here to avoid file access
  ElemType _elem_type;
  unsigned int _nsd;
  unsigned int _nodes_per_elem;
};


/**
 * This class is used to write data in the HDF5 format.
 * See H5IO_gf.h for example usage.
 */
class H5Writer {
 public:
  explicit H5Writer(bool parallel);
  virtual ~H5Writer();

  /**
   * Open a file to write to.
   * @param filename Path to the file to write to (include the .h5 extension).
   * @param append Whether or not to append to existing files.
   */
  void open_file(const char* filename, bool append);

  /**
   * Create a new "root group" in the file.
   * SHOULD BE CALLED ON ALL PROCESSES PERFORMING THE WRITE,
   * WITH IDENTICAL PARAMETERS!
   * @param group_name Name of the group we're creating.
   * @param grid_type Type of the grid. NSD is calculated from this.
   * @param attribute_names Vector of attributes we're saving.
   * @param total_node_count Total number of nodes.
   * @param total_elem_count Total number of elements.
   */
  void create_group(const char* group_name, int nsd, ElemType elem_type,
                    const std::vector<std::string>& attribute_names,
                    unsigned int total_node_count,
                    unsigned int total_elem_count);


  // Set collective to true if you plan on doing one big write,
  // instead of many smaller writes.
  // If each process does not perform the same number of writes when using
  // collective mode, the program will freeze!

  /**
   * Write a block of contiguous node data to the destination file.
   * @param node_offset Index into the file where this block will be written.
   * @param node_count Number of nodes in this block.
   * @param node_coords Array of node coordinates. Should be
                        (nsd() * node_count) doubles large.
   * @param node_attribs Array of attributes to write. Should be
   *                     (attribs().size() * node_count) doubles large.
   * @param collective If true, HDF5 will try to use MPI-IO to orchestrate
                       writes across every process in an optimal order.
                       May improve write throughput.
   */
  void write_node_data(PetscInt node_offset, PetscInt node_count,
                       const double* node_coords, const double* node_attribs,
                       bool collective = false);

  /**
   * Write a block of element data.
   * @param elem_offset Index into the file where this block will be written.
   * @param elem_count Number of elements in this block.
   * @param elem_data Array of element connectivity data. Should be
                      (nodes_per_elem * elem_count) PetscInts large.
   * @param collective If true, HDF5 will try to use MPI-IO to orchestrate
                       writes across every process in an optimal order.
                       May improve write throughput.
   */
  void write_elem_data(PetscInt elem_offset, PetscInt elem_count,
                       const PetscInt* elem_data, bool collective = false);

  /**
   * Used to write a block of elements that is not
   * necessarily contiguous in the destination file.
   * @param elem_offsets Array of element IDs. Determines where
   *                     in the file a particular element will be written.
   * @param elem_count Number of elements to write.
   * @param elem_data Array of element connectivity data. Should be
   *                  (nodes_per_elem * elem_count) PetscInts large.
   * @param collective If true, HDF5 will try to use MPI-IO to orchestrate
   *                   writes across every process in an optimal order.
   *                   May improve write throughput.
   */
  void write_elem_data_multi(PetscInt* elem_offsets, PetscInt elem_count,
                             const PetscInt* elem_data,
                             bool collective = false);

  /**
   * Used to write a block of nodes that is not
   * necessarily contiguous in the destination file.
   * This is useful for setting up large collective writes, which gives the
   * best throughput. Recommended over write_node_data() on parallel systems.
   * @param node_offsets Array of node IDs. Determines where
                         in the file a particular node will be written.
   * @param node_count Number of nodes to write.
   * @param node_coords Coordinates of each node.
   * @param node_attribs Array of arrays of coordinates to write.
   * @param collective If true, HDF5 will try to use MPI-IO to orchestrate
                       writes across every process in an optimal order.
                       May improve write throughput.
   */
  void write_node_data_multi(PetscInt* node_offsets, PetscInt node_count,
                       const double* node_coords, const double** node_attribs,
                       bool collective = false);

  /**
   * Generate an XDMF file for ParaView from all data in the open HDF5 file.
   * Should only be called after open_file().
   * Can be called without create_group(), since it uses all data in the file.
   * You should only call this on one process.
   */
  void write_xdmf();

  /**
   * Close an open file. Called automatically by destructor.
   */
  void close();

 private:
  bool _parallel;
  std::string _file_path;
  H5OutputGroup _out_group;

  hid_t _file_id;
};

/**
 * This class is used to read data from the HDF5 format.
 * See H5IO_gf.h for example usage.
 * For parallel systems, multiple processes can keep the same file open
 * and read from different areas simultaneously.
 */
class H5Reader {
 public:
  explicit H5Reader(bool parallel);
  virtual ~H5Reader();

  void open_file(const char* filename);
  void open_group(const char* group_name);

  inline ElemType elem_type() const { return _in_group.elem_type(); }
  inline unsigned int nsd() const { return _in_group.nsd(); }
  inline unsigned int nodes_per_elem() const {
    return _in_group.nodes_per_elem();
  }
  inline const std::vector<H5Attrib>& attribs() { return _in_group.attribs(); }
  inline PetscInt node_count() const { return _in_group.node_count(); }
  inline PetscInt elem_count() const { return _in_group.elem_count(); }

  /**
   * Read a block of node data (coordinates + attributes).
   * Will read nsd() doubles per node_count into node_coords, and
   * attribs.size() doubles per node_count into node_attribs.
   * @param node_offset Node ID to start reading from (zero-based).
   * @param node_count Number of nodes to read.
                       (will read from node_offset to node_offset + node_count)
   * @param node_coords Where to store the node coordinates. Should be at least
                        (nsd() * node_count) doubles large.
   * @param node_attribs Where to store the node attributes. Should be at least
                         (attribs.size() * node_count) doubles large.
   */
  void read_node_data(PetscInt node_offset, PetscInt node_count,
                       double* node_coords, double* node_attribs);
  /**
   * Read a block of element connectivity data from the file.
   * @param elem_offset Element ID to start reading from (zero-based).
   * @param elem_count Number of elements to read.
                       (will read from elem_offset to elem_offset + elem_count)
   * @param elem_data Where to store the element conenctivity data.
                      Should be at least (nodes_per_elem() * elem_count) large.
   */
  void read_elem_data(PetscInt elem_offset, PetscInt elem_count,
                       PetscInt* elem_data);

  /**
   * Close the open file.
   * Automatically called by the destructor.
   */
  void close();

 private:
  bool _parallel;
  H5OutputGroup _in_group;
  hid_t _file_id;
};

}  // namespace TALYFEMLIB

#endif

