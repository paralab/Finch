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
#ifdef ENABLE_HDF5
#include <talyfem/file_io/h5_io.h>

#include <assert.h>
#include <string>
#include <vector>

#include <talyfem/grid/grid_common.h>  // for GridType and ElemType


namespace TALYFEMLIB {

// format we store data in
// you could potentially change these to be smaller to save space -
// hdf5 should automatically perform type conversion during writing
// (part of the write call takes an "input data format" parameter)
#define NODE_COORDINATE_TYPE H5T_NATIVE_DOUBLE
#define NODE_DATA_TYPE H5T_NATIVE_DOUBLE

#ifdef PETSC_USE_64BIT_INDICES
  #define ELEM_CONNECTIVITY_INPUT_TYPE H5T_NATIVE_LONG
  #define ELEM_CONNECTIVITY_DATA_TYPE H5T_NATIVE_LONG
#else
  #define ELEM_CONNECTIVITY_INPUT_TYPE H5T_NATIVE_INT
  #define ELEM_CONNECTIVITY_DATA_TYPE H5T_NATIVE_INT
#endif

#define UNUSED(x) (void)(x)

// used to get all the names of subgroups in a group
typedef std::vector<std::string> group_vec_t;
herr_t list_h5_groups(hid_t group, const char* name, void* vec_raw) {
  group_vec_t* vector = reinterpret_cast<group_vec_t*>(vec_raw);
  vector->push_back(name);
  return 0;
}


// H5OutputGroup
H5OutputGroup::H5OutputGroup() {
  _root = -1;
  _coord_dset = -1;
  _elem_dset = -1;
}

H5OutputGroup::~H5OutputGroup() {
  close();
}

void H5OutputGroup::open(hid_t file_id, const char* group_name) {
  close();

  // if no group name is supplied, pick the first group available in the file
  group_vec_t group_names;
  if (group_name == NULL || group_name[0] == '\0') {
    H5Giterate(file_id, ".", 0, &list_h5_groups, &group_names);
    if (group_names.size() == 0) {
      throw FileIOException()
        << "No groups present in file, no group name specified to load!";
    }
    // use the first group found
    group_name = group_names.at(0).c_str();
  }

  // open the group
  _root = H5Gopen(file_id, group_name, H5P_DEFAULT);

  // open the coordinate and element datasets
  _coord_dset = H5Dopen(_root, "node_coords", H5P_DEFAULT);
  _elem_dset = H5Dopen(_root, "elem_connectivity", H5P_DEFAULT);

  // open the attribute datasets
  // get a list of all the groups in the "node_data" group
  // TODO: verify order on this is correct (should be in order of creation)
  group_vec_t data_groups;
  H5Giterate(_root, "node_data", 0, &list_h5_groups, &data_groups);

  std::string attrib_path_base = "node_data/";
  for (unsigned int i = 0; i < data_groups.size(); i++) {
    std::string path = attrib_path_base + data_groups.at(i);
    _attribs.push_back(H5Attrib(
        data_groups.at(i),
        H5Dopen(_root, path.c_str(), H5P_DEFAULT)));
  }

  // get _elem_type from the "grid_type" space
  hid_t nsd_attr = H5Aopen(_root, "nsd", H5P_DEFAULT);
  H5Aread(nsd_attr, H5T_NATIVE_INT, &_nsd);
  H5Aclose(nsd_attr);

  hid_t elem_type_attr = H5Aopen(_root, "elem_type", H5P_DEFAULT);
  H5Aread(elem_type_attr, H5T_NATIVE_INT, &_elem_type);
  H5Aclose(elem_type_attr);

  _nodes_per_elem = get_nodes_in_element(_elem_type);
}

void H5OutputGroup::create_group(hid_t file_id, const char* group_name,
                  int nsd, ElemType elem_type,
                  const std::vector<std::string>& attribute_names,
                  unsigned int total_node_count,
                  unsigned int total_elem_count) {
  close();

  // create our "root" group
  _root = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (_root < 0) {
    throw FileIOException() <<
        "Cannot create HDF5 group " << group_name <<
        " - it probably already exists. "
        "Try removing the old file or using a different title?";
  }

  // set the "grid_type" attribute
  hsize_t grid_type_dim[] = { 1 };
  hid_t space = H5Screate_simple(1, grid_type_dim, NULL);

  hid_t nsd_attr = H5Acreate(_root, "nsd", H5T_NATIVE_INT, space,
                                   H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(nsd_attr, H5T_NATIVE_INT, &nsd);
  H5Aclose(nsd_attr);

  hid_t elem_type_attr = H5Acreate(_root, "elem_type", H5T_NATIVE_INT, space,
                                   H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(elem_type_attr, H5T_NATIVE_INT, &elem_type);
  H5Aclose(elem_type_attr);

  H5Sclose(space);

  // find nsd and nodes per element from grid_type
  unsigned int nodes_per_elem = get_nodes_in_element(elem_type);
  create_datasets(nsd, attribute_names, total_node_count, nodes_per_elem,
                  total_elem_count);

  // cache grid_type, nsd, and nodes_per_elem since they're used in writing
  _elem_type = elem_type;
  _nsd = nsd;
  _nodes_per_elem = nodes_per_elem;
}

void H5OutputGroup::close() {
  if (_root < 0)
    return;

  // close datasets
  H5Dclose(_coord_dset);
  _coord_dset = -1;

  for (unsigned int i = 0; i < _attribs.size(); i++) {
    H5Dclose(_attribs.at(i).dset);
  }
  _attribs.clear();

  H5Dclose(_elem_dset);
  _elem_dset = -1;

  H5Gclose(_root);
  _root = -1;
}

void H5OutputGroup::create_datasets(unsigned int nsd,
    const std::vector<std::string>& attrib_names, unsigned int total_nodes,
    unsigned int nodes_per_elem, unsigned int total_elems) {
  hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);

  // first, we create the dataset for coordinates (nsd * number of nodes large)
  hsize_t dims[2] = { total_nodes, nsd };
  hid_t filespace = H5Screate_simple(2, dims, NULL);

  _coord_dset = H5Dcreate(_root, "node_coords", NODE_COORDINATE_TYPE, filespace,
                          H5P_DEFAULT, plist_id, H5P_DEFAULT);
  assert(_coord_dset >= 0);
  H5Sclose(filespace);

  // create the "node_data" group, where we store node attributes
  hid_t node_data_group = H5Gcreate(_root, "node_data", H5P_DEFAULT,
                                    H5P_DEFAULT, H5P_DEFAULT);
  assert(node_data_group >= 0);
  H5Gclose(node_data_group);

  // next, we create the datasets for each attribute
  hsize_t attrib_dim[1] = { total_nodes };
  filespace = H5Screate_simple(1, attrib_dim, NULL);
  assert(filespace >= 0);

  std::string attrib_base = "node_data/";
  for (unsigned int i = 0; i < attrib_names.size(); i++) {
    if (attrib_names.at(i).find('/') != std::string::npos)
      throw TALYException() << "Invalid variable name: '" << attrib_names.at(i) << "' - cannot contain /";

    std::string path = attrib_base + attrib_names.at(i);
    hid_t attrib_dset = H5Dcreate(_root, path.c_str(), NODE_DATA_TYPE,
                                  filespace, H5P_DEFAULT, plist_id,
                                  H5P_DEFAULT);
    assert(attrib_dset >= 0);
    _attribs.push_back(H5Attrib(attrib_names.at(i), attrib_dset));
  }
  H5Sclose(filespace);

  // finally, we create the dataset for element connectivity data
  hsize_t elem_dim[2] = { total_elems, nodes_per_elem };
  filespace = H5Screate_simple(2, elem_dim, NULL);
  assert(filespace >= 0);
  _elem_dset = H5Dcreate(_root, "elem_connectivity",
                         ELEM_CONNECTIVITY_DATA_TYPE, filespace, H5P_DEFAULT,
                         plist_id, H5P_DEFAULT);
  assert(_elem_dset >= 0);
  H5Sclose(filespace);

  H5Pclose(plist_id);
}


PetscInt H5OutputGroup::node_count() const {
  hid_t dspace = H5Dget_space(_coord_dset);

  hsize_t dims[2];
  H5Sget_simple_extent_dims(dspace, dims, NULL);
  H5Sclose(dspace);

  return dims[0];
}

PetscInt H5OutputGroup::elem_count() const {
  hid_t dspace = H5Dget_space(_elem_dset);

  hsize_t dims[2];
  H5Sget_simple_extent_dims(dspace, dims, NULL);
  H5Sclose(dspace);

  return dims[0];
}





H5Writer::H5Writer(bool parallel) {
  _parallel = parallel;
  _file_id = -1;
}

H5Writer::~H5Writer() {
  close();
}

void H5Writer::open_file(const char* filename, bool append) {
  // if a file is already open, close it
  if (_file_id >= 0) {
    close();
  }

  // create file access property list
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  if (_parallel) {
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  }

  // mute errors for H5Fis_hdf5 call (in case file doesn't exist)
  H5E_auto_t error_func;
  void* error_cl_data;
  H5Eget_auto(H5E_DEFAULT, &error_func, &error_cl_data);
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  // check file exists and is an hdf5 file before trying to append
  htri_t hdf5_status = H5Fis_hdf5(filename);

  // unmute error messages
  H5Eset_auto(H5E_DEFAULT, error_func, error_cl_data);

  // if we aren't appending to an existing HDF5 file, create a new blank file
  // (using truncate mode)
  if (!append || hdf5_status <= 0) {
    _file_id = H5Fcreate(filename, H5F_ACC_TRUNC,
                         H5P_DEFAULT, plist_id);
  } else {
    // otherwise, open the existing file
    _file_id = H5Fopen(filename, H5F_ACC_RDWR, plist_id);
  }

  H5Pclose(plist_id);  // free the properties list

  // remember the file path for when we need to write XDMF files
  _file_path = filename;
}

void H5Writer::create_group(const char* group_name, int nsd, ElemType elem_type,
                            const std::vector<std::string>& attribute_names,
                            unsigned int total_node_count,
                            unsigned int total_elem_count) {
  if (_file_id < 0) {
    throw FileIOException() << "Cannot create group before opening file!";
  }
  _out_group.create_group(_file_id, group_name, nsd, elem_type, attribute_names,
                          total_node_count, total_elem_count);
}

void H5Writer::close() {
  _out_group.close();

  // close file
  if (_file_id >= 0) {
    herr_t res = H5Fclose(_file_id);
    _file_id = -1;
    _file_path = "";
  }
}

void H5Writer::write_node_data_multi(PetscInt* node_offsets,
                                     PetscInt node_count,
                                     const double* node_coords,
                                     const double** node_attribs,
                                     bool collective) {
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (_parallel) {
    H5Pset_dxpl_mpio(plist_id,
                     collective ? H5FD_MPIO_COLLECTIVE : H5FD_MPIO_INDEPENDENT);
  }

  const unsigned int& nsd = _out_group.nsd();

  // first write coordinate data

  // create memory space
  hsize_t memspace_coord_count[2] = {
    static_cast<hsize_t>(node_count),
    static_cast<hsize_t>(nsd)
  };

  hid_t memspace = H5Screate_simple(2, memspace_coord_count, NULL);
  assert(memspace >= 0);

  // set up the filespace (where we're writing to)
  hid_t filespace = H5Dget_space(_out_group.coord_dset());
  assert(filespace >= 0);

  std::vector <hsize_t> coord_offsets;
  coord_offsets.resize(node_count * nsd * 2);
  for (PetscInt i = 0; i < node_count; i++) {
    for (unsigned int d = 0; d < nsd; d++) {
      coord_offsets[(i * nsd + d) * 2 + 0] = node_offsets[i];
      coord_offsets[(i * nsd + d) * 2 + 1] = d;
    }
  }
  herr_t status = H5Sselect_elements(filespace, H5S_SELECT_SET,
                                     node_count * nsd, &coord_offsets[0]);
  UNUSED(status);
  assert(status >= 0);

  status = H5Dwrite(_out_group.coord_dset(), H5T_NATIVE_DOUBLE, memspace,
                    filespace, plist_id, node_coords);
  assert(status >= 0);

  H5Sclose(memspace);
  H5Sclose(filespace);

  // then write node attribute data
  hsize_t memspace_attrib_count[1] = { static_cast<hsize_t>(node_count) };
  memspace = H5Screate_simple(1, memspace_attrib_count, NULL);
  assert(memspace >= 0);

  // set up the filespace
  std::vector<hsize_t> attrib_offsets;
  attrib_offsets.resize(node_count);
  for (PetscInt node = 0; node < node_count; node++) {
    attrib_offsets[node] = node_offsets[node];
  }
  for (unsigned int i = 0; i < _out_group.attribs().size(); i++) {
    const hid_t& attrib_dset = _out_group.attribs().at(i).dset;
    filespace = H5Dget_space(attrib_dset);
    assert(filespace >= 0);
    status = H5Sselect_elements(filespace, H5S_SELECT_SET, node_count,
                                &attrib_offsets[0]);
    assert(status >= 0);

    // write this attribute
    status = H5Dwrite(attrib_dset, H5T_NATIVE_DOUBLE, memspace, filespace,
                      plist_id, node_attribs[i]);
    assert(status >= 0);
    H5Sclose(filespace);
  }

  H5Sclose(memspace);

  H5Pclose(plist_id);
}

void H5Writer::write_node_data(PetscInt node_offset, PetscInt node_count,
                               const double* node_coords,
                               const double* node_attribs, bool collective) {
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (_parallel) {
    H5Pset_dxpl_mpio(plist_id,
                     collective ? H5FD_MPIO_COLLECTIVE : H5FD_MPIO_INDEPENDENT);
  }

  const unsigned int& nsd = _out_group.nsd();

  hsize_t coord_count[2];
  hsize_t coord_offset[2];
  coord_count[0] = node_count;
  coord_count[1] = nsd;
  coord_offset[0] = node_offset;
  coord_offset[1] = 0;

  // first write coordinate data
  hid_t memspace = H5Screate_simple(2, coord_count, NULL);

  // select hyperslab
  hid_t filespace = H5Dget_space(_out_group.coord_dset());
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, coord_offset, NULL,
                      coord_count, NULL);

  herr_t status = H5Dwrite(_out_group.coord_dset(), H5T_NATIVE_DOUBLE, memspace,
                           filespace, plist_id, node_coords);
  UNUSED(status);
  assert(status >= 0);

  H5Sclose(memspace);
  H5Sclose(filespace);

  // then write node attribute data
  hsize_t attrib_offset[1] = { static_cast<hsize_t>(node_offset) };
  hsize_t attrib_count[1] = { static_cast<hsize_t>(node_count) };
  for (unsigned int i = 0; i < _out_group.attribs().size(); i++) {
    memspace = H5Screate_simple(1, attrib_count, NULL);

    const hid_t& attrib_dset = _out_group.attribs().at(i).dset;

    // select hyperslab
    filespace = H5Dget_space(attrib_dset);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, attrib_offset, NULL,
                        attrib_count, NULL);

    status = H5Dwrite(attrib_dset, H5T_NATIVE_DOUBLE, memspace, filespace,
                      plist_id, node_attribs + i);
    assert(status >= 0);

    H5Sclose(memspace);
    H5Sclose(filespace);
  }

  H5Pclose(plist_id);
}

void H5Writer::write_elem_data(PetscInt elem_offset, PetscInt elem_count,
                               const PetscInt* elem_data, bool collective) {
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (_parallel)
    H5Pset_dxpl_mpio(plist_id,
                     collective ? H5FD_MPIO_COLLECTIVE : H5FD_MPIO_INDEPENDENT);

  hsize_t count[2] = {
    static_cast<hsize_t>(elem_count),
    static_cast<hsize_t>(_out_group.nodes_per_elem())
  };
  hsize_t offset[2] = {
    static_cast<hsize_t>(elem_offset),
    0
  };

  hid_t memspace = H5Screate_simple(2, count, NULL);

  // select hyperslab
  hid_t filespace = H5Dget_space(_out_group.elem_dset());
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

  herr_t status = H5Dwrite(_out_group.elem_dset(), ELEM_CONNECTIVITY_INPUT_TYPE,
                           memspace, filespace, plist_id, elem_data);
  UNUSED(status);
  assert(status >= 0);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Pclose(plist_id);
}


void H5Writer::write_elem_data_multi(PetscInt* elem_offsets, PetscInt elem_count,
                                     const PetscInt* elem_data, bool collective) {
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (_parallel) {
    H5Pset_dxpl_mpio(plist_id,
                     collective ? H5FD_MPIO_COLLECTIVE : H5FD_MPIO_INDEPENDENT);
  }

  const unsigned int& nodes_per_elem = _out_group.nodes_per_elem();

  hsize_t count[2] = {
      static_cast<hsize_t>(elem_count),
      static_cast<hsize_t>(nodes_per_elem)
  };

  hid_t memspace = H5Screate_simple(2, count, NULL);
  assert(memspace >= 0);

  // set up the filespace (where we're writing to)
  hid_t filespace = H5Dget_space(_out_group.elem_dset());
  assert(filespace >= 0);

  std::vector <hsize_t> offsets;
  offsets.resize(elem_count * nodes_per_elem * 2);
  for (PetscInt i = 0; i < elem_count; i++) {
    for (unsigned int d = 0; d < nodes_per_elem; d++) {
      offsets[(i * nodes_per_elem + d) * 2 + 0] = elem_offsets[i];
      offsets[(i * nodes_per_elem + d) * 2 + 1] = d;
    }
  }
  herr_t status = H5Sselect_elements(filespace, H5S_SELECT_SET,
                                     elem_count * nodes_per_elem, &offsets[0]);
  UNUSED(status);
  assert(status >= 0);

  status = H5Dwrite(_out_group.elem_dset(), ELEM_CONNECTIVITY_INPUT_TYPE,
                    memspace, filespace, plist_id, elem_data);
  assert(status >= 0);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Pclose(plist_id);
}

// The best resource I found on the XDMF format is here:
// http://www.xdmf.org/index.php/XDMF_Model_and_Format
// Basically, it's an XML document, laid out like this:
// <Domain>
//  (for each output we want to display, similar to Tecplot's zone concept:)
//  <Grid ...>
//   <Topology ...>
//    ... element connectivity data ...
//   </Topology>
//   <Geometry ...>
//    ... node coordinate data ...
//   </Geometry>
//   (and then, for EACH attribute we want to display:)
//   <Attribute ...>
//    ... node attribute data ...
//   </Attribute>
//  </Grid>
// </Domain>
// A lot of attribute data has to be filled in at write time; most notably
// the dimensions for the data we want to view (which is redundant,
// as the HDF5 file already has them).
// IMPORTANT NOTE: Dimensions in XDMF files are written slowest varying
// dimension first!
// For example, "ZYX".

/*
  This currently doesn't work for the 1D case because, as far as I can tell,
  XDMF only supports 2D and 3D data. In order to "fake" the 1D data being 2D,
  we need to make a grid that specifies the Y value as node data.  This means
  creating multiple <Grid>s, one for each piece of node data.

  If anyone needs to add this in the future, here are some hints:
   - use TopologyType="Polyline".
   - use the non-interlaced GeometryType="X_Y_Z" (note underscores).
     This will let you specify the X/Y/Z values as separate arrays, like this:
      "<DataItem Format=\"HDF\" Dimensions=\"NumberOfNodes\">\n"   // X
      "%s:/%s/node_coords\n"
      "</DataItem>\n"
      "<DataItem Format=\"HDF\" Dimensions=\"NumberOfNodes\">\n"   // Y
      "%s:/%s/node_data/%s\n"
      "</DataItem>\n"
      "<DataItem Format=\"XML\" Dimensions=\"NumberOfNodes\">\n"   // Z
      "0.0\n"
      "0.0\n"
      ...(repeat %d times)...
      "</DataItem>\n"
   - don't bother outputting <Attribute>s
     (all it will do is color the line, the Y values are what's important)
     Z must be present even if it's all zeroes
*/

void H5Writer::write_xdmf() {
  if (_file_id < 0) {
    throw FileIOException() <<
        "write_xdmf should be called while the HDF5 file is still open!";
  }

  std::string path = _file_path;
  path += ".xdmf";
  FILE* f = fopen(path.c_str(), "w");
  if (!f) {
    throw FileIOException() <<
        "Could not open \"" << path << "\" to write XDMF file!";
  }

  fprintf(f,
    "<?xml version=\"1.0\" ?>\n"
    "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
    "<Xdmf Version=\"2.0\">\n"
    "<Domain>\n"
    "<Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n");

  // each group has its own grid
  group_vec_t groups;
  H5Giterate(_file_id, ".", 0, &list_h5_groups, &groups);

  for (unsigned int i = 0; i < groups.size(); i++) {
    const std::string& group_name = groups.at(i);

    H5OutputGroup group;
    group.open(_file_id, group_name.c_str());

    const ElemType elem_type = group.elem_type();
    const unsigned int nsd = group.nsd();
    const PetscInt node_count = group.node_count();
    const unsigned int nodes_per_elem = group.nodes_per_elem();
    const PetscInt elem_count = group.elem_count();

    // find the TopologyType name from grid type
    const char* topology_type = NULL;
    switch (elem_type) {
      case kElem3dHexahedral: topology_type = "Hexahedron"; break;
      case kElem2dBox: topology_type = "Quadrilateral"; break;
      case kElem2dTriangle: topology_type = "Triangle"; break;
      case kElem3dTetrahedral: topology_type = "Tetrahedron"; break;
      case kElem1d: topology_type = "Polyline"; break;
      default: throw NotImplementedException() << "HDF5: unknown element type";
    }

    const char* geometry_type = NULL;
    switch (nsd) {
      case 3: geometry_type = "XYZ"; break;
      case 2: geometry_type = "XY"; break;
      default: PrintWarning("Unsupported NSD for XDMF output "
          "(nsd = ", nsd, "), ignoring group ", group_name.c_str());
        continue;
    }

    PetscInt geometry_dimX = nsd;
    PetscInt geometry_dimY = node_count;

    // start actually writing to the file

    fprintf(f, "<Grid Name=\"%s\">\n", group_name.c_str());

    // need to replace
    // TopologyType, NumberOfElements,
    // number of elements, nodes per element,
    // file name, group name
    const char* topology =
      "<Topology TopologyType=\"%s\" NumberOfElements=\"%d\">\n"
      "<DataItem Format=\"HDF\" Dimensions=\"%d %d\">\n"
      "%s:/%s/elem_connectivity\n"
      "</DataItem>\n"
      "</Topology>\n";
    fprintf(f, topology,
            topology_type, elem_count, elem_count, nodes_per_elem,
            _file_path.c_str(), group_name.c_str());

    // need to replace
    // geometry type (XY or XYZ),
    // coordinate dimensions (1 if nsd == 3), nodes, nsd,
    // file name, group name
    const char* geometry =
      "<Geometry GeometryType=\"%s\">\n"
      "<DataItem Format=\"HDF\" Dimensions=\"%d %d\">\n"
      "%s:/%s/node_coords\n"
      "</DataItem>\n"
      "</Geometry>\n";

    fprintf(f, geometry,
          geometry_type,
          geometry_dimY, geometry_dimX, /* this order is not a typo! */
          _file_path.c_str(), group_name.c_str());

    // need to replace
    // attribute name, number of nodes,
    // file name, group name, attribute name
    const char* attrib_text =
      "<Attribute Name=\"%s\" Center=\"Node\">\n"
      "<DataItem Format=\"HDF\" Dimensions=\"%d\">\n"
      "%s:/%s/node_data/%s\n"
      "</DataItem>\n"
      "</Attribute>\n";

    for (unsigned int j = 0; j < group.attribs().size(); j++) {
      const std::string& attrib_name = group.attribs().at(j).name;
      fprintf(f, attrib_text,
              attrib_name.c_str(), node_count, _file_path.c_str(),
              group_name.c_str(),
              attrib_name.c_str());
    }

    if (group_name.find("t=") == 0) {
      std::string ts = group_name.substr(2);
      fprintf(f, "<Time Value=\"%s\" />\n", ts.c_str());
    }

    fprintf(f, "</Grid>\n");
  }

  fprintf(f,
    "</Grid>\n"
    "</Domain>\n"
    "</Xdmf>\n");

  fclose(f);
}



// Loading

H5Reader::H5Reader(bool parallel) {
  _parallel = parallel;
  _file_id = -1;
}

H5Reader::~H5Reader() {
  close();
}

void H5Reader::open_file(const char* filename) {
  // if a file is already open, close it
  if (_file_id >= 0)
    close();

  // create file access property list
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  if (_parallel)
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

  _file_id = H5Fopen(filename, H5F_ACC_RDONLY, plist_id);

  H5Pclose(plist_id);  // free the properties list
}

void H5Reader::open_group(const char* group_name) {
  if (_file_id < 0) {
    throw FileIOException() << "File must be open before opening a group!";
  }
  _in_group.open(_file_id, group_name);
}

void H5Reader::read_node_data(PetscInt node_offset, PetscInt node_count,
                              double* node_coords_out,
                              double* node_attribs_out) {
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (_parallel)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

  const unsigned int& nsd = _in_group.nsd();

  hsize_t coord_count[2];
  hsize_t coord_offset[2];
  coord_count[0] = node_count;
  coord_count[1] = nsd;
  coord_offset[0] = node_offset;
  coord_offset[1] = 0;

  // first write coordinate data
  hid_t memspace = H5Screate_simple(2, coord_count, NULL);

  // select hyperslab
  hid_t filespace = H5Dget_space(_in_group.coord_dset());
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, coord_offset, NULL,
                      coord_count, NULL);

  herr_t status = H5Dread(_in_group.coord_dset(), H5T_NATIVE_DOUBLE, memspace,
                          filespace, plist_id, node_coords_out);
  UNUSED(status);
  assert(status >= 0);

  H5Sclose(memspace);
  H5Sclose(filespace);

  // then write node attribute data
  hsize_t attrib_offset[1] = { static_cast<hsize_t>(node_offset) };
  hsize_t attrib_count[1] = { static_cast<hsize_t>(node_count) };
  for (unsigned int i = 0; i < _in_group.attribs().size(); i++) {
    memspace = H5Screate_simple(1, attrib_count, NULL);

    const hid_t& attrib_dset = _in_group.attribs().at(i).dset;

    // select hyperslab
    filespace = H5Dget_space(attrib_dset);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, attrib_offset, NULL,
                        attrib_count, NULL);

    status = H5Dread(attrib_dset, H5T_NATIVE_DOUBLE, memspace, filespace,
                     plist_id, node_attribs_out + i);
    assert(status >= 0);

    H5Sclose(memspace);
    H5Sclose(filespace);
  }

  H5Pclose(plist_id);
}

void H5Reader::read_elem_data(PetscInt elem_offset, PetscInt elem_count,
                              PetscInt* elem_data_out) {
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (_parallel)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

  hsize_t count[2] = {
    static_cast<hsize_t>(elem_count),
    static_cast<hsize_t>(_in_group.nodes_per_elem())
  };
  hsize_t offset[2] = {
    static_cast<hsize_t>(elem_offset),
    0
  };

  hid_t memspace = H5Screate_simple(2, count, NULL);

  // select hyperslab
  hid_t filespace = H5Dget_space(_in_group.elem_dset());
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

  herr_t status = H5Dread(_in_group.elem_dset(), ELEM_CONNECTIVITY_INPUT_TYPE,
                          memspace, filespace, plist_id, elem_data_out);
  UNUSED(status);
  assert(status >= 0);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Pclose(plist_id);
}

void H5Reader::close() {
  _in_group.close();

  // close file
  if (_file_id >= 0) {
    H5Fclose(_file_id);
    _file_id = -1;
  }
}

}  // namespace TALYFEMLIB
#endif
