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
#include <talyfem/file_io/gmsh_io.h>

#include <assert.h>

#include <functional>  // for std::function
#include <limits>  // for std::numeric_limits
#include <string>
#include <vector>

#include <talyfem/common/exceptions.h>
#include <talyfem/grid/elem_common.h>
#include <talyfem/grid/nodeindicator.h>

namespace gmsh {

TALYFEMLIB::ElemType gmsh_elm_to_taly_elm(gmsh::ElementType t) {
  switch (t) {
    case LINE_ORDER1:
    case LINE_ORDER2:
    case LINE_ORDER3:
    case LINE_ORDER4:
      return TALYFEMLIB::kElem1d;

    case BOX_ORDER1:
    case BOX_ORDER2:
    case BOX_ORDER3:
      return TALYFEMLIB::kElem2dBox;

    case TRIANGLE_ORDER1:
    case TRIANGLE_ORDER2:
    case TRIANGLE_9NODE:
    case TRIANGLE_10NODE:
        return TALYFEMLIB::kElem2dTriangle;

    case HEXAHEDRON_ORDER1:
    case HEXAHEDRON_ORDER2:
    case HEXAHEDRON_ORDER3:
    case HEXAHEDRON_ORDER4:
      return TALYFEMLIB::kElem3dHexahedral;

    case TETRAHEDRON_ORDER1:
    case TETRAHEDRON_ORDER2:
    case TETRAHEDRON_ORDER3:
    case TETRAHEDRON_ORDER4:
      return TALYFEMLIB::kElem3dTetrahedral;

    default:
      throw TALYFEMLIB::NotImplementedException() << "gmsh_elm_to_taly_elm not "
            "implemented for gmsh ElementType " << t;
  }
}

int nodes_per_element(ElementType type) {
  switch (type) {
    case POINT: return 1;

    case LINE_ORDER1: return 2;
    case LINE_ORDER2: return 3;
    case LINE_ORDER3: return 4;
    case LINE_ORDER4: return 5;

    case TRIANGLE_ORDER1: return 3;
    case TRIANGLE_ORDER2: return 6;
    case TRIANGLE_9NODE: return 9;
    case TRIANGLE_10NODE: return 10;

    case BOX_ORDER1: return 4;
    case BOX_ORDER2: return 9;
    case BOX_ORDER3: return 16;

    case TETRAHEDRON_ORDER1: return 4;
    case TETRAHEDRON_ORDER2: return 10;
    case TETRAHEDRON_ORDER3: return 20;
    case TETRAHEDRON_ORDER4: return 56;

    case HEXAHEDRON_ORDER1: return 8;
    case HEXAHEDRON_ORDER2: return 27;
    case HEXAHEDRON_ORDER3: return 64;
    case HEXAHEDRON_ORDER4: return 125;

    default:
      throw TALYFEMLIB::NotImplementedException() << "Unknown nodes per element for gmsh "
                                         "element type '" << type << "'";
  }
}

int nsd_of_element(ElementType type) {
  switch (type) {
    case LINE_ORDER1:
    case LINE_ORDER2:
    case LINE_ORDER3:
    case LINE_ORDER4:
      return 1;

    case TRIANGLE_ORDER1:
    case TRIANGLE_ORDER2:
    case TRIANGLE_9NODE:
    case TRIANGLE_10NODE:
    case BOX_ORDER1:
    case BOX_ORDER2:
    case BOX_ORDER3:
      return 2;

    case TETRAHEDRON_ORDER1:
    case TETRAHEDRON_ORDER2:
    case TETRAHEDRON_ORDER3:
    case TETRAHEDRON_ORDER4:
    case HEXAHEDRON_ORDER1:
    case HEXAHEDRON_ORDER2:
    case HEXAHEDRON_ORDER3:
    case HEXAHEDRON_ORDER4:
      return 3;

    default:
      throw TALYFEMLIB::NotImplementedException() << "Unknown nsd for gmsh element type "
                                         "'" << type << "'";
  }
}

/**
 * Gets the ElementType for a "surface" for type.
 * This is not defined by the Gmsh standard.
 * Used to calculate "primary_surf_type".
 */
ElementType surface_type_of_element(ElementType type) {
  switch (type) {
    case LINE_ORDER1:
    case LINE_ORDER2:
    case LINE_ORDER3:
      return POINT;

    case TRIANGLE_ORDER1:
    case BOX_ORDER1:
      return LINE_ORDER1;

    case TRIANGLE_ORDER2:
    case BOX_ORDER2:
      return LINE_ORDER2;

    case TRIANGLE_9NODE:
    case TRIANGLE_10NODE:
      return LINE_ORDER3;

    case BOX_ORDER3:
      return LINE_ORDER3;

    case TETRAHEDRON_ORDER1:
      return TRIANGLE_ORDER1;
    case TETRAHEDRON_ORDER2:
      return TRIANGLE_ORDER2;
    case TETRAHEDRON_ORDER3:
      return TRIANGLE_10NODE;

    case HEXAHEDRON_ORDER1:
      return BOX_ORDER1;
    case HEXAHEDRON_ORDER2:
      return BOX_ORDER2;
    case HEXAHEDRON_ORDER3:
      return BOX_ORDER3;

    default:
      throw TALYFEMLIB::NotImplementedException() << "Unknown surface type for gmsh "
                                         "element type '" << type << "'";
  }
}

int get_order(ElementType type) {
  switch (type) {
    case LINE_ORDER1:
    case TRIANGLE_ORDER1:
    case BOX_ORDER1:
    case TETRAHEDRON_ORDER1:
    case HEXAHEDRON_ORDER1:
      return 1;

    case LINE_ORDER2:
    case TRIANGLE_ORDER2:
    case BOX_ORDER2:
    case TETRAHEDRON_ORDER2:
    case HEXAHEDRON_ORDER2:
      return 2;

    case LINE_ORDER3:
    case TRIANGLE_9NODE:
    case TRIANGLE_10NODE:
    case BOX_ORDER3:
    case TETRAHEDRON_ORDER3:
    case HEXAHEDRON_ORDER3:
      return 3;

    case LINE_ORDER4:
    case TETRAHEDRON_ORDER4:
    case HEXAHEDRON_ORDER4:
      return 4;

    default:
      throw TALYFEMLIB::NotImplementedException() << "Unknown order of gmsh element type "
                                         "'" << type << "'";
  }
}

ElementType to_order_1(ElementType type) {
  switch (type) {
    case POINT:
      return POINT;

    case LINE_ORDER1:
    case LINE_ORDER2:
    case LINE_ORDER3:
    case LINE_ORDER4:
      return LINE_ORDER1;

    case TRIANGLE_ORDER1:
    case TRIANGLE_ORDER2:
    case TRIANGLE_9NODE:
    case TRIANGLE_10NODE:
      return TRIANGLE_ORDER1;

    case BOX_ORDER1:
    case BOX_ORDER2:
    case BOX_ORDER3:
      return BOX_ORDER1;

    case TETRAHEDRON_ORDER1:
    case TETRAHEDRON_ORDER2:
    case TETRAHEDRON_ORDER3:
    case TETRAHEDRON_ORDER4:
      return TETRAHEDRON_ORDER1;

    case HEXAHEDRON_ORDER1:
    case HEXAHEDRON_ORDER2:
    case HEXAHEDRON_ORDER3:
    case HEXAHEDRON_ORDER4:
      return HEXAHEDRON_ORDER1;

    default:
      throw TALYFEMLIB::NotImplementedException() << "Unknown order 1 element type for "
                                         "gmsh element type '" << type << "'";
  }
}

void reorder_connectivity_for_taly(ElementType type,
                                   std::vector<NodeID>& connectivity) {
  const int* order_map = NULL;
  if (type == HEXAHEDRON_ORDER1) {
    static const int hex_linear_map[8] = {
      4, 5, 6, 7,
      0, 1, 2, 3
    };
    order_map = hex_linear_map;
  } else if (type == HEXAHEDRON_ORDER2) {
    /*        gmsh layout                             taly layout
     *
     *      3----13----2                           7-------19--------6
     *      |\         |\                          |\       |\       |\
     *      |15    24  | 14                        |11-------24-------10
     *      9  \ 20    11 \                        | |\     | |\     | |\
     *      |   7----19+---6                       | | 3-------14--------2
     *      |22 |  26  | 23|                      20-|-|---21-|-|---18 | |
     *      0---+-8----1   |                       |\| |    |\| |    |\| |
     *       \ 17    25 \  18      ==>             |25-|-----26-|-----23 |
     *       10 |  21    12|                       | |\|    | |\|    | |\|
     *         \|         \|                       | |15-------16-------13
     *          4----16----5                       4-|-|---17-|-|----5 | |
     *                                              \| |     \| |     \| |
     *                                               8-|-----22-|------9 |
     *     (Z+ = out of page)                         \|       \|       \|
     *                                                 0-------12--------1
     */
    static const int hex_quad_map[27] = {
      4, 5, 6, 7,
      0, 1, 2, 3,
      10, 12, 14, 15,
      16, 18, 19, 17,
      25, 8, 11, 13,
      9, 20, 21, 23,
      24, 22, 26
    };
    order_map = hex_quad_map;
  } else if (type == HEXAHEDRON_ORDER3) {
    /*    gmsh
     *    X+ = right, Y+ = up, Z+ = out of page (towards viewer)
     *
     *         Y = 1/4 (bottom face)
     *
     *      --------------------
     *     |\                   |\
     *     | \                  | \
     *     |  \                 |  \
     *     |   \                |   \
     *     |    \               |    \
     *     |     \--------------------\
     *     |     |              |     |
     *     |     |              |     |
     *     |     |              |     |            Y = 2/4
     *     |     |              |     |       21-----51-----50----23
     *     |     |              |     |        \                    \
     *     6-----|-30-----31----7     |        46    62      63     42
     *      \    |               \    |          \                    \
     *      29   |54      55     27   |          47    61      60     41
     *        \  |                 \  |            \                    \
     *        28 |  53      52     26 |             17------38-----39---13
     *          \|                   \|
     *           5-----25------24-----4
     *
     *
     *
     *              Y = 3/4                          Y = 4/4 (top face)
     *         20-----48-----49----22               2-----18-----19------3
     *          \                    \               \                    \
     *          45    58      59     43              15     34      33    11
     *            \                    \               \                    \
     *            44    57      56     40              14     35      32    10
     *              \                    \               \                    \
     *               16------37----36----12               1-----9-------8-----0
     */
    /*
     *         taly
     *
     *      Y = 1/4 (bottom face)                        Y = 2/4
     *     4------28-----29-----5                35-----36-----37----30
     *      \                    \                \                    \
     *      12    52      53     13               59    60      61     54
     *        \                    \                \                    \
     *         8    40      41     9                47    48      49     42
     *          \                    \                \                    \
     *           0------16-----17-----1                23------24----25----18
     *
     *
     *           Y = 3/4                           Y = 4/4 (top face)
     *     34-----39-----38----31                7------33-----32-----6
     *      \                    \                \                    \
     *      58    63      62     55               15    57      56     14
     *        \                    \                \                    \
     *        46    51      50     43               11    45      44     10
     *          \                    \                \                    \
     *           22------27----26----19                3-------21----20-----2
     *
     */
    static const int hex_cubic_map[64] = {
      5, 4, 0, 1,      // 0
      6, 7, 3, 2,      // 4
      28, 26, 10, 14,  // 8
      29, 27, 11, 15,  // 12
      25, 24, 13, 12,  // 16
      8, 9, 16, 17,    // 20
      38, 39, 36, 37,  // 24
      30, 31, 23, 22,  // 28
      19, 18, 20, 21,  // 32
      51, 50, 49, 48,  // 36
      53, 52, 41, 40,  // 40
      32, 35, 44, 47,  // 44
      61, 60, 56, 57,  // 48
      54, 55, 42, 43,  // 52
      33, 34, 45, 46,
      62, 63, 59, 58
    };
    order_map = hex_cubic_map;
  }

  // apply the reordering, if any
  if (order_map != NULL) {
    std::vector<NodeID> new_conn(connectivity.size());
    for (unsigned int i = 0; i < new_conn.size(); i++) {
      new_conn[i] = connectivity[order_map[i]];
    }
    connectivity = std::move(new_conn);
  }
}

TALYFEMLIB::SurfaceIndicator::IndicatorType tags_to_surface_indicators(
    const gmsh::Element& elem) {
  return elem.tags.empty() ? 0 : TALYFEMLIB::INDICATOR_NUM(elem.tags.at(0));
}

Reader::Reader()
  : cur_sec_(NULL), nodes_sec_(NULL), elements_sec_(NULL),
    n_nodes_(0), n_elements_(0), nsd_(-1),
    primary_elm_type_(NONE), primary_surf_type_(NONE) {
}

Reader::~Reader() {
  close();
}

void Reader::open(const std::string& path) {
  if (file_.is_open())
    throw TALYFEMLIB::FileIOException() << "File already open!";

  file_.open(path);
  if (!file_) {
    throw TALYFEMLIB::FileIOException() << "Could not open Gmsh file '" << path << "'. "
                            << "Does the file exist? Check permissions?";
  }

  build_section_list();
  check_section_list();
}

void Reader::close() {
  if (!file_.is_open())
    return;

  file_.close();
  cur_sec_ = NULL;
  sections_.clear();
  elm_type_counts_.clear();
  nodes_sec_ = NULL;
  elements_sec_ = NULL;
  nsd_ = -1;
  primary_elm_type_ = NONE;
  primary_surf_type_ = NONE;
}

// parse_func is an optional function which will be called on each line in the
// section with the current line number (zero-indexed).
// parse_func can read **on the current line** as it wants - it must not
// consume the newline.
void read_section(std::ifstream& file, const std::string& section_name,
                  const std::function<void(unsigned int line_no)>& parse_func) {
  unsigned int line_no = 0;
  while (file.good()) {
    char c;
    if (!file.get(c))
      break;

    if (c == '$') {
      std::string end_section;
      std::getline(file, end_section);
      if (end_section != std::string("End") + section_name) {
        throw TALYFEMLIB::FileIOException() << "Expected $End" << section_name << ", "
                                   "got $" << end_section;
      }

      return;  // done
    } else {
      if (parse_func) {
        // give parse_func a chance to parse the line
        file.unget();
        parse_func(line_no);
        assert(file.good());
      }

      // skip the rest of the line
      file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      line_no++;
    }
  }

  throw TALYFEMLIB::FileIOException() << "Section '" << section_name << "' "
                             "ended unexpectedly";
}

void Reader::build_section_list() {
  assert(sections_.empty());
  assert(elm_type_counts_.empty());

  std::string line;
  while (!file_.eof() && file_.good()) {
    std::getline(file_, line);

    // skip empty lines (not explicitly permitted by the standard...)
    if (line.empty())
      continue;

    if (line.at(0) != '$') {
      throw TALYFEMLIB::FileIOException() << "Expected section header "
                                 "(got '" << line << "')";
    }

    const std::streamoff start = file_.tellg();
    const std::string section_name = line.substr(1);

    if (section_name == "Elements") {
      read_section(file_, "Elements", [this](unsigned int line_no) {
        if (line_no >= 1) {
          ElmID id;
          int type;
          file_ >> id >> type;
          this->elm_type_counts_[(ElementType) type] += 1;
        }
      });
    } else {
      read_section(file_, section_name, NULL);
    }
    std::streamoff end = file_.tellg();

    sections_.push_back(Section { section_name, start, (end - start) });
  }
}

void Reader::seek_section(Section* sec) {
  if (sec) {
    file_.clear();
    file_.seekg(sec->start);
  }

  cur_sec_ = sec;
}

Section* require_one_section(std::vector<Section>& vec, const char* name) {
  Section* result = NULL;

  for (auto it = vec.begin(); it != vec.end(); it++) {
    if (it->name == name) {
      if (result == NULL) {
        result = &*it;
      } else {
        throw TALYFEMLIB::FileIOException() << "Too many $" << name << " sections!";
      }
    }
  }

  if (!result)
    throw TALYFEMLIB::FileIOException() << "Missing required $" << name << " section";

  return result;
}

template <typename T>
void check(std::ifstream& file, const T& expected, const char* name = "value") {
  T val;
  file >> val;

  if (!file.good()) {
    throw TALYFEMLIB::FileIOException() << "Failed to read " << name << " for validation "
                            << "(expected '" << expected << "')";
  }

  if (val != expected) {
    throw TALYFEMLIB::FileIOException() << "Expected " << name << " '" << expected
                            << "', but got '" << val << "'";
  }
}

void Reader::check_section_list() {
  // validate mesh_fmt_sec
  Section* mesh_fmt_sec = require_one_section(sections_, "MeshFormat");
  seek_section(mesh_fmt_sec);
  check<std::string>(file_, "2.2", "MeshFormat version number");
  check<int>(file_, 0, "MeshFormat file type");
  check<int>(file_, sizeof(double), "MeshFormat data size");

  // read the number of nodes
  nodes_sec_ = require_one_section(sections_, "Nodes");
  seek_section(nodes_sec_);
  file_ >> n_nodes_;
  if (!file_.good() || n_nodes_ <= 0)
    throw TALYFEMLIB::FileIOException() << "Failed to read or invalid number of nodes";
  nodes_sec_->start = file_.tellg();

  // read the number of elements
  elements_sec_ = require_one_section(sections_, "Elements");
  seek_section(elements_sec_);
  file_ >> n_elements_;
  if (!file_.good() || n_elements_ <= 0)
    throw TALYFEMLIB::FileIOException() << "Failed to read or invalid number of elements";
  elements_sec_->start = file_.tellg();

  seek_section(NULL);

  // calculate nsd, primary element type, and primary surface type
  primary_elm_type_ = elm_type_counts_.begin()->first;
  for (auto it = elm_type_counts_.begin(); it != elm_type_counts_.end(); it++) {
    if (it->first == POINT)
      continue;

    if (nodes_per_element(it->first) > nodes_per_element(primary_elm_type_)) {
      primary_elm_type_ = it->first;
    }
  }
  nsd_ = nsd_of_element(primary_elm_type_);
  primary_surf_type_ = surface_type_of_element(primary_elm_type_);
}

Node Reader::read_node() {
  if (cur_sec_ != nodes_sec_)
    seek_section(nodes_sec_);

  Node node;
  file_ >> node.number >> node.coords[0] >> node.coords[1] >> node.coords[2];

  if (!file_.good())
    throw TALYFEMLIB::FileIOException() << "Unexpected error reading node";

  return node;
}

Element Reader::read_element() {
  if (cur_sec_ != elements_sec_)
    seek_section(elements_sec_);

  Element elm;
  file_ >> elm.number;

  int type_raw;
  file_ >> type_raw;
  elm.type = static_cast<ElementType>(type_raw);
  if (!file_.good())  // need to make sure elm.type is valid for connectivity
    throw TALYFEMLIB::FileIOException() << "Unexpected error reading element";

  int n_tags = 0;
  file_ >> n_tags;

  // need to make sure n_tags is valid before we start reading tags
  if (!file_.good() && n_tags >= 0) {
    throw TALYFEMLIB::FileIOException() << "Unexpected error reading element "
                            << elm.number << "'s tags";
  }

  elm.tags.resize(n_tags);
  for (int i = 0; i < n_tags; i++) {
    file_ >> elm.tags[i];
  }

  int num_nodes = nodes_per_element(elm.type);
  elm.connectivity.resize(num_nodes);
  for (int i = 0; i < num_nodes; i++) {
    file_ >> elm.connectivity[i];
  }

  // reorder connectivity to match what Taly uses
  reorder_connectivity_for_taly(elm.type, elm.connectivity);

  if (!file_.good()) {
    throw TALYFEMLIB::FileIOException() << "Error reading element " << elm.number;
  }

  return elm;
}

}  // namespace gmsh
