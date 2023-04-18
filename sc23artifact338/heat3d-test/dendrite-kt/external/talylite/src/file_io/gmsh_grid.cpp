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
#include <talyfem/file_io/gmsh_grid.h>

#include <algorithm>
#include <map>
#include <vector>

#include <talyfem/file_io/gmsh_io.h>
#include <talyfem/file_io/common.h>
#include <talyfem/grid/elem-types.h>  // for make_elem_of_type()
#include <talyfem/grid/node.h>
#include <talyfem/grid/nodeid_types.h>
#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/common/exceptions.h>
#include <talyfem/domain_decomposition/mesh_partition.h>
#include <talyfem/grid/grid_common.h>

namespace TALYFEMLIB {

bool is_gmsh_file(const char* path) {
  return has_ending(path, ".msh");
}

// This version uses ~1/3 as much memory per surface, but at the cost of
// being much harder to understand and less robust: https://pastebin.com/cqWpaEBV
class GmshSurfaceMapper {
 public:
  void add(const gmsh::Element& e) {
    // sort so we have a consistent ordering (since it's hard to make gmsh
    // always give counter-clockwise surfaces)
    // We only sort the order 1 nodes (vertices), because Taly only supports
    // order-1 surfaces.
    std::vector<gmsh::NodeID> sorted_connectivity = e.connectivity;
    std::sort(sorted_connectivity.begin(),
              sorted_connectivity.begin() + sorted_connectivity.size());

    // make 1-indexed
    for (unsigned int j = 0; j < sorted_connectivity.size(); j++) {
      sorted_connectivity.at(j) -= 1;
    }

    auto flags = gmsh::tags_to_surface_indicators(e);
    map_.insert({ {sorted_connectivity}, flags});
  }

  std::vector<SurfaceIndicator> find(const ELEM* elem) {
    std::vector<SurfaceIndicator> out;

    const int* full_surf_array = elem->GetSurfaceCheckArray();
    const int nodes_per_surf = elem->GetNodesPerSurface();
    const int num_surfaces = elem->GetSurfaceCount();

    Key key;
    key.nodes.resize(nodes_per_surf);

    // loop through all possible surfaces on this element
    for (int i = 0; i < num_surfaces; i++) {
      // surf_id_ptr points to the surface ID, followed by a list of
      // ElemNodeIDs for this surface (the +1s are to skip the surface ID)
      const int* surf_id_ptr = full_surf_array
                               + i * elem->GetSurfaceCheckArrayRowLength();
      const int* elm_node_ids = surf_id_ptr + 1;
      for (int j = 0; j < nodes_per_surf; j++) {
        key.nodes.at(j) = elem->ElemToLocalNodeID(elm_node_ids[j]);
      }

      // sort for consistent comparison
      std::sort(key.nodes.begin(), key.nodes.end());

      // so, does a surface exist for these nodes?
      auto it = map_.find(key);

      // no surface found
      if (it == map_.end())
        continue;

      // surface exists, create it
      {
        SurfaceIndicator surf_indicator(*surf_id_ptr);
        surf_indicator.set_indicators(it->second);
        out.push_back(surf_indicator);
      }
    }

    return out;
  }

  inline size_t n_surfaces() const {
    return map_.size();
  }

 protected:

  struct Key {
    std::vector<gmsh::NodeID> nodes;

    inline bool operator<(const Key& rhs) const {
      const Key& lhs = *this;

      if (lhs.nodes.size() < rhs.nodes.size())
        return true;
      if (lhs.nodes.size() > rhs.nodes.size())
        return false;

      for (std::vector<gmsh::NodeID>::size_type i = 0;
           i < lhs.nodes.size(); i++) {
        gmsh::NodeID a = lhs.nodes.at(i);
        gmsh::NodeID b = rhs.nodes.at(i);
        if (a == b)
          continue;

        return (a < b);
      }
      return false;
    }
  };

  std::map<Key, SurfaceIndicator::IndicatorType> map_;
};

void gmsh_load_grid(GRID* grid, const char* path, const char* load_id,
                    bool load_indicators) {
  gmsh::Reader reader;
  reader.open(path);

  const auto& elm_types = reader.elm_type_counts();
  std::vector< std::vector<gmsh::ElementType> > type_groups = {
      {gmsh::LINE_ORDER1},
      {gmsh::BOX_ORDER1, gmsh::TRIANGLE_ORDER1},
      {gmsh::HEXAHEDRON_ORDER1, gmsh::TETRAHEDRON_ORDER1},
  };

  std::vector<gmsh::ElementType>* accepted_types = NULL;
  for (auto it = type_groups.begin(); it != type_groups.end(); it++) {
    // is this element type in the group?
    if (std::find(it->begin(), it->end(), gmsh::to_order_1(reader.primary_elm_type())) != it->end()) {
      accepted_types = &*it;
    }
  }

  assert(accepted_types != NULL);

  std::vector<gmsh::ElementType> accepted_surf_types;
  for (auto it = accepted_types->begin(); it != accepted_types->end(); it++) {
    accepted_surf_types.push_back(gmsh::surface_type_of_element(*it));
  }

  const PetscInt n_nodes = reader.n_nodes();

  // lambda that returns true if a given type should be included in the taly mesh
  auto is_vol_elm = [&accepted_types] (gmsh::ElementType t) -> bool {
    return (std::find(accepted_types->begin(), accepted_types->end(), gmsh::to_order_1(t)) != accepted_types->end());
  };

  // lambda that returns true if a given type is a surface elm type
  auto is_surf_elm = [&accepted_surf_types] (gmsh::ElementType t) -> bool {
    return (std::find(accepted_surf_types.begin(), accepted_surf_types.end(), gmsh::to_order_1(t)) != accepted_surf_types.end());
  };

  int max_order = 1;
  int n_total_elements = 0;
  for (auto it = elm_types.begin(); it != elm_types.end(); it++) {
    if (is_vol_elm(it->first)) {
      n_total_elements += it->second;
      max_order = std::max(max_order, gmsh::get_order(it->first));
    }
  }

  // always 3, as our coords could be in any dimension...
  const int nsd = 3;

  grid->set_nsd(nsd);
  grid->set_basis_order(max_order);
  grid->set_n_nodes(n_nodes);
  grid->set_n_elements(n_total_elements);
  grid->redimArrays(grid->n_nodes(), grid->n_elements());

  // create nodes
  for (LocalNodeID i = 0; i < n_nodes; i++) {
    const gmsh::Node data = reader.read_node();

    NODE* node = new NODE();
    node->setCoor(data.coords[0], data.coords[1], data.coords[2]);

    // gmsh is 1-indexed, we are 0-indexed
    grid->node_array_[data.number - 1] = node;
  }

  GmshSurfaceMapper surfaces;

  // create elements
  unsigned int next_elm_id = 0;
  for (unsigned int i = 0; i < reader.n_elements(); i++) {
    gmsh::Element data = reader.read_element();
    if (is_vol_elm(data.type)) {
      // make connectivity one-indexed
      for (unsigned int j = 0; j < data.connectivity.size(); j++) {
        data.connectivity[j] -= 1;
      }

      ELEM* elem = make_elem_of_type(gmsh::gmsh_elm_to_taly_elm(data.type));
      elem->redim(data.connectivity.size(),
                  phys_to_local(data.connectivity).data());

      elem->set_elm_id(next_elm_id);  // ignoring gmsh element IDs...
      grid->elm_array_[next_elm_id] = elem;
      elem->Validate(grid);
      next_elm_id++;
    } else if (is_surf_elm(data.type) && load_indicators) {
      surfaces.add(data);

      // convert surface tags to flags so we can add the indicators to the node
      // right away
      SurfaceIndicator::IndicatorType add_indicators;
      add_indicators = tags_to_surface_indicators(data);

      for (unsigned int j = 0; j < data.connectivity.size(); j++) {
        LocalNodeID node_id = data.connectivity[j] - 1;
        grid->node_array_[node_id]->addIndicators(add_indicators);
      }
    }
  }

  if (load_indicators) {
    unsigned int n_added_surfaces = 0;
    for (int e = 0; e < grid->n_elements(); e++) {
      ELEM* elem = grid->elm_array_[e];
      elem->surface_indicator_ = surfaces.find(elem);
      n_added_surfaces += elem->surface_indicator_.size();
      elem->CalculateSurfaceNormals(grid);
    }

    grid->SetCaredSurfaceIndicator();

    if (n_added_surfaces != surfaces.n_surfaces()) {
      PrintWarning("Only ", n_added_surfaces, " of ", surfaces.n_surfaces(),
                   " Gmsh surfaces matched volume element surfaces.");
      PrintWarning("Ignoring Gmsh surface data and falling back on "
                   "automatic surface generation based on node indicators "
                   "(ELEM::GenSurfaceIndicator).");
      PrintWarning("If your indicators are not side-specific "
                   "(that is, the same indicator can appear on both the "
                   "left surface and bottom surface), this may generate "
                   "false surfaces in the corners of your mesh.");
      PrintWarning("If your code uses Integrands4side(), you may "
                   "get incorrect results.");

      grid->GenElmSurfaceIndicator();
    }
  }
}

void gmsh_load_grid_dd(GRID* grid, const char* path, const char* load_id,
                       bool load_indicators) {
  CMeshPartition part;
  PetscErrorCode err = part.LoadFromGmshFile(path);
  if (err != 0)
    throw FileIOException() << "Error loading mesh from Gmsh file!";

  part.TransferToGrid(grid);
  part.PartitionFree();
}

}  // namespace TALYFEMLIB

