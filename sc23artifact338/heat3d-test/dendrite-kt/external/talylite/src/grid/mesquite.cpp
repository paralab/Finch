#ifdef ENABLE_MESQUITE

#include <talyfem/grid/mesquite.h>
#include <talyfem/grid/node.h>
#include <talyfem/grid/elem.h>

namespace TALYFEMLIB {

int TalyMesqMesh::get_geometric_dimension(Mesquite2::MsqError &err) {
  return grid_->nsd();
}

void TalyMesqMesh::get_all_elements(std::vector<Mesquite2::Mesh::ElementHandle> &elements,
                                    Mesquite2::MsqError &err) {
  elements.resize(grid_->n_elements());

  // should be able to hold an ElemID
  static_assert(sizeof(Mesquite2::Mesh::ElementHandle) >= sizeof(int),
                "Mesquite ElementHandle too small to hold TalyFEM element index");

  for (int i = 0; i < elements.size(); i++) {
    elements.at(i) = elem_id_to_handle(i);
  }
}

void TalyMesqMesh::get_all_vertices(std::vector<Mesquite2::Mesh::VertexHandle> &vertices,
                                    Mesquite2::MsqError &err) {
  vertices.resize(grid_->n_nodes());

  // should be able to hold a LocalNodeID
  static_assert(sizeof(Mesquite2::Mesh::VertexHandle) >= sizeof(LocalNodeID),
                "Mesquite VertexHandle too small to hold TalyFEM node index - 32-bit system w/ 64-bit PETSc indices?");

  for (LocalNodeID i = 0; i < vertices.size(); i++) {
    vertices.at(i) = node_id_to_handle(i);
  }
}

void TalyMesqMesh::vertices_get_fixed_flag(Mesquite2::Mesh::VertexHandle const* vert_array,
                                           std::vector<bool> &fixed_flag_array, size_t num_vtx,
                                           Mesquite2::MsqError &err) {
  fixed_flag_array.resize(num_vtx);
  for (size_t i = 0; i < num_vtx; i++) {
    LocalNodeID node_idx = vert_handle_to_id(vert_array[i]);
    bool fixed = ((grid_->GetNode(node_idx)->indicators() & fixed_mask_) != 0);
    fixed_flag_array.at(i) = fixed;
  }
}

void TalyMesqMesh::vertices_get_slaved_flag(Mesquite2::Mesh::VertexHandle const* vert_array,
                                            std::vector<bool> &slaved_flag_array, size_t num_vtx,
                                            Mesquite2::MsqError &err) {
  slaved_flag_array.resize(num_vtx);
  for (size_t i = 0; i < num_vtx; i++) {
    LocalNodeID node_idx = vert_handle_to_id(vert_array[i]);
    slaved_flag_array.at(i) = (bool) slaved_.at(node_idx);
  }
}

void TalyMesqMesh::vertices_get_coordinates(Mesquite2::Mesh::VertexHandle const* vert_array,
                                            Mesquite2::MsqVertex* coordinates, size_t num_vtx,
                                            Mesquite2::MsqError &err) {
  for (size_t i = 0; i < num_vtx; i++) {
    LocalNodeID node_idx = vert_handle_to_id(vert_array[i]);
    const ZEROPTV &node_coords = grid_->GetNode(node_idx)->location();
    coordinates[i].set(node_coords.data());
  }
}

void TalyMesqMesh::vertex_set_coordinates(Mesquite2::Mesh::VertexHandle vertex,
                                          const Mesquite2::Vector3D &coordinates,
                                          Mesquite2::MsqError &err) {
  LocalNodeID node_idx = vert_handle_to_id(vertex);
  ZEROPTV ptv = ZEROPTV(coordinates.x(), coordinates.y(), coordinates.z());
  grid_->GetNode(node_idx)->location() = ptv;
}

void TalyMesqMesh::vertex_set_byte(Mesquite2::Mesh::VertexHandle vertex, unsigned char byte,
                                   Mesquite2::MsqError &err) {
  LocalNodeID node_idx = vert_handle_to_id(vertex);
  mesquite_vert_byte_.at(node_idx) = byte;
}

void TalyMesqMesh::vertices_set_byte(Mesquite2::Mesh::VertexHandle const* vert_array,
                                     const unsigned char* byte_array, size_t array_size,
                                     Mesquite2::MsqError &err) {
  for (size_t i = 0; i < array_size; i++) {
    LocalNodeID node_idx = vert_handle_to_id(vert_array[i]);
    mesquite_vert_byte_.at(node_idx) = byte_array[i];
  }
}

void TalyMesqMesh::vertex_get_byte(Mesquite2::Mesh::VertexHandle const vertex, unsigned char* byte,
                                   Mesquite2::MsqError &err) {
  LocalNodeID node_idx = vert_handle_to_id(vertex);
  *byte = mesquite_vert_byte_.at(node_idx);
}

void TalyMesqMesh::vertices_get_byte(Mesquite2::Mesh::VertexHandle const* vertex, unsigned char* byte_array,
                                     size_t array_size, Mesquite2::MsqError &err) {
  for (size_t i = 0; i < array_size; i++) {
    LocalNodeID node_idx = vert_handle_to_id(vertex[i]);
    byte_array[i] = mesquite_vert_byte_.at(node_idx);
  }
}

void TalyMesqMesh::vertices_get_attached_elements(Mesquite2::Mesh::VertexHandle const* vertex_array,
                                                  size_t num_vertex,
                                                  std::vector<Mesquite2::Mesh::ElementHandle> &elements,
                                                  std::vector<size_t> &offsets, Mesquite2::MsqError &err) {
  offsets.clear();
  offsets.reserve(num_vertex + 1);

  elements.clear();
  elements.reserve(num_vertex * 2);  // a minimum guess
  for (size_t i = 0; i < num_vertex; i++) {
    offsets.push_back(elements.size());

    LocalNodeID node_idx = vert_handle_to_id(vertex_array[i]);
    const auto& this_elems = elems_containing_node_.at(node_idx);

    elements.reserve(elements.size() + this_elems.size());
    for (int j = 0; j < this_elems.size(); j++) {
      elements.push_back(elem_id_to_handle(this_elems.at(j)));
    }
  }
  offsets.push_back(elements.size());
}

void TalyMesqMesh::elements_get_attached_vertices(Mesquite2::Mesh::ElementHandle const* elem_handles,
                                                  size_t num_elems,
                                                  std::vector<Mesquite2::Mesh::VertexHandle> &vert_handles,
                                                  std::vector<size_t> &offsets, Mesquite2::MsqError &err) {
  offsets.clear();
  offsets.reserve(num_elems + 1);

  vert_handles.clear();
  vert_handles.reserve(num_elems * 2);  // a guess

  for (size_t i = 0; i < num_elems; i++) {
    offsets.push_back(vert_handles.size());

    int elem_idx = elem_handle_to_id(elem_handles[i]);
    const ELEM* elem = grid_->GetElm(elem_idx);

    vert_handles.reserve(vert_handles.size() + elem->n_nodes());
    for (ElemNodeID j = 0; j < elem->n_nodes(); j++) {
      LocalNodeID node_idx = elem->ElemToLocalNodeID(j);
      vert_handles.push_back(node_id_to_handle(node_idx));
    }
  }
  offsets.push_back(vert_handles.size());
}

Mesquite2::EntityTopology taly_elemtype_to_mesquite(ElemType type) {
  switch (type) {
    case kElem3dHexahedral: return Mesquite2::HEXAHEDRON;
    case kElem2dBox: return Mesquite2::QUADRILATERAL;
    case kElem3dTetrahedral: return Mesquite2::TETRAHEDRON;
    case kElem2dTriangle: return Mesquite2::TRIANGLE;

    case kElem1d:
    default: throw NotImplementedException() << "Unhandled taly to mesquite type " << type;
  }
}

void TalyMesqMesh::elements_get_topologies(Mesquite2::Mesh::ElementHandle const* element_handle_array,
                                           Mesquite2::EntityTopology* element_topologies,
                                           size_t num_elements, Mesquite2::MsqError &err) {
  for (size_t i = 0; i < num_elements; i++) {
    int elem_id = elem_handle_to_id(element_handle_array[i]);
    const ELEM* elem = grid_->GetElm(elem_id);
    element_topologies[i] = taly_elemtype_to_mesquite(elem->elmType());
  }
}

Mesquite2::TagHandle
TalyMesqMesh::tag_create(const std::string &tag_name, Mesquite2::Mesh::TagType type, unsigned length,
                                     const void* default_value, Mesquite2::MsqError &err) {
  throw NotImplementedException() << "tag functions not implemented";
}

void TalyMesqMesh::tag_delete(Mesquite2::TagHandle handle, Mesquite2::MsqError &err) {
  throw NotImplementedException() << "tag functions not implemented";
}

Mesquite2::TagHandle TalyMesqMesh::tag_get(const std::string &name, Mesquite2::MsqError &err) {
  throw NotImplementedException() << "tag functions not implemented";
}

void TalyMesqMesh::tag_properties(Mesquite2::TagHandle handle, std::string &name_out,
                                              Mesquite2::Mesh::TagType &type_out, unsigned &length_out,
                                              Mesquite2::MsqError &err) {
  throw NotImplementedException() << "tag functions not implemented";
}

void TalyMesqMesh::tag_set_element_data(Mesquite2::TagHandle handle, size_t num_elems,
                                                    Mesquite2::Mesh::ElementHandle const* elem_array,
                                                    const void* tag_data, Mesquite2::MsqError &err) {
  throw NotImplementedException() << "tag functions not implemented";
}

void TalyMesqMesh::tag_set_vertex_data(Mesquite2::TagHandle handle, size_t num_elems,
                                                   Mesquite2::Mesh::VertexHandle const* node_array,
                                                   const void* tag_data, Mesquite2::MsqError &err) {
  throw NotImplementedException() << "tag functions not implemented";
}

void TalyMesqMesh::tag_get_element_data(Mesquite2::TagHandle handle, size_t num_elems,
                                                    Mesquite2::Mesh::ElementHandle const* elem_array, void* tag_data,
                                                    Mesquite2::MsqError &err) {
  throw NotImplementedException() << "tag functions not implemented";
}

void TalyMesqMesh::tag_get_vertex_data(Mesquite2::TagHandle handle, size_t num_elems,
                                                   Mesquite2::Mesh::VertexHandle const* node_array, void* tag_data,
                                                   Mesquite2::MsqError &err) {
  throw NotImplementedException() << "tag functions not implemented";
}

TalyMesqMesh::TalyMesqMesh(GRID* grid) {
  grid_ = grid;
  assert(grid_ != NULL);

  fixed_mask_ = static_cast<NodeIndicator>(~0x0);

  redim();
}

/**
 * Returns the index of the first "higher order" node in a particular element type.
 * This is the same as the number of nodes in the lowest-order element.
 * Used to decide which nodes are considered "slaves" to Mesquite.
 * @param type element type
 * @return index of the first "higher order" node for this element type
 *         (regardless of if this element is higher order or not)
 */
ElemNodeID first_higher_order_idx(ElemType type) {
  switch (type) {
    case kElem3dHexahedral: return 8;
    case kElem2dBox: return 4;
    case kElem1d: return 2;
    case kElem3dTetrahedral: return 4;
    case kElem2dTriangle: return 3;
    default: throw NotImplementedException() << "first higher order node for elemType " << type << " not implemented";
  }
}

void TalyMesqMesh::redim() {
  mesquite_vert_byte_.resize(grid_->n_nodes(), 0);
  slaved_.resize(grid_->n_nodes(), false);
  elems_containing_node_.resize(grid_->n_nodes());

  // loop over all elements
  for (int elm_id = 0; elm_id < grid_->n_elements(); elm_id++) {
    const ELEM* elem = grid_->GetElm(elm_id);

    // fill in slaved_
    // loop over the higher-order nodes in the element and mark them as slaved
    ElemNodeID start = first_higher_order_idx(elem->elmType());
    for (ElemNodeID i = start; i < elem->n_nodes(); i++) {
      LocalNodeID node_idx = elem->ElemToLocalNodeID(i);
      slaved_.at(node_idx) = true;
    }

    // fill in elems_containing_vert_
    // loop over all the nodes in this element and add this element to their elems_containing list
    for (ElemNodeID i = 0; i < elem->n_nodes(); i++) {
      LocalNodeID node_idx = elem->ElemToLocalNodeID(i);
      elems_containing_node_.at(node_idx).push_back(elm_id);
    }
  }
}

}  // namespace TALYFEMLIB

#endif