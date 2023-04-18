#pragma once

#ifdef ENABLE_MESQUITE

#include <Mesquite_all_headers.hpp>
#include <talyfem/grid/grid_types/grid.h>

namespace TALYFEMLIB {

/**
 * Adapter for Mesquite mesh operations on a TalyFEM GRID object.
 * This allows us to perform Mesquite mesh optimizations directly on the TalyFEM GRID data,
 * with no copying.
 *
 * The Mesquite ElementHandle (internally a void* pointer) is treated as an element ID (int).
 * The Mesquite VertexHandle (also internally a void*) is treated as a LocalNodeID
 * (aka PetscInt, which is either a 32-bit or 64-bit int depending on how PETSc was built).
 * No extra memory is allocated for indices - we just cast the indices to void* and back when necessary.
 */
class TalyMesqMesh : public Mesquite::Mesh {
 public:
  int get_geometric_dimension(Mesquite2::MsqError &err) override;

  void get_all_elements(std::vector<ElementHandle> &elements, Mesquite2::MsqError &err) override;

  void get_all_vertices(std::vector<VertexHandle> &vertices, Mesquite2::MsqError &err) override;

  void vertices_get_fixed_flag(VertexHandle const* vert_array, std::vector<bool> &fixed_flag_array, size_t num_vtx,
                               Mesquite2::MsqError &err) override;

  void vertices_get_slaved_flag(VertexHandle const* vert_array, std::vector<bool> &slaved_flag_array, size_t num_vtx,
                                Mesquite2::MsqError &err) override;

  void vertices_get_coordinates(VertexHandle const* vert_array, Mesquite2::MsqVertex* coordinates, size_t num_vtx,
                                Mesquite2::MsqError &err) override;

  void vertex_set_coordinates(VertexHandle vertex, const Mesquite2::Vector3D &coordinates,
                              Mesquite2::MsqError &err) override;

  void vertex_set_byte(VertexHandle vertex, unsigned char byte, Mesquite2::MsqError &err) override;

  void vertices_set_byte(const VertexHandle* vert_array, const unsigned char* byte_array, size_t array_size,
                         Mesquite2::MsqError &err) override;

  void vertex_get_byte(const VertexHandle vertex, unsigned char* byte, Mesquite2::MsqError &err) override;

  void vertices_get_byte(const VertexHandle* vertex, unsigned char* byte_array, size_t array_size,
                         Mesquite2::MsqError &err) override;

  void vertices_get_attached_elements(const VertexHandle* vertex_array, size_t num_vertex,
                                      std::vector<ElementHandle> &elements, std::vector<size_t> &offsets,
                                      Mesquite2::MsqError &err) override;

  void elements_get_attached_vertices(const ElementHandle* elem_handles, size_t num_elems,
                                      std::vector<VertexHandle> &vert_handles, std::vector<size_t> &offsets,
                                      Mesquite2::MsqError &err) override;

  void elements_get_topologies(const ElementHandle* element_handle_array, Mesquite2::EntityTopology* element_topologies,
                               size_t num_elements, Mesquite2::MsqError &err) override;

  Mesquite2::TagHandle tag_create(const std::string &tag_name, TagType type, unsigned length, const void* default_value,
                                  Mesquite2::MsqError &err) override;

  void tag_delete(Mesquite2::TagHandle handle, Mesquite2::MsqError &err) override;

  Mesquite2::TagHandle tag_get(const std::string &name, Mesquite2::MsqError &err) override;

  void tag_properties(Mesquite2::TagHandle handle, std::string &name_out, TagType &type_out, unsigned &length_out,
                      Mesquite2::MsqError &err) override;

  void tag_set_element_data(Mesquite2::TagHandle handle, size_t num_elems, const ElementHandle* elem_array,
                            const void* tag_data, Mesquite2::MsqError &err) override;

  void tag_set_vertex_data(Mesquite2::TagHandle handle, size_t num_elems, const VertexHandle* node_array,
                           const void* tag_data, Mesquite2::MsqError &err) override;

  void
  tag_get_element_data(Mesquite2::TagHandle handle, size_t num_elems, const ElementHandle* elem_array, void* tag_data,
                       Mesquite2::MsqError &err) override;

  void
  tag_get_vertex_data(Mesquite2::TagHandle handle, size_t num_elems, const VertexHandle* node_array, void* tag_data,
                      Mesquite2::MsqError &err) override;

  void release_entity_handles(const EntityHandle* handle_array, size_t num_handles, Mesquite2::MsqError &err) override {}

  void release() override {}

  // ---

  explicit TalyMesqMesh(GRID* grid);
  virtual ~TalyMesqMesh() = default;

  void redim();

  /**
   * Change the bitmask that specifies which nodes are considered "fixed"
   * (not affected by the Mesquite mesh optimization routines).
   * Default is any node indicator means the node is fixed.
   * Use something like "INDICATOR_NUM(1) | INDICATOR_NUM(2) | ..."
   * to build a list of indicators.
   * @param boundary_flag_mask if (node & mask) != 0, the node is considered fixed
   */
  inline void set_fixed_mask(NodeIndicator boundary_flag_mask) {
    fixed_mask_ = boundary_flag_mask;
  }

 private:
  inline static LocalNodeID vert_handle_to_id(Mesquite2::Mesh::VertexHandle v) {
    return *((LocalNodeID*) &v);
  }

  inline static Mesquite2::Mesh::VertexHandle node_id_to_handle(LocalNodeID id) {
    return reinterpret_cast<Mesquite2::Mesh::VertexHandle>(id);
    //return (Mesquite2::Mesh::VertexHandle) ((intptr_t) id);
  }

  inline static Mesquite2::Mesh::ElementHandle elem_id_to_handle(int id) {
    //return (Mesquite2::Mesh::ElementHandle) ((intptr_t) id);
    return reinterpret_cast<Mesquite2::Mesh::ElementHandle>(id);
  }

  inline static int elem_handle_to_id(Mesquite2::Mesh::ElementHandle e) {
    return *((int*) &e);
  }

  GRID* grid_;  ///< grid we are using

  NodeIndicator fixed_mask_;  ///< bitmask that specifies which node indicator flags mean the node should be "fixed"
                              ///< (not affected by mesh optimization routines) - default is any (~0x0)
  std::vector<bool> slaved_;  ///< slaved_[LocalNodeID] = true means the given node ID is a "higher order" vertex,
                              ///< part of a higher-order element that can be recomputed after corners are moved
                              ///< taly doesn't track this information, so it is filled in at redim time here

  std::vector<unsigned char> mesquite_vert_byte_;  ///< per-vertex byte of Mesquite data

  // not the most memory efficient choice, but simplest
  std::vector< std::vector<int> > elems_containing_node_;  ///< each entry corresponds to a node, and contains
                                                           ///< the list of elements that contain that node
                                                           ///< used in vertices_get_attached_elements
};

}

#endif