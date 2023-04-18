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
#ifndef GRID_GRIDFIELD_H_
#define GRID_GRIDFIELD_H_

#include <talyfem/domain_decomposition/cross_domain_communicator.h>
#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/femelm.h>
#include <talyfem/grid/nodeid_types.h>
#include <talyfem/grid/shareinfo.h>


namespace TALYFEMLIB {

DEFINE_HAS_MEMBER(has_UpdateDataStructures, UpdateDataStructures)

/**
 * Grid field class to manipulate data stored at node points.
 */
template<class NodeData>
class GridField {
 public:
  GridField()
      : p_grid_(NULL),
        append_output_(false),
        communicator_(),
        p_node_data_array_(NULL),
        node_count_(0) {
  }

  virtual ~GridField() {
    CleanNodeData();
  }

  /**
   * Free memory used to store data at nodal points
   */
  virtual void CleanNodeData() {
    if (p_node_data_array_) {
      for (LocalNodeID i = 0; i < node_count_; i++) {
        delete p_node_data_array_[i];
      }
      delete [] p_node_data_array_;
    }
    p_node_data_array_ = NULL;
  }

  /**
   * Set the grid for this object
   *
   * This can also set up the cross domain communication if domain
   * decomposition is being used.
   * Note that the base implementation of this wipes the node data
   *
   * @param p_grid pointer to grid to use for this object
   * @param set_up_comm whether to set up the cross domain communicator
   */
  virtual void redimGrid(GRID* p_grid, bool set_up_comm = true) {
    CleanNodeData();
    p_grid_ = p_grid;
    if ((p_grid_->parallel_type_ == kWithDomainDecomp) && set_up_comm) {
      communicator_.Initialize(p_grid);
    }
  }

  /**
   * Set up the data structures to store node data in this object
   *
   * This creates the array of node data objects
   *
   * @param p_grid grid to use (if NULL, use object's grid)
   */
  virtual void redimNodeData(GRID* p_grid = NULL) {
    if (p_grid) {  // if grid is given, use it
      p_grid_ = p_grid;
    } else {  // if no grid is given, store existing grid in the passed pointer
      p_grid = p_grid_;
    }
    CleanNodeData();
    node_count_ = p_grid_->n_nodes();

    p_node_data_array_ = new NodeData*[node_count_];
    for (LocalNodeID i = 0; i < node_count_; i++) {
      p_node_data_array_[i] = new NodeData();
    }
  }

  /**
   * Set up the element data structures to store element data in this object
   *
   * This creates the array of element data objects. Currently, this does
   * not do anything in the base implementation. The default values of -1
   * for parameters have no significiance.
   *
   * @param p_grid grid to use (if NULL, use object's grid)
   * @param gpno number of gauss points ?
   * @param surfaceno number of surfaces ?
   * @param surfacegpno number of surface gauss points ?
   */
  virtual void redimElmData(GRID* p_grid = NULL, int gpno = -1,
                            int surfaceno = -1, int surfacegpno = -1) {
    if (p_grid) {  // if grid is given, use it
      p_grid_ = p_grid;
    } else {  // if no grid is given, store existing grid in the passed pointer
      p_grid = p_grid_;
    }
  }

  /**
   * Returns the node data object for the given node id
   *
   * @param node_id index of node whose node data is desired
   * @return node data object for the given node id
   */
  NodeData& GetNodeData(LocalNodeID node_id) {
    assert(node_id >= 0 && node_id < p_grid_->n_nodes());
    return *(p_node_data_array_[node_id]);
  }

  /**
   * Returns the node data object for the given node id
   *
   * @param node_id index of node whose node data is desired
   * @return node data object for the given node id
   */
  const NodeData& GetNodeData(LocalNodeID node_id) const {
    assert(node_id >= 0 && node_id < p_grid_->n_nodes());
    return *(p_node_data_array_[node_id]);
  }

  /**
   * Returns the FE averaged value of the given data index
   *
   * @param fe finite element to average value over
   * @param index index of data value (in node data object) to average
   * @return value averaged over the element
   */
  double valueFEM(const FEMElm& fe, int index) const {
    double sum = 0.0;

    if (fe.basis_function() == BASIS_HERMITE) {
      // Hermite is a special case - some basis functions are derivatives.
      // This assumes the node data is stored such that u = index
      // and du_x = index + 1, du_y = index + 2, etc.
      // 1D: -, x
      // 2D: -, x, y, xy
      // 3D: -, x, y, z, xy, xz, yz, xyz
      const int nbf = fe.nbf();
      const int nbf_per_node = fe.nbf_per_node();
      for (int a = 0; a < nbf; a++) {
        const int node_idx = fe.elem()->ElemToLocalNodeID(a / nbf_per_node);
        const int deriv = a % nbf_per_node;
        sum += fe.N(a) * GetNodeData(node_idx).value(index + deriv);
      }
    } else {
      const int nbf = fe.nbf();
      for (ElemNodeID a = 0; a < nbf; a++) {
        const int node_index = fe.elem()->ElemToLocalNodeID(a);
        sum += fe.N(a) * GetNodeData(node_index).value(index);
      }
    }
    return sum;
  }

  /**
   * Returns the FE averaged derivative of value of the given data index
   *
   * @param fe finite element to average value over
   * @param index index of data value (in node data object) to average
   * @param dir specifies which direction the derivative is for (0=x, 1=y, 2=z)
   * @return derivative value averaged over the element
   */
  double valueDerivativeFEM(const FEMElm& fe, int index, int dir) const {
    double sum = 0.0;

    if (fe.basis_function() == BASIS_HERMITE) {
      const int nbf = fe.nbf();
      const int nbf_per_node = fe.nbf_per_node();
      for (int a = 0; a < nbf; a++) {
        const int node_idx = fe.elem()->ElemToLocalNodeID(a / nbf_per_node);
        const int deriv = a % nbf_per_node;
        sum += fe.dN(a, dir) * GetNodeData(node_idx).value(index + deriv);
      }
    } else {
      const int nbf = fe.nbf();
      for (ElemNodeID a = 0; a < nbf; a++) {
        const int node_index = fe.elem()->ElemToLocalNodeID(a);
        sum += fe.dN(a, dir) * GetNodeData(node_index).value(index);
      }
    }
    return sum;
  }

  /**
   * Returns the FE averaged 2nd derivative of value of the given data index.
   * @param fe finite element to average value over
   * @param index index of data value (in node data object) to average
   * @param dir1 first direction in which to take derivative
   * @param dir2 second direction in which to take derivative
   * @return 2nd derivative value (wrt dir1 and dir2) averaged over the element
   */
  double value2DerivativeFEM(const FEMElm& fe, int index, int dir1, int dir2) const {
    double sum = 0.0;

    if (fe.basis_function() == BASIS_HERMITE) {
      throw NotImplementedException() << "value2DerivativeFEM not "
                                         "implemented for hermite";
    }

    const int nbf = fe.nbf();
    for (ElemNodeID a = 0; a < nbf; a++) {
      const int node_index = fe.elem()->ElemToLocalNodeID(a);
      sum += fe.d2N(a, dir1, dir2) * GetNodeData(node_index).value(index);
    }
    return sum;
  }

  /**
   * Updates the node data structures
   *
   * The meaning of updating the data structure is dependent on the
   * UpdateDataStructures function in the NodeData class
   */
  template <bool Exists = has_UpdateDataStructures<NodeData>::value>
  typename std::enable_if<Exists, void>::type UpdateDataStructures() {
    const int n_nodes = p_grid_->n_nodes();
    for (LocalNodeID i = 0; i < n_nodes; i++) {
      GetNodeData(i).UpdateDataStructures();
    }
  }

  /**
   * Copies data from one place in the node data to another in the same node
   *
   * This is used, for example, when overwriting the values from a previous
   * time step with values from the current timestep in preparating for the
   * next calculation.
   *
   * @param source_index where in the node data to move the values from
   * @param destination_index where in the node data to move the values to
   */
  virtual void CopyNodeDataFromTo(int source_index, int destination_index) {
    const int n_nodes = p_grid_->n_nodes();
    for (LocalNodeID i = 0; i < n_nodes; i++) {
      NodeData& node_data = GetNodeData(i);
      node_data.value(destination_index) = node_data.value(source_index);
    }
  }

  /**
   * Stores data from given index of all NodeData objects into the given array.
   *
   * @param array where to put the data (must be length = n_nodes)
   * @param source_index where to get data from in NodeData objects
   */
  void NodeDataToArray(double *array, int source_index) const {
    const int n_nodes = p_grid_->n_nodes();
    for (LocalNodeID i = 0; i < n_nodes; i++) {
      NodeData& node_data = GetNodeData(i);
      array[i] = node_data.value(source_index);
    }
  }

  /**
   * Puts data from given array into specified location in all NodeData objects.
   *
   * @param array array with data to store (must be length = n_nodes)
   * @param stride distance between values in array
   * @param destination_index where to put data in NodeData objects
   */
  void NodeDataFromArray(const double *array, int destination_index,
                         int stride = 1) {
    const int n_nodes = p_grid_->n_nodes();
    for (LocalNodeID i = 0; i < n_nodes; i++) {
      NodeData& node_data = GetNodeData(i);
      node_data.value(destination_index) = array[i*stride];
    }
  }

  /**
   * Syncs nodal values across domain boundaries
   *
   * This performs the synchonization by copying the value from the process
   * that "owns" the node to every other process that has the node.
   *
   * @param index which value in the node data to sync
   */
  void Communicate(int index) {
    // There is nothing to do if there is no domain decomposition
    if (p_grid_->parallel_type_ == kNoDomainDecomp) { return; }

    // set up the data to send
    const int n_send_nodes = communicator_.n_send_nodes();
    const int *send_nodes = communicator_.send_node_indices();
    double *send_data = communicator_.send_data();
    for (int i = 0; i < n_send_nodes; i++) {
      int node_id = send_nodes[i];
      send_data[i] = GetNodeData(node_id).value(index);
    }

    communicator_.SendRecvData();  // exchange data

    // copy the received data
    const int n_recv_nodes = communicator_.n_recv_nodes();
    const PetscInt* recv_nodes = communicator_.recv_node_indices();
    const double *recv_data = communicator_.recv_data();
    for (int i = 0; i < n_recv_nodes; i++) {
      int node_id = recv_nodes[i];
      GetNodeData(node_id).value(index) = recv_data[i];
    }
  }

  /**
   * Returns a node data value averaged over the entire system
   *
   * This correctly accounts for nodes that are duplicated across node
   * boundaries. It does *NOT* account for nodes duplicated via periodic
   * boundaries. When using periodic bounds, the average will be biased.
   *
   * @param index index of the node data item to average
   * @return the node data value averaged over the entire system
   */
  double CalcNodalAverage(int index) {
    double sum = CalcNodalSum(index);
    return sum / static_cast<double>(p_grid_->n_nodes());
  }

  /**
   * Returns the sum of node data value over the entire system
   *
   * This correctly accounts for nodes that are duplicated across node
   * boundaries. It does *NOT* account for nodes duplicated via periodic
   * boundaries. When using periodic bounds, the sum will double count
   * periodic values.
   *
   * @param index index of the node data item to sum
   * @return sum of node data value over the entire system
   */
  double CalcNodalSum(int index) {
    double sum = 0.0;
    if (p_grid_->parallel_type_ == kNoDomainDecomp) {
      // simply loop over all nodes and sum the values
      const int n_nodes = p_grid_->n_nodes();
      for (int i_node = 0; i_node < n_nodes; i_node++) {
        const double nodal_value = GetNodeData(i_node).value(index);
        sum += nodal_value;
      }
    } else {
      // loop over all nodes, add the value to the sum. If the node is shared
      // among multiple processes, scale the value by 1/(n_sharing_processes)
      const int n_nodes = p_grid_->n_nodes();
      for (int i_node = 0; i_node < n_nodes; i_node++) {
        const int n_processes = p_grid_->node_belong_(i_node).GetShareCount();
        const double nodal_value = GetNodeData(i_node).value(index);
        sum += nodal_value / static_cast<double>(n_processes);
      }
      // now sum the values over the whole system
      MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    return sum;
  }

  /**
   * Calculates the inner product of 2 nodal variables.
   * @param indexx first variable index
   * @param y second variable's GridField
   * @param indexy second variable index
   * @param basis_rel_order order of integration
   * @returns inner product of this.GetNodeData().value(x)
              and y.GetNodeData().value(y)
   */
  virtual double InnerProduct(int indexx, const GridField<NodeData>& y,
                              int indexy,
                              int basis_rel_order = 0) const {
    const GridField<NodeData>& x = *this;
    double sum = 0.0;
    const int n_elements = p_grid_->n_elements();

    FEMElm fe(p_grid_, BASIS_FIRST_DERIVATIVE);
    for (int elm_id = 0; elm_id < n_elements; elm_id++) {
      fe.refill(elm_id, basis_rel_order);
      while (fe.next_itg_pt()) {
        const double thetax = x.valueFEM(fe, indexx);
        const double thetay = y.valueFEM(fe, indexy);
        const double w = fe.detJxW();
        sum += thetax * thetay * w;
      }
    }
    return sum;
  }

  /**
   * @param indexx first node variable index
   * @param y second node variable's GridField
   * @param indexy second node variable index
   * @param basis_rel_order order of integration
   * @returns the inner product of gradients between two nodal variables
   */
  virtual double InnerProductGradient(int indexx, const GridField<NodeData>& y,
                                      int indexy,
                                      int basis_rel_order = 0) const {
    const GridField<NodeData>& x = *this;
    double sum = 0.0;
    const int nsd = p_grid_->nsd();
    const int n_elements = p_grid_->n_elements();

    FEMElm fe(p_grid_, BASIS_FIRST_DERIVATIVE);
    for (int elm_id = 0; elm_id < n_elements; elm_id++) {
      fe.refill(elm_id, basis_rel_order);
      const int n_nodes = fe.elem()->n_nodes();
      while (fe.next_itg_pt()) {
        const double w = fe.detJxW();
        for (ElemNodeID a = 0; a < n_nodes; a++) {
          const LocalNodeID node_a = fe.elem()->ElemToLocalNodeID(a);
          const double phiai = x.GetNodeData(node_a).value(indexx);
          for (ElemNodeID b = 0; b < n_nodes; b++) {
            const LocalNodeID node_b = fe.elem()->ElemToLocalNodeID(b);
            const double phibj = y.GetNodeData(node_b).value(indexy);
            for (int i_dim = 0; i_dim < nsd; i_dim++) {
              sum += fe.dN(a, i_dim) * fe.dN(b, i_dim) * phiai * phibj * w;
            }
          }
        }
      }
    }
    return sum;
  }

  /**
   * Gets the local point corresponding to a global point in a given element
   *
   * @param ptv_global global point we want to find
   * @param[out] ptv_local coordinates of local point in element
   * @param elm_id id of element to find point in
   */
  void GetLocalPtv(const ZEROPTV& ptv_global, ZEROPTV& ptv_local,
                   int elm_id) const {
    p_grid_->GetLocalPtv(ptv_global, ptv_local, elm_id);
  }

  /**
   * Given a global point, find the element containing that point, convert
   * the point to isoparametric space, evaluate the basis functions using that
   * point as the integration point, then return valueFEM().
   * This is very similar to GRID::findElmAndLocPt.
   * @param ptv_global point in global coordinates
   * @param index NodeData variable index to calculate valueFEM with
   * @param elm_id element to look in (-1 for automatic)
   * @returns valueFEM evaluated at the given global point
   */
  double ValueAtPoint(const ZEROPTV& ptv_global, int index,
                      int elm_id = -1) const {
    if (elm_id < 0 || elm_id >= p_grid_->n_elements()) {
      const ELEM* elm = p_grid_->kd_tree().elm_containing_pt(ptv_global);
      if (elm == NULL) {
        throw TALYException() << "No element contains point " << ptv_global;
      }
      elm_id = elm->elm_id();
    }

    ZEROPTV ptv_local;
    GetLocalPtv(ptv_global, ptv_local, elm_id);
    return ValueAtLocPoint(ptv_local, index, elm_id);
  }

  /**
   * Given a local point, evaluate the basis function at that point
   * and return valueFEM() for the given index.
   * @param ptv_local local point
   * @param index NodeData variable index to calculate valueFEM with
   * @param elm_id element to calculate basis functions for
   * @param basis_rel_order relative order of integration (0 by default)
   * @returns valueFEM() of variable index at isoparametric point ptv_local
   */
  double ValueAtLocPoint(const ZEROPTV& ptv_local, int index,
                         int elm_id,
                         int basis_rel_order = 0) const {
    FEMElm fe(p_grid_, BASIS_POSITION | BASIS_FIRST_DERIVATIVE);
    fe.refill(elm_id, basis_rel_order);
    fe.calc_at(ptv_local);
    return valueFEM(fe, index);
  }

  GRID* p_grid_;  ///< pointer to underlying grid
  bool append_output_;  ///< whether to append output when writing grid field

 protected:
  CrossDomainCommunicator communicator_;  ///< communicator for sending data
                                          ///< across domain boundaries
  NodeData** p_node_data_array_;  ///< node data
  int node_count_;  ///< number of nodes. This is only used during destruction
                    ///< in order to allow the object to be destroyed without
                    ///< needing information from the grid object.
};

}  // namespace TALYFEMLIB

#endif  // GRID_GRIDFIELD_H_
