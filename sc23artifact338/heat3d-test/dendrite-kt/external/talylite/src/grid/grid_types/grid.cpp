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
#include <talyfem/grid/grid_types/grid.h>

#include <algorithm>  // for std::min
#include <vector>

#include <talyfem/input_data/input_data.h>
#include <talyfem/domain_decomposition/mesh_partition.h>
#include <talyfem/file_io/file_io.h>
#include <talyfem/grid/node.h>
#include <talyfem/grid/nodedata.h>
#include <talyfem/grid/segment.h>
#include <talyfem/grid/shareinfo.h>
#include <talyfem/grid/elem-types.h>  // for make_elem_of_type
#include <talyfem/grid/kdtree.h>

namespace TALYFEMLIB {

void ReOrient(GRID* pGrid, ELEM* pElm) {
  if (pElm->elmType() == kElem2dTriangle) {
    ZEROPTV PAB, PBC;
    LocalNodeID A = pElm->ElemToLocalNodeID(0);
    LocalNodeID B = pElm->ElemToLocalNodeID(1);
    LocalNodeID C = pElm->ElemToLocalNodeID(2);
    for (int i = 0; i < 2; i++) {
      PAB(i) = pGrid->GetNode(B)->getCoor(i) - pGrid->GetNode(A)->getCoor(i);
      PBC(i) = pGrid->GetNode(C)->getCoor(i) - pGrid->GetNode(B)->getCoor(i);
    }
    ZEROPTV cross;
    cross.crossProduct(PAB, PBC);
    if (cross.z() >= 0)
      return;
    LocalNodeID tmp = pElm->node_id_array(0);
    pElm->set_node_id_array(0, pElm->node_id_array(1));
    pElm->set_node_id_array(1, tmp);
  } else if (pElm->elmType() == kElem3dTetrahedral) {
    PrintStatus("Checking element: ", pElm->elm_id());
    ZEROPTV PAB, PAC, PAD;
    LocalNodeID A = pElm->ElemToLocalNodeID(0);
    LocalNodeID B = pElm->ElemToLocalNodeID(1);
    LocalNodeID C = pElm->ElemToLocalNodeID(2);
    LocalNodeID D = pElm->ElemToLocalNodeID(3);
    for (int i = 0; i < 3; i++) {
      PAB(i) = pGrid->GetNode(B)->getCoor(i) - pGrid->GetNode(A)->getCoor(i);
      PAC(i) = pGrid->GetNode(C)->getCoor(i) - pGrid->GetNode(A)->getCoor(i);
      PAD(i) = pGrid->GetNode(D)->getCoor(i) - pGrid->GetNode(A)->getCoor(i);
    }
    ZEROPTV cross;
    cross.crossProduct(PAC, PAD);

    PrintStatus("+-->Orientation: ", PAB.innerProduct(cross));

    if (PAB.innerProduct(cross) >= 0.0) {
      return;
    }
    LocalNodeID tmp = pElm->node_id_array(2);
    pElm->set_node_id_array(2, pElm->node_id_array(3));
    pElm->set_node_id_array(3, tmp);
  } else {
    return;
  }
}

GRID::GRID(int basis_function_order)
    : elm_array_(NULL),
      node_array_(NULL),
      n_faces_(0),
      parallel_type_(kNoDomainDecomp),
      with_neighbors_(false),
      n_nodes_(0),
      n_elements_(0),
      nsd_(0),
      n_owned_nodes_(-1),
      n_total_nodes_(-1),
      n_nodes_per_direction_(NULL),
      n_elems_per_direction_(NULL),
      grid_type_(kGrid3dBox),
      grid_id_(0),
      n_subgrids_(1),
      basis_order_(basis_function_order),
      kd_tree_(this) {
  set_nsd(3);
}

GRID::~GRID() {
  Cleanup();
  // these are not in cleanup() because redim calls cleanup and these
  // should remain in place during redim
  if (n_nodes_per_direction_ != NULL) { delete [] n_nodes_per_direction_; }
  if (n_elems_per_direction_ != NULL) { delete [] n_elems_per_direction_; }
}

GRID& GRID::operator=(const GRID& grid) {
  // copy
  redimArrays(grid.n_nodes(), grid.n_elements());

  for (LocalNodeID i_node = 0; i_node < n_nodes(); i_node++) {
    NODE* pNode = grid.node_array_[i_node];

    NODE* pNewNode = new NODE();
    double x = pNode->x();
    double y = pNode->y();
    double z = pNode->z();
    pNewNode->setCoor(x, y, z);

    pNewNode->setIndicators(pNode->indicators());
    pNewNode->elem_id_ = pNode->elem_id_;
    node_array_[i_node] = pNewNode;
  }

  for (int elm_id = 0; elm_id < grid.n_elements(); elm_id++) {
    ELEM* pElm = grid.elm_array_[elm_id];

    ELEM* pNewElm = make_elem_of_type(pElm->elmType());
    pNewElm->redim(pElm->n_nodes(), pElm->node_id_array());
    pNewElm->surface_indicator_ = pElm->surface_indicator_;
    pNewElm->set_elm_id(pElm->elm_id());
    elm_array_[elm_id] = pNewElm;
  }

  cared_surface_indicator_ = grid.cared_surface_indicator_;
  set_grid_type(grid.grid_type());
  set_nsd(grid.nsd());

  boundary_ = grid.boundary_;

  parallel_type_ = grid.parallel_type_;
  set_n_total_nodes(grid.n_total_nodes());
  set_grid_id(grid.grid_id());

  solution_map_ = grid.solution_map_;
  physical_map_ = grid.physical_map_;

  n_shared_nodes_ = grid.n_shared_nodes_;
  shared_nodes_ = grid.shared_nodes_;

  n_owned_nodes_ = grid.n_owned_nodes_;
  n_total_nodes_ = grid.n_total_nodes_;

  node_belong_ = grid.node_belong_;
  kd_tree_.rebuild();

  return *this;
}

GRID& GRID::AddGrid(GRID& grid, ZEROARRAY<int>& overlapNodeID,
                    ZEROARRAY<int>& orginNodeID) {
  ELEM** old_elm_array = elm_array_;
  NODE** old_node_array = node_array_;
  int oldelmno = n_elements();
  int oldnodeno = n_nodes();

  elm_array_ = NULL;
  node_array_ = NULL;
  n_elements_ = n_nodes_ = 0;

  // copy
  redimArrays(grid.n_nodes() + oldnodeno - overlapNodeID.size(),
               grid.n_elements() + oldelmno);

  ZEROARRAY<int> nodeidmap;
  nodeidmap.redim(grid.n_nodes());
  for (LocalNodeID i_node = 0; i_node < grid.n_nodes(); i_node++) {
    nodeidmap(i_node) = oldnodeno + i_node - overlapNodeID.size();
  }
  for (int i_node = 0; i_node < overlapNodeID.size(); i_node++) {
    nodeidmap(overlapNodeID(i_node)) = orginNodeID(i_node);
  }
  for (LocalNodeID i_node = 0; i_node < n_nodes(); i_node++) {
    node_array_[i_node] = NULL;
  }
  for (int i_node = 0; i_node < oldnodeno; i_node++) {
    NODE* pNewNode = old_node_array[i_node];
    node_array_[i_node] = pNewNode;
  }
  for (LocalNodeID i_node = 0; i_node < grid.n_nodes(); i_node++) {
    if (node_array_[nodeidmap(i_node)]) {
      continue;
    }
    NODE* pNode = grid.node_array_[i_node];

    NODE* pNewNode = new NODE();
    double x = pNode->x();
    double y = pNode->y();
    double z = pNode->z();
    pNewNode->setCoor(x, y, z);

    pNewNode->setIndicators(pNode->indicators());
    pNewNode->elem_id_ = pNode->elem_id_ + oldelmno;
    node_array_[nodeidmap(i_node)] = pNewNode;
  }

  for (int elm_id = 0; elm_id < oldelmno; elm_id++) {
    ELEM* pNewElm = old_elm_array[elm_id];
    elm_array_[elm_id] = pNewElm;
  }
  for (int elm_id = 0; elm_id < grid.n_elements(); elm_id++) {
    PetscInt newe = oldelmno + elm_id;
    ELEM* pElm = grid.elm_array_[elm_id];

    ELEM* pNewElm = make_elem_of_type(pElm->elmType());
    pNewElm->redim(pElm->n_nodes(), pElm->node_id_array());
    for (ElemNodeID i = 0; i < pNewElm->n_nodes(); i++) {
      pNewElm->set_node_id_array(i, nodeidmap(pElm->node_id_array(i)));
    }
    pNewElm->surface_indicator_ = pElm->surface_indicator_;
    pNewElm->set_elm_id(pElm->elm_id() + oldelmno);
    elm_array_[newe] = pNewElm;
  }

  for (int c = 0; c < grid.cared_surface_indicator_.size(); c++) {
    AddCaredSurfaceIndicator(grid.cared_surface_indicator_(c));
  }
  set_nsd(grid.nsd());
  set_grid_type(grid.grid_type());
  parallel_type_ = grid.parallel_type_;
  kd_tree_.rebuild();

  return *this;
}

bool GRID::BoNode(LocalNodeID node_id) const {
  assert(node_id >= 0 && node_id < n_nodes());
  return node_array_[node_id]->HasIndicators();
}

bool GRID::BoNodeFlags(LocalNodeID node_id,
                       NodeIndicator indicatorFlags) const {
  assert(node_id >= 0 && node_id < n_nodes());
  return node_array_[node_id]->BoNodeFlags(indicatorFlags);
}

void GRID::BoundingBox(int elm_id, ZEROPTV& min_val, ZEROPTV& max_val) const {
  for (int dir = 0; dir < nsd(); dir++) {
    min_val(dir) = 1e20;
    max_val(dir) = -1e20;
  }

  const ElemNodeID n_nodes_per_elem = GetNumNodesInElm(elm_id);
  for (ElemNodeID i = 0; i < n_nodes_per_elem; i++) {
    for (int k = 0; k < nsd(); k++) {
      double c = GetCoord(GetLocalNodeID(elm_id, i), k);
      if (c < min_val(k)) { min_val(k) = c; }
      if (c > max_val(k)) { max_val(k) = c; }
    }
  }
}

int GRID::CalcNumNodesOthers() const {
  if (n_subgrids() == 1) { return 0; }  // no other grids
  int sum = 0;
  for (LocalNodeID i_node = 0; i_node < n_nodes(); i_node++) {
    if (!(node_belong_(i_node).is_owned)) {
      sum++;
    }
  }
  return sum;
}

void GRID::CreateGrid(const InputData* pIdata) {
  if (pIdata->ifDD) {
    this->redimDD(pIdata->L, pIdata->Nelem);
    PrintStatus("Grid generated, Domain Decomposition requested");
  } else {  // no DD
    this->redim(pIdata->L, pIdata->Nelem);
    PrintStatus("Grid generated, No Domain Decomposition requested");
  }
}

void GRID::FindElmAndLocPt(const ZEROPTV& ptvg, int &elm_id, ZEROPTV& ptvl,
                           bool search) const {
  if (elm_id < 0 || elm_id >= n_elements() || search) {
    const ELEM* elm = kd_tree_.elm_containing_pt(ptvg);
    if (elm == NULL) {
      throw TALYException() << "No element contains point '" << ptvg << "'";
    }
    elm_id = elm->elm_id();
  }

  GetLocalPtv(ptvg, ptvl, elm_id);
}

void GRID::FindElmAndLocPt(const ZEROPTV& ptvg, int &elm_id, ZEROPTV& ptvl,
                           const ZEROARRAY<int>& arrayElms) const {
  bool found_elm = false;
  for (int i = 0; i < arrayElms.size(); i++) {
    elm_id = arrayElms(i);
    if (GetElm(elm_id)->IsInnerPoint(this, ptvg)) {
      found_elm = true;
      break;
    }
  }

  if (!found_elm) {
    throw TALYException() << "No element contains point '" << ptvg << "' from "
                             "the given list (of " << arrayElms.size() <<
                             " elements)";
  }

  GetLocalPtv(ptvg, ptvl, elm_id);
}

void GRID::GenElmSurfaceIndicator() {
  for (int i = 0; i < n_elements(); i++) {
    elm_array_[i]->GenSurfaceIndicator(this, cared_surface_indicator_);
  }
}

double GRID::GetCoord(LocalNodeID node_id, int dir) const {
  assert(node_id >= 0 && node_id < n_nodes());
  return node_array_[node_id]->getCoor(dir);
}

void GRID::GetCoord(ZEROPTV& point, LocalNodeID node_id) const {
  for (int dir = 0; dir < nsd(); dir++) {
    point(dir) = GetCoord(node_id, dir);
  }
}

ELEM* GRID::GetElm(int e) {
  assert(e >= 0 && e < n_elements_);
  return elm_array_[e];
}

const ELEM* GRID::GetElm(int e) const {
  assert(e >= 0 && e < n_elements_);
  return elm_array_[e];
}

LocalNodeID GRID::GetLocalNodeID(int elm_id, ElemNodeID node_id) const {
  return GetElm(elm_id)->ElemToLocalNodeID(node_id);
}

void GRID::GetLocalPtv(const ZEROPTV &ptvg, ZEROPTV &ptvl, int elm_id) const {

  FEMElm fe(this, BASIS_POSITION | BASIS_FIRST_DERIVATIVE);

  // Always use linear basis functions, even for quadratic elements.
  // Nonlinear bases are not guaranteed to converge in the numerical
  // approximation below, and higher order BFs are more expensive to compute.
  // We don't support curved elements, and higher-order elements can be safely reduced
  // to linear (since the node points are in the same order), so this is an acceptable optimization.
  fe.refill(elm_id, BASIS_LINEAR, 0);

  // use (0, 0, 0) as our initial guess
  ptvl = ZEROPTV(0.0, 0.0, 0.0);

  // maximum number of iterations (after which we give up & throw an exception)
  static const int MAX_STEPS = 100;

  // upper limit on L2 error for ptvl
  static const double MIN_ERROR = 1e-7;

  int n_steps = 0;
  double err = 0.0;
  while (n_steps < MAX_STEPS) {

    // calculate the basis function values at our current guess
    fe.calc_at(ptvl);

    // calculate (L2 error)^2 for this guess, and at the same time build
    // the vector from our guess (mapped to physical space) to ptvg
    ZEROPTV DX;
    err = 0.0;
    for (int i = 0; i < fe.nsd(); i++) {
      DX(i) = (ptvg(i) - fe.position()(i));
      err += DX(i) * DX(i);
    }

    // if the error on our guess is acceptable, return it
    // (we use MIN_ERROR^2 here to avoid calculating sqrt(err))
    if (err < MIN_ERROR*MIN_ERROR)
      return;

    // otherwise, move our guess by the inverse of the Jacobian times
    // our delta x vector (calculated above)
    double jaccInv = 1.0 / fe.jacc();
    for (int i = 0; i < fe.nsd(); i++) {
      for (int j = 0; j < fe.nsd(); j++) {
        ptvl(i) += fe.cof(j, i) * DX(j) * jaccInv;
      }
    }

    n_steps++;
  }

  // if we get here, our loop above failed to converge
  throw TALYException() << "GetLocalPtv did not converge after "
                        << MAX_STEPS << " iterations (final error: "
                        << sqrt(err) << ")";
}

int GRID::GetNumOwnedNodes() {
  if (parallel_type_ == kNoDomainDecomp) {
    n_owned_nodes_ = 0;
    for (LocalNodeID n = 0; n < this->n_nodes(); n++) {
      int elm_id = node_array_[n]->elem_id_;
      if (IsMyElement(elm_id)) {
        n_owned_nodes_++;
      }
    }
  }
  return n_owned_nodes_;
}

NODE* GRID::GetNode(LocalNodeID n) {
  assert(n >= 0 && n < n_nodes());
  return node_array_[n];
}

const NODE* GRID::GetNode(LocalNodeID n) const {
  assert(n >= 0 && n < n_nodes());
  return node_array_[n];
}

void GRID::GetNodeCoords(ZEROARRAY<double>& position) const {
  assert(position.size() >= n_nodes() * nsd());
  for (LocalNodeID i_node = 0; i_node < n_nodes(); i_node++) {
    for (int dir = 0; dir < nsd(); dir++) {
      double coordinate_value = node_array_[i_node]->getCoor(dir);
      int idx = i_node * nsd() + dir;
      position(idx) = coordinate_value;
    }
  }
}

bool GRID::IsInnerPoint(const ZEROPTV& point) const {
  // this is very inefficient, likely could be substantially improved
  for (int elm_id = 0; elm_id < n_elements(); elm_id++) {
    if (GetElm(elm_id)->IsInnerPoint(this, point)) {
      return 1;
    }
  }
  return 0;
}

bool GRID::IsMyElement(int elm_id) const {
  if (parallel_type_ == kNoDomainDecomp) {
    int mpi_rank = GetMPIRank();  // TODO: replace this with grid_id?
    if (GetGridIdOfElementOwner(elm_id) == mpi_rank) {
      return true;
    }
    return false;
  }
  return true;
}

void GRID::LoadFromALBERTAFileDD(const char* filename) {
  CMeshPartition pmesh;
  pmesh.LoadFromALBERTAFile(filename, 0);
  pmesh.TransferToGrid(this);
  // ~ pmesh.GetISNodeCopies(&isCmpGlbNodes, &isShrNodes,
  // ~                       &bISGlbAndShrNodesCreated);
  pmesh.PartitionFree();
}

void GRID::LoadFromALBERTAFile(const char* filename,
                               void (*ptr)(NODE* pNode, int A),
                               bool clock_wise) {
  FILE* fp = fopen(filename, "r");
  if (!fp) {
    throw TALYException() << "ALBERTA file " << filename << " does not exist!";
  }

  char strDIM[16], strDIM_OF_WORLD[16], strNOfVertices[16], strNOfElements[16];
  PetscInt DIM, DIM_OF_WORLD;

  DIM = atoi(getParameter(fp, "DIM", strDIM, ":"));
  DIM_OF_WORLD = atoi(getParameter(fp, "DIM_OF_WORLD", strDIM_OF_WORLD, ":"));
  n_nodes_ = atoi(getParameter(fp, "number of vertices", strNOfVertices, ":"));
  n_elements_ = atoi(getParameter(fp, "number of elements",
                                  strNOfElements, ":"));

  PrintStatus("DIM: ", DIM);
  PrintStatus("DIM_OF_WORLD: ", DIM_OF_WORLD);
  PrintStatus("number of vertices: ", n_nodes());
  PrintStatus("number of elements: ", n_elements());

  set_grid_type(static_cast<GridType>(DIM));
  set_nsd(DIM_OF_WORLD);
  redimArrays(n_nodes(), n_elements());

  while (fp) {
    char value[1024];
    char* position[100];
    getLine(fp, value);

    if (strncmp(value, "vertex coordinates:", 20) == 0) {
      PrintStatus("Found vertex coordinates");
      for (PetscInt i = 0; i < n_nodes(); ++i) {
        getLine(fp, value);
        PetscReal x, y, z;
        divideParameters(value, position, "\t, ");
        x = atof(position[0]);
        y = atof(position[1]);
        z = atof(position[2]);
        // ~ PrintStatus(x, y, z);

        NODE* pNewNode = new NODE();
        pNewNode->setCoor(x, y, z);

        if (ptr) {
          ptr(pNewNode, i + 1);  // TODO: should there be a +1 here?
        }

        node_array_[i] = pNewNode;
      }
      break;
    }
  }

  SetCaredSurfaceIndicator();

  PetscInt node_no = grid_type() + 1;
  LocalNodeID *nodeIDArray = new LocalNodeID[node_no];
  while (fp) {
    char value[1024];
    char* position[100];
    getLine(fp, value);

    if (strncmp(value, "element vertices:", 20) == 0) {
      PrintStatus("Found element vertices");
      for (PetscInt i = 0; i < n_elements(); ++i) {
        getLine(fp, value);
        divideParameters(value, position, "\t, ");

        for (PetscInt j = 0; j < node_no; ++j) {
          nodeIDArray[j] = atoi(position[j]);
        }

        ElemType type;
        switch (grid_type()) {
          case kGrid2dTriangle:
            type = kElem2dTriangle;
            break;
          case kGrid3dTet:
            type = kElem3dTetrahedral;
            break;
          default:
            throw NotImplementedException() <<
                "Unsupported element type for ALBERTA!";
        }

        ELEM* pNewElm = make_elem_of_type(type);
        if (grid_type() == kGrid2dTriangle) {
          pNewElm->redim(3, nodeIDArray, true);
        } else if (grid_type() == kGrid3dTet) {
          pNewElm->redim(4, nodeIDArray, true);
        }
        pNewElm->GenSurfaceIndicator(this, cared_surface_indicator_);
        pNewElm->set_elm_id(i);

        // ~ ReOrient(this, pNewElm);
        elm_array_[i] = pNewElm;
      }
      break;
    }
  }
  delete [] nodeIDArray;

  while (fp) {
    char value[1024];
    char* position[100];
    getLine(fp, value);

    if (strncmp(value, "element neighbours:", 20) == 0) {
      PrintStatus("Found element neighbours");
      for (PetscInt i = 0; i < n_elements(); ++i) {
        getLine(fp, value);
        divideParameters(value, position, "\t, ");
        // TODO: is this the same as the previously defined node_no?
        PetscInt node_no_lcl = grid_type() + 1;

        for (PetscInt j = 0; j < node_no_lcl; ++j) {
          elm_array_[i]->set_elem_id_array(j, atoi(position[j]) + 1);
        }
      }
      break;
    }
  }

  PrintStatus("Generate element faces");
  PetscInt glbFaceID = 0;
  for (PetscInt i = 0; i < n_elements(); ++i) {
    for (PetscInt face = 0; face < elm_array_[i]->GetSurfaceCount(); ++face) {
      PetscInt ngbrID = elm_array_[i]->elem_id_array(face);

      if (ngbrID == 0)
        continue;  // element has no neighbor on this face

      if (elm_array_[i]->face_id_array(face) != 0)
        continue;  // already filled (from neighbor)
      int n_surface = elm_array_[i]->GetSurfaceCount();
      for (PetscInt faceNgbr = 0; faceNgbr < n_surface; ++faceNgbr) {
        if (i + 1 == elm_array_[ngbrID - 1]->elem_id_array(faceNgbr)) {
          PetscInt new_val = ++glbFaceID;
          elm_array_[i]->set_face_id_array(face, new_val);
          elm_array_[ngbrID - 1]->set_face_id_array(faceNgbr, new_val);
          break;
        }
      }
    }
  }
  n_faces_ = glbFaceID;
}

void GRID::LoadFromGridFile(const char* filename) {
  FILE* fp = fopen(filename, "r");
  static char value[1024];
  char* position[100];

  set_nsd(atoi(getParameter(fp, "Number of space dim.", value)));
  n_elements_ = atoi(getParameter(fp, "Number of elements", value));
  n_nodes_ = atoi(getParameter(fp, "Number of nodes", value));
  redimArrays(n_nodes(), n_elements());

  for (;;) {
    getLine(fp, value);
    if (strchr(value, '#')) {
      break;
    }
  }
  for (LocalNodeID A = 0; A < n_nodes(); A++) {
    getLine(fp, value);
    if (nsd() == 2) {
      // ~ int nodeID;
      double x, y;
      int indicatorno;

      divideParameters(value, position, "(,)[]\t ");
      // ~ nodeID = atoi(position[0]);
      x = atof(position[1]);
      y = atof(position[2]);
      indicatorno = atoi(position[3]);

      NODE* pNewNode = new NODE();
      pNewNode->setCoor(x, y, 0);
      for (int i = 0; i < indicatorno; i++) {
        pNewNode->addIndicatorNum(atoi(position[4 + i]));
      }
      node_array_[A] = pNewNode;
    } else {
      // ~ int nodeID;
      double x, y, z;
      int indicatorno;

      divideParameters(value, position, "(,)[]\t ");
      // ~ nodeID = atoi(position[0]);
      x = atof(position[1]);
      y = atof(position[2]);
      z = atof(position[3]);
      indicatorno = atoi(position[4]);

      NODE* pNewNode = new NODE();
      pNewNode->setCoor(x, y, z);
      for (int i = 0; i < indicatorno; i++) {
        pNewNode->addIndicatorNum(atoi(position[5 + i]));
      }
      node_array_[A] = pNewNode;
    }
  }

  SetCaredSurfaceIndicator();

  for (;;) {
    getLine(fp, value);
    if (strchr(value, '#')) {
      break;
    }
  }
  for (int e = 0; e < n_elements(); e++) {
    getLine(fp, value);
    // ~ int elm_id;

    divideParameters(value, position, "(,)[]\t ");

    // ~ elm_id =  atoi( position[0] );
    if (strcmp(position[1], "ElmB8n3D") == 0) {
      LocalNodeID nodeIDArray[8];
      nodeIDArray[0] = atoi(position[3]) - 1;
      nodeIDArray[1] = atoi(position[4]) - 1;
      nodeIDArray[3] = atoi(position[5]) - 1;
      nodeIDArray[2] = atoi(position[6]) - 1;
      nodeIDArray[4] = atoi(position[7]) - 1;
      nodeIDArray[5] = atoi(position[8]) - 1;
      nodeIDArray[7] = atoi(position[9]) - 1;
      nodeIDArray[6] = atoi(position[10]) - 1;

      ELEM* pNewElm = new ELEM3dHexahedral();
      pNewElm->redim(8, nodeIDArray);
      pNewElm->GenSurfaceIndicator(this, cared_surface_indicator_);
      pNewElm->set_elm_id(e);
      elm_array_[e] = pNewElm;
      set_grid_type(kGrid3dBox);
    } else if (strcmp(position[1], "ElmB4n2D") == 0) {
      LocalNodeID nodeIDArray[4];
      nodeIDArray[0] = atoi(position[3]) - 1;
      nodeIDArray[1] = atoi(position[4]) - 1;
      nodeIDArray[3] = atoi(position[6]) - 1;
      nodeIDArray[2] = atoi(position[5]) - 1;

      ELEM* pNewElm = new ELEM2dBox();
      pNewElm->redim(4, nodeIDArray);
      pNewElm->GenSurfaceIndicator(this, cared_surface_indicator_);
      pNewElm->set_elm_id(e);
      ReOrient(this, pNewElm);
      elm_array_[e] = pNewElm;
      set_grid_type(kGrid2dBox);
    } else if (strcmp(position[1], "ElmT3n2D") == 0) {
      LocalNodeID nodeIDArray[3];
      nodeIDArray[0] = atoi(position[3]) - 1;
      nodeIDArray[1] = atoi(position[4]) - 1;
      nodeIDArray[2] = atoi(position[5]) - 1;

      ELEM* pNewElm = new ELEM2dTriangle();
      pNewElm->redim(3, nodeIDArray);
      pNewElm->GenSurfaceIndicator(this, cared_surface_indicator_);
      pNewElm->set_elm_id(e);
      ReOrient(this, pNewElm);
      elm_array_[e] = pNewElm;
      set_grid_type(kGrid2dTriangle);
    } else if (strcmp(position[1], "ElmT4n3D") == 0) {
      LocalNodeID nodeIDArray[4];
      nodeIDArray[0] = atoi(position[3]) - 1;
      nodeIDArray[1] = atoi(position[4]) - 1;
      nodeIDArray[2] = atoi(position[5]) - 1;
      nodeIDArray[3] = atoi(position[6]) - 1;

      ELEM* pNewElm = new ELEM3dTetrahedral();
      pNewElm->redim(4, nodeIDArray);
      pNewElm->GenSurfaceIndicator(this, cared_surface_indicator_);
      pNewElm->set_elm_id(e);
      elm_array_[e] = pNewElm;
      set_grid_type(kGrid3dTet);
    } else {
      throw NotImplementedException() << "Unsupported element type " <<
                                         position[1];
    }
  }
  fclose(fp);
}

// File naming scheme: [filename].[mpi_rank]

// Format:
// [num_total_nodes] [num_local_nodes] [num_owned_nodes] [num_local_elements]\n

// [x] [y] [z] [boundary_indicators] \n        (z is included but 0 in 2D case)
// ... (repeat for num_local_nodes entries) ...

// [local_elem_id] [node0] [node1] [...] \n
// ... (repeat for num_local_elements entries) ...

// [global ID for local node 0] [global ID for local node 1] ...
//      (repeat for num_local_nodes entries) ... \n
// [solution ID for local node 0] [solution ID for local node 1] ...
//      (repeat for num_local_nodes entries) ... \n

// [num_shared_nodes on process 0] [num_shared_nodes on process 1] ...
//      (repeat for mpi_size entries) ... \n

// [local ID of shared node] ... (repeat for num_shared_nodes entries) ... \n

// [num_processes_shared_with]
//    [first process shared with] [local node id on that process]
//    [second process shared with] [local node id on that process]
//    ... (repeat for num_processes_shared_with entries) ... \n
// ... (repeat for num_shared_nodes entries) ...

int GRID::LoadFromParallel(const char* prefix, GridType new_grid_type,
                           const int* elems_per_direction, int order) {
  char filename[1024];
  int ferr;
  int mpi_rank = GetMPIRank();
  int mpi_size = GetMPISize();
  snprintf(filename, sizeof(filename), "%s.%d", prefix, mpi_rank);
  FILE* fp = fopen(filename, "r");

  if (!fp) {
    throw FileIOException() << "Parallel file " << filename << " not found!";
  }

  PhysicalNodeID total_node_count;
  int owned_node_count;
  // TODO: is %d always correct for total_node_count?
  ferr = fscanf(fp, "%" PETSCINT_F "%d%d%d", &total_node_count, &n_nodes_,
                &owned_node_count, &n_elements_);
  redimArrays(n_nodes(), n_elements());

  /// set parallel info
  // (some of them will be cleaned in redim, so reset here)
  parallel_type_ = kWithDomainDecomp;
  set_n_subgrids(mpi_size);
  set_grid_id(mpi_rank);
  set_n_total_nodes(total_node_count);
  set_n_owned_nodes(owned_node_count);
  set_grid_type(new_grid_type);

  /// read  all nodes
  for (LocalNodeID A = 0; A < n_nodes(); A++) {
    double x = 0.0, y = 0.0, z = 0.0;
    ferr = fscanf(fp, "%lf%lf%lf", &x, &y, &z);
    // make sure we read all the expected values, 3 coordinates
    assert(ferr == 3);

    NodeIndicator boundaryIndicators;
    ferr = fscanf(fp, NODE_INDICATOR_FORMAT, &boundaryIndicators);

    NODE* pNewNode = new NODE();
    pNewNode->setCoor(x, y, z);
    pNewNode->setIndicators(boundaryIndicators);
    node_array_[A] = pNewNode;
  }
  SetCaredSurfaceIndicator();

  // read  all elements
  ElemType elmType = grid_to_elem_type(grid_type());
  for (int e = 0; e < n_elements(); e++) {
    int type = -1;
    LocalNodeID nodeIDArray[27];
    ELEM* pNewElm = make_elem_of_type(elmType);
    if (grid_type() == kGrid3dBox) {
      switch (order) {
        case 1: {
#ifdef PETSC_USE_64BIT_INDICES
          ferr = fscanf(fp, "%d%d%d%d%d%d%d%d%d", &type,
              nodeIDArray+0, nodeIDArray+1, nodeIDArray+2, nodeIDArray+3,
              nodeIDArray+4, nodeIDArray+5, nodeIDArray+6, nodeIDArray+7);
#else
          ferr = fscanf(fp, "%d%d%d%d%d%d%d%d%d", &type, nodeIDArray + 0,
                        nodeIDArray + 1, nodeIDArray + 2, nodeIDArray + 3,
                        nodeIDArray + 4, nodeIDArray + 5, nodeIDArray + 6,
                        nodeIDArray + 7);
#endif
          // make sure we read all the expected values, 8 nodes + 1 type
          assert(ferr == 9);

          pNewElm->redim(8, nodeIDArray);
          break;
        }
        case 2: {
#ifdef PETSC_USE_64BIT_INDICES
          ferr = fscanf(fp,
              "%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
              &type,
              nodeIDArray+0, nodeIDArray+1, nodeIDArray+2, nodeIDArray+3,
              nodeIDArray+4, nodeIDArray+5, nodeIDArray+6, nodeIDArray+7,
              nodeIDArray+8, nodeIDArray+9, nodeIDArray+10, nodeIDArray+11,
              nodeIDArray+12, nodeIDArray+13, nodeIDArray+14, nodeIDArray+15,
              nodeIDArray+16, nodeIDArray+17, nodeIDArray+18, nodeIDArray+19,
              nodeIDArray+20, nodeIDArray+21, nodeIDArray+22, nodeIDArray+23,
              nodeIDArray+24, nodeIDArray+25, nodeIDArray+26);
#else
          ferr = fscanf(fp,
              "%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
              &type,
              nodeIDArray+0, nodeIDArray+1, nodeIDArray+2, nodeIDArray+3,
              nodeIDArray+4, nodeIDArray+5, nodeIDArray+6, nodeIDArray+7,
              nodeIDArray+8, nodeIDArray+9, nodeIDArray+10, nodeIDArray+11,
              nodeIDArray+12, nodeIDArray+13, nodeIDArray+14, nodeIDArray+15,
              nodeIDArray+16, nodeIDArray+17, nodeIDArray+18, nodeIDArray+19,
              nodeIDArray+20, nodeIDArray+21, nodeIDArray+22, nodeIDArray+23,
              nodeIDArray+24, nodeIDArray+25, nodeIDArray+26);
#endif
          // make sure we read all the expected values, 27 nodes + 1 type
          assert(ferr == 28);
          pNewElm->redim(27, nodeIDArray);
          break;
        }
        default: {
          throw NotImplementedException() <<
              "GRID::loadFromParallel does not support the given basis order "
              "for 3D grid box";
        }
      }
      set_nsd(3);
    } else if (grid_type() == kGrid2dBox) {
      if (order != 1) {
        throw NotImplementedException() <<
            "GRID::loadFromParallel only supports basis 1 for "
            "2D grid box";
      }
#ifdef PETSC_USE_64BIT_INDICES
      ferr = fscanf(fp, "%d%d%d%d%d", &type,
          nodeIDArray+0, nodeIDArray+1, nodeIDArray+2, nodeIDArray+3);
#else
      ferr = fscanf(fp, "%d%d%d%d%d", &type, nodeIDArray + 0, nodeIDArray + 1,
                    nodeIDArray + 2, nodeIDArray + 3);
#endif
      // make sure we read all the expected values, 4 nodes + 1 type
      assert(ferr == 5);

      pNewElm->redim(4, nodeIDArray);
      set_nsd(2);
    } else if (grid_type() == kGrid2dTriangle) {
      if (order != 1) {
        throw NotImplementedException() <<
            "GRID::loadFromParallel only supports basis 1 for "
            "2D triangle grid box";
      }
#ifdef PETSC_USE_64BIT_INDICES
      ferr = fscanf(fp, "%d%d%d%d", &type,
          nodeIDArray+0, nodeIDArray+1, nodeIDArray+2);
#else
      ferr = fscanf(fp, "%d%d%d%d", &type, nodeIDArray + 0, nodeIDArray + 1,
                    nodeIDArray + 2);
#endif
      // make sure we read all the expected values, 3 nodes + 1 type
      assert(ferr == 4);

      pNewElm->redim(3, nodeIDArray);
      set_nsd(2);
    } else if (grid_type() == kGrid3dTet) {
      if (order != 1) {
        throw NotImplementedException() <<
            "GRID::loadFromParallel only supports basis 1 for "
            "3D tet grid";
      }
#ifdef PETSC_USE_64BIT_INDICES
      ferr = fscanf(fp, "%d%d%d%d%d", &type,
          nodeIDArray+0, nodeIDArray+1, nodeIDArray+2, nodeIDArray+3);
#else
      ferr = fscanf(fp, "%d%d%d%d%d", &type, nodeIDArray + 0, nodeIDArray + 1,
                    nodeIDArray + 2, nodeIDArray + 3);
#endif
      // make sure we read all the expected values, 4 nodes + 1 type
      assert(ferr == 5);

      pNewElm->redim(4, nodeIDArray);
      set_nsd(3);
    } else {
      throw NotImplementedException() <<
          "Unhandled gridType in GRID::loadFromParallel!";
    }

    pNewElm->GenSurfaceIndicator(this, cared_surface_indicator_);
    pNewElm->set_elm_id(type);
    elm_array_[e] = pNewElm;
  }

  // read the global node IDs
  physical_map_.redim(n_nodes());
  for (LocalNodeID A = 0; A < n_nodes(); A++) {
    PhysicalNodeID gid;
    ferr = fscanf(fp, "%" PETSCINT_F, &gid);
    // make sure we read all the expected values, 1 id
    assert(ferr == 1);

    physical_map_(A) = gid;
  }

  // read solution node IDs
  solution_map_.redim(n_nodes());
  for (LocalNodeID A = 0; A < n_nodes(); A++) {
    SolutionNodeID sid;
    ferr = fscanf(fp, "%" PETSCINT_F, &sid);
    // make sure we read all the expected values, 1 id
    assert(ferr == 1);

    solution_map_(A) = sid;
  }

  ////////////////// read cmu node nos
  n_shared_nodes_.redim(n_subgrids());
  for (int g = 0; g < n_subgrids(); g++) {
    int cmunodeno;
    ferr = fscanf(fp, "%d", &cmunodeno);
    // make sure we read all the expected values, 1 count
    assert(ferr == 1);

    n_shared_nodes_(g) = cmunodeno;
  }
  shared_nodes_.redim(n_shared_nodes_(grid_id()));
  for (int i = 0; i < shared_nodes_.size(); i++) {
    int id;
    ferr = fscanf(fp, "%d", &id);
    // make sure we read all the expected values, 1 id
    assert(ferr == 1);

    shared_nodes_(i) = id;
  }

  ////////////////// read node belong
  node_belong_.redim(n_nodes());
  for (LocalNodeID A = 0; A < n_nodes(); A++) {
    node_belong_(A).is_owned = true;  // assume this process owns this
  }
  for (int i = 0; i < shared_nodes_.size(); i++) {
    int A = shared_nodes_(i);
    // 2 1 6 0 4
    int no;
    int procs[1024];
    int cmids[1024];
    ferr = fscanf(fp, "%d", &no);
    // make sure we read all the expected values, 1 value
    assert(ferr == 1);

    procs[0] = mpi_rank;
    cmids[0] = i;
    for (int j = 1; j <= no; j++) {
      ferr = fscanf(fp, "%d%d", procs + j, cmids + j);
      // make sure we read all the expected values, 1 value
      assert(ferr == 2);
    }
    no++;
    for (int j = 0; j < no; j++) {
      for (int k = j + 1; k < no; k++) {
        if (procs[k] < procs[j]) {
          std::swap(procs[j], procs[k]);
          std::swap(cmids[j], cmids[k]);
        }
      }
    }
    for (int j = 0; j < no; j++) {
      ShareDetails sharedetails(procs[j], cmids[j]);
      node_belong_(A).share_data.appendData(sharedetails);
    }
    if (node_belong_(A).share_data(0).grid_id() < grid_id()) {
      node_belong_(A).is_owned = false;  // process does not own this
    }
  }
  fclose(fp);

  // set up the number of elements in each direction, if appropriate
  if (elems_per_direction != NULL) {
    // include false because we don't want to overwrite the number of nodes
    // and elements in the system.
    SetNodeElemCounts(elems_per_direction, false);
  }
  return ferr;
}

void GRID::LoadTriangleFormat(const char* filename) {
  char nodefile[1024];
  char polyfile[1024];
  char elefile[1024];

  char buf[1024];
  char* position[100];
  int edgeno;
  snprintf(nodefile, sizeof(nodefile), "%s.node", filename);
  snprintf(polyfile, sizeof(polyfile), "%s.poly", filename);
  snprintf(elefile, sizeof(elefile), "%s.ele", filename);
  FILE* fpnode = fopen(nodefile, "r");
  FILE* fppoly = fopen(polyfile, "r");
  FILE* fpele = fopen(elefile, "r");

  getLine(fpnode, buf);
  divideParameters(buf, position, "(,)[]\t ");
  n_nodes_ = atoi(position[0]);
  set_nsd(atoi(position[1]));

  getLine(fpele, buf);
  divideParameters(buf, position, "(,)[]\t ");
  n_elements_ = atoi(position[0]);

  getLine(fppoly, buf);
  getLine(fppoly, buf);
  divideParameters(buf, position, "(,)[]\t ");
  edgeno = atoi(position[0]);

  redimArrays(n_nodes(), n_elements());
  for (LocalNodeID A = 0; A < n_nodes(); A++) {
    getLine(fpnode, buf);
    divideParameters(buf, position, "(,)[]\t ");

//      int nodeID = atoi(position[0]);
    double x = atof(position[1]);
    double y = atof(position[2]);
    double z = 0;
    NODE* pNewNode = new NODE();
    pNewNode->setCoor(x, y, z);
    node_array_[A] = pNewNode;
  }

  for (int e = 0; e < edgeno; e++) {
    getLine(fppoly, buf);
    int no = divideParameters(buf, position, "(,)[]\t ");
//      int edgeID = atoi(position[0]);
    int A = atoi(position[1]);
    int B = atoi(position[2]);
    int indicator = 0;
    if (no >= 4) {
      indicator = atoi(position[3]);
    }
    node_array_[A - 1]->addIndicatorNum(indicator);
    node_array_[B - 1]->addIndicatorNum(indicator);
  }
  SetCaredSurfaceIndicator();

  for (int e = 0; e < n_elements(); e++) {
    getLine(fpele, buf);
    divideParameters(buf, position, "(,)[]\t ");
    LocalNodeID nodeIDArray[3];
    nodeIDArray[0] = atoi(position[1]) - 1;
    nodeIDArray[1] = atoi(position[2]) - 1;
    nodeIDArray[2] = atoi(position[3]) - 1;

    ELEM* pNewElm = new ELEM2dTriangle();
    pNewElm->redim(3, nodeIDArray);
    pNewElm->GenSurfaceIndicator(this, cared_surface_indicator_);
    pNewElm->set_elm_id(e);
    elm_array_[e] = pNewElm;
    ReOrient(this, pNewElm);
  }

  set_grid_type(kGrid2dTriangle);

  fclose(fpnode);
  fclose(fppoly);
  fclose(fpele);
  return;
}

void GRID::Move(const ZEROARRAY<double>& solution, double dt) {
  for (LocalNodeID i_node = 0; i_node < n_nodes(); i_node++) {
    for (int dir = 0; dir < nsd(); dir++) {
      node_array_[i_node]->setCoor(dir,
          node_array_[i_node]->getCoor(dir) +
          solution(i_node * nsd() + dir) * dt);
    }
  }
}

void GRID::Move(double dx, double dy, double dz) {
  const double du[3] = { dx, dy, dz };
  for (LocalNodeID i_node = 0; i_node < n_nodes(); i_node++) {
    for (int dir = 0; dir < nsd(); dir++) {
      node_array_[i_node]->setCoor(dir,
                                  node_array_[i_node]->getCoor(dir) + du[dir]);
    }
  }
}

void GRID::Scale(double dx, double dy, double dz) {
  const double scale[3] = { dx, dy, dz };
  for (LocalNodeID i_node = 0; i_node < n_nodes(); i_node++) {
    for (int dir = 0; dir < nsd(); dir++) {
      double coord_old = node_array_[i_node]->getCoor(dir);
      node_array_[i_node]->setCoor(dir, coord_old * scale[dir]);
    }
  }
}

void GRID::MoveTo(const ZEROARRAY<double>& solution, GRID& new_grid,
                  double dt) {
  // copy the current grid into new_grid
  new_grid = *this;

  // move new_grid
  new_grid.Move(solution, dt);
}

void GRID::PrintElmSurfaceIndicator() const {
  for (int i = 0; i < n_elements(); i++) {
    elm_array_[i]->PrintSurfaceIndicator(i + 1);
  }
}

void GRID::PrintGridFile(const char* filename) const {
  if ((grid_type() != kGrid2dTriangle) && (grid_type() != kGrid2dBox)) {
    throw NotImplementedException() <<
        "Can only support print of kGrid2dTriangle and kGrid2dBox";
  }

  FILE* fp = fopen(filename, "w");

  fprintf(fp, "\n\nFinite element mesh (GRID):\n\n");
  fprintf(fp, "Number of space dim. = %d\n", nsd());
  fprintf(fp, "Number of elements = %d\n", n_elements());
  fprintf(fp, "Number of nodes = %d\n\n", n_nodes());
  fprintf(fp, "All elements are of the same type : dpTRUE\n");
  fprintf(fp, "Max number of nodes in an element: 3\n");
  fprintf(fp, "Only one subdomain               : dpTRUE\n");
  fprintf(fp, "Lattice data                     ? 1\n");
  fprintf(fp, "d=2 domain=[0,3.6]x[0,2] indices=[1:21]x[1:11]\n");
  fprintf(fp, "%d Boundary indicators: \n", cared_surface_indicator_.size());
  for (int i = 0; i < cared_surface_indicator_.size(); i++) {
    fprintf(fp, "Ind=%d\n", cared_surface_indicator_(i));
  }

  fprintf(fp, "\nNodal coordinates and nodal boundary indicators,\n");
  fprintf(fp, "the columns contain:\n");
  fprintf(fp, "- node number\n");
  fprintf(fp, "- coordinates\n");
  fprintf(fp, "- no of boundary indicators that are set (ON)\n");
  fprintf(fp, "- the boundary indicators that are set (ON) if any.\n");
  fprintf(fp, "# \n");

  for (LocalNodeID A = 0; A < n_nodes(); A++) {
    if (nsd() == 2) {
      NODE* pNode = node_array_[A];
      fprintf(fp, "%d ( %e, %e) [%d] ", A, pNode->x(), pNode->y(),
              pNode->getIndicatorNo());
      for (unsigned int i = 0; i < MAX_NODE_INDICATORS; i++) {
        if (pNode->indicators() & INDICATOR_NUM(i))
          fprintf(fp, "%d ", i);
      }
      fprintf(fp, "\n");
    } else {
      throw NotImplementedException() << "PrintGridFile does not support 3D";
    }
  }

  fprintf(fp, "\nElement types and connectivity\n");
  fprintf(fp, "the columns contain:\n");
  fprintf(fp, "- element number\n");
  fprintf(fp, "- element type\n");
  fprintf(fp, "- subdomain number\n");
  fprintf(fp, "- the global node numbers of the nodes in the element.\n");
  fprintf(fp, "# \n");

  for (int e = 0; e < n_elements(); e++) {
    if (nsd() == 2) {
      ELEM* pElm = elm_array_[e];
      switch (pElm->elmType()) {
        case kElem2dTriangle:
#ifdef PETSC_USE_64BIT_INDICES
          fprintf(fp, "%d ElmT3n2D 1 %d %d %d\n", e,
                  pElm->node_id_array(0)+1, pElm->node_id_array(1)+1,
                  pElm->node_id_array(2)+1);
#else
          fprintf(fp, "%d ElmT3n2D 1 %d %d %d\n", e, pElm->node_id_array(0) + 1,
                  pElm->node_id_array(1) + 1, pElm->node_id_array(2) + 1);
#endif
          break;
        case kElem2dBox:
#ifdef PETSC_USE_64BIT_INDICES
          fprintf(fp, "%d ElmB4n2D 1 %d %d %d %d\n", e,
                  pElm->node_id_array(0)+1, pElm->node_id_array(1)+1,
                  pElm->node_id_array(2)+1, pElm->node_id_array(3)+1);
#else
          fprintf(fp, "%d ElmB4n2D 1 %d %d %d %d\n", e,
                  pElm->node_id_array(0) + 1, pElm->node_id_array(1) + 1,
                  pElm->node_id_array(2) + 1, pElm->node_id_array(3) + 1);
#endif
          break;
        default:
          throw NotImplementedException() <<
              "Can only support print of kGrid2dTriangle and kGrid2dBox";
      }
    } else {
      throw NotImplementedException() << "PrintGridFile does not support 3D";
    }
  }
  fclose(fp);
}

int GRID::GetGridIdOfElementOwner(int elm_id) const {
  if (parallel_type_ == kNoDomainDecomp) {
    int mpi_size = GetMPISize();  // TODO: replace this with n_subgrids?

    // this is the minimum number of elements on each process. The first process
    // will have more than this if the number of elements is not an even
    // multiple of the number of processes.
    int min_elements_per_process = n_elements() / mpi_size;
    // this is the proposed grid id, it may be corrected below if this ends
    // up being a negative value (due to the calculation method, this is
    // possible if the element is one of the last ones in the list of elements)
    int proposed_id = mpi_size - 1 - (elm_id / min_elements_per_process);
    return (proposed_id < 0) ? 0 : proposed_id;
  }
  throw NotImplementedException() <<
      "GetGridIdOfElementOwner does not support domain decomposition";
}

void GRID::ReadGrid(const InputData* pIdata) {
  if (pIdata->inputFilenameGrid.empty()) {
    PrintWarning("Grid not generated nor loaded!");
  }

  std::vector< FileIO<NODEData>* > io = make_io_options<NODEData>();
  for (unsigned int i = 0; i < io.size(); i++) {
    if (io.at(i)->is_this_type(pIdata->inputFilenameGrid.c_str())) {
      try {
        io.at(i)->load_grid(this, pIdata);
        PrintStatus("Mesh loaded successfully");
      } catch(...) {
        free_vec_contents(io);
        throw;
      }
      free_vec_contents(io);
      return;
    }
  }

  free_vec_contents(io);
  throw FileIOException() << "Could not match file type for grid input file \""
      << pIdata->inputFilenameGrid << "\".";
}

void GRID::redimArrays(int nodeno, int elmno) {
  Cleanup();

  n_elements_ = elmno;

  elm_array_ = static_cast<ELEM**>(malloc(sizeof(ELEM*) * n_elements()));
  for (int i = 0; i < n_elements(); i++) {
    elm_array_[i] = NULL;
  }

  n_nodes_ = nodeno;
  node_array_ = static_cast<NODE**>(malloc(sizeof(NODE*) * nodeno));
  for (LocalNodeID i = 0; i < n_nodes(); i++) {
    node_array_[i] = NULL;
  }

  set_n_subgrids(1);
  set_grid_id(0);
  set_n_total_nodes(n_nodes());
  set_n_owned_nodes(n_nodes());
}

void GRID::redimIndicator(const ZEROARRAY<ZEROARRAY<int> >& newIndicator,
                          const ZEROARRAY<int>& newIndValue,
                          bool reset_surface_indicators) {
  for (LocalNodeID i_node = 0; i_node < n_nodes(); i_node++) {
    NODE* pNode = node_array_[i_node];

    NodeIndicator newIndicators = 0;
    for (int i = 0; i < newIndicator.size(); i++) {
      for (unsigned int j = 0; j < MAX_NODE_INDICATORS; j++) {
        if (pNode->BoNode(j) && newIndicator(i).contains(j)) {
          newIndicators |= INDICATOR_NUM(newIndValue(i));
          break;
        }
      }
    }
    pNode->setIndicators(newIndicators);
  }

  if (!reset_surface_indicators) {
    return;
  }
  SetCaredSurfaceIndicator();
  for (int e = 0; e < n_elements(); e++) {
    ELEM* pElm = elm_array_[e];
    pElm->GenSurfaceIndicator(this, cared_surface_indicator_);
  }
}

void GRID::ReOrderNodes() {
  ResetNodeOwner();
  GetNumOwnedNodes();
  ZEROARRAY<int> newNodeIDArray;  // starts from 0
  newNodeIDArray.redim(n_nodes());
  int curNodeID = 0;
  int mpi_size = GetMPISize();
  for (int process_id = 0; process_id < mpi_size; process_id++) {
    for (int elm_id = 0; elm_id < n_elements(); elm_id++) {
      if (GetGridIdOfElementOwner(elm_id) != process_id) {
        continue;
      }
      ELEM* pElm = GetElm(elm_id);
      for (ElemNodeID localNodeID = 0; localNodeID < pElm->n_nodes();
           localNodeID++) {
        LocalNodeID nodeID = pElm->node_id_array(localNodeID);  // start from 0
        if (node_array_[nodeID]->elem_id_ == elm_id) {
          newNodeIDArray(nodeID) = curNodeID;
          curNodeID++;
        }
      }
    }
  }

  newNodeIDArray.print(5);

  for (int e = 0; e < n_elements(); e++) {
    ELEM* pElm = GetElm(e);
    for (ElemNodeID localNodeID = 0; localNodeID < pElm->n_nodes();
         localNodeID++) {
      LocalNodeID nodeID = pElm->node_id_array(localNodeID);  // start from 0
      int newNodeID = newNodeIDArray(nodeID);
      pElm->set_node_id_array(localNodeID, newNodeID);
    }
  }

  ZEROARRAY<NODE*> newNodeArray;
  newNodeArray.redim(n_nodes());
  for (LocalNodeID i = 0; i < n_nodes(); i++) {
    newNodeArray(i) = node_array_[i];
  }
  for (LocalNodeID i = 0; i < n_nodes(); i++) {
    node_array_[newNodeIDArray(i)] = newNodeArray(i);
  }
}

void GRID::SetCaredSurfaceIndicator() {
  cared_surface_indicator_.cleanup();
  for (int i_node = 0; i_node < n_nodes(); i_node++) {
    NODE* pNode = node_array_[i_node];
    // TODO: if surface indicators are also bitmask, this can be reaplces by
    // a single logic operation
    for (unsigned int j = 0; j < MAX_NODE_INDICATORS; j++) {
      if (pNode->BoNode(j)) {
        AddCaredSurfaceIndicator(j);
      }
    }
  }
}

void GRID::SetNodeCoords(const ZEROARRAY<double>& position) {
  assert(position.size() >= n_nodes() * nsd());
  for (LocalNodeID i_node = 0; i_node < n_nodes(); i_node++) {
    for (int dir = 0; dir < nsd(); dir++) {
      int idx = i_node * nsd() + dir;
      double coordinate_value = position(idx);

      node_array_[i_node]->setCoor(dir, coordinate_value);
    }
  }
}

// ******************************
//        Setter functions
// ******************************

void GRID::set_nsd(int new_nsd) {
  if (new_nsd < 1 || new_nsd > 3) {
    throw TALYException() << "Invalid spatial dimension: " << new_nsd;
  }
  nsd_ = new_nsd;
  if (n_nodes_per_direction_ != NULL) { delete [] n_nodes_per_direction_; }
  if (n_elems_per_direction_ != NULL) { delete [] n_elems_per_direction_; }
  n_nodes_per_direction_ = new int[new_nsd];
  n_elems_per_direction_ = new int[new_nsd];
}

// ******************************
//      Protected functions
// ******************************

void GRID::CreateElementsBasis1() {
  throw NotImplementedException() <<
      "Basis order 1 not supported in grid redim";
}

void GRID::CreateElementsBasis2() {
  throw NotImplementedException() <<
      "Basis order 2 not supported in grid redim";
}

void GRID::CreateElementsBasis3() {
  throw NotImplementedException() <<
      "Basis order 3 not supported in grid redim";
}

void GRID::CreateNodes(const double* dimensions) {
  throw NotImplementedException() << "CreateNodes not supported in grid redim";
}

void GRID::redim(const double* dimensions, const int* n_elems) {
  ValidateParams(dimensions, n_elems);

  SetNodeElemCounts(n_elems);

  redimArrays(n_nodes(), n_elements());

  CreateNodes(dimensions);

  switch (basis_order()) {
    case 1:
      CreateElementsBasis1();
      break;
    case 2:
      CreateElementsBasis2();
      break;
    case 3:
      CreateElementsBasis3();
      break;
    default:
      throw NotImplementedException() <<
          "Unsupported basis order for redim of grid type: " << GridTypeName();
  }
  SetCaredSurfaceIndicator();
  GenElmSurfaceIndicator();
}

// default implementation, throw not implemented
void GRID::redimDD(const double* dimensions, const int* n_elems) {
  throw NotImplementedException() << "Grid NOT generated, "
      "Domain Decomposition requested. "
      "(option not implemented)";
}

void GRID::SetNodeElemCounts(const int* n_elems, bool set_totals) {
  PetscInt n_nodes_product = 1;
  PetscInt n_elems_product = 1;
  for (int i = 0; i < nsd(); i++) {
    int n_nodes_this_direction = n_elems[i] + 1;
    int n_elems_this_direction = n_elems[i] / basis_order();
    // set # of nodes and elems for the ith direction
    set_n_nodes_per_direction(i, n_nodes_this_direction);
    set_n_elems_per_direction(i, n_elems[i]);  // requested number of elements
    // calculate the total number of nodes and elements
    n_nodes_product *= n_nodes_this_direction;
    n_elems_product *= n_elems_this_direction;  // true number of elements
  }
  if (set_totals) {
    set_n_nodes(n_nodes_product);
    set_n_elements(n_elems_product);
  }
}

void GRID::ValidateParams(const double* dimensions, const int* n_elems) const {
  for (int i = 0; i < nsd(); i++) {
    assert(dimensions[i] > 0.0);
    assert(n_elems[i] > 0);
    if (n_elems[i] % basis_order() != 0) {
      throw TALYException() << "Incorrect element count. "
                               "Must be a multiple of the basis function.";
    }
  }
}

// ******************************
//      Private functions
// ******************************

void GRID::AddCaredSurfaceIndicator(int indicator) {
  if (!cared_surface_indicator_.contains(indicator)) {
    cared_surface_indicator_.appendData(indicator);
  }
}

void GRID::CalcBoundarySegmentLength() {
  double previousTotalLen = 0;

  for (BoundaryList_type::iterator it = boundary_.begin();
       it != boundary_.end(); it++) {
    it->calcLen(this);
    it->previousTotalLen = previousTotalLen;
    previousTotalLen += it->len;
  }
}

void GRID::Cleanup() {

  if (elm_array_) {
    for (int i = 0; i < n_elements_; i++) {
      if (elm_array_[i]) {
        delete elm_array_[i];
      }
    }
    free(elm_array_);
    elm_array_ = NULL;
    n_elements_ = 0;
  }

  if (node_array_) {
    for (LocalNodeID i = 0; i < n_nodes(); i++) {
      if (node_array_[i]) {
        delete node_array_[i];
      }
    }
    free(node_array_);
    node_array_ = NULL;
    n_nodes_ = 0;
  }
}

void GRID::GenBoundary(int nStartNode, int nEndNode) {
  nStartNode--;
  if ((grid_type() != kGrid2dTriangle) && (grid_type() != kGrid2dBox)) {
    throw NotImplementedException() << "GRID::GenBoundary only supports "
        "2dTriangle and 2dBox";
  }
  if (nEndNode == -1) {
    nEndNode = nStartNode;
  } else {
    nEndNode--;
  }

  boundary_.clear();
  int nCurNode = nStartNode;
  int nNextNode = -1;

  std::vector<ELEM*> elmList;
  std::vector<int> preNodeList;
  std::vector<int> nextNodeList;

  while (true) {
    elmList.clear();
    preNodeList.clear();
    nextNodeList.clear();

    for (int e = 0; e < n_elements(); e++) {
      int preNode, nextNode;
      if (GetElm(e)->Contains(nCurNode, preNode, nextNode)) {
        ELEM* pElm = GetElm(e);
        elmList.push_back(pElm);
        preNodeList.push_back(preNode);
        nextNodeList.push_back(nextNode);
      }
    }

    ELEM* pElm = NULL;

    std::vector<int>::iterator next = nextNodeList.begin();
    std::vector<ELEM*>::iterator elm = elmList.begin();
    while (next != nextNodeList.end() && elm != elmList.end()) {
      pElm = *elm;
      nNextNode = *next;

      bool ok = true;
      for (std::vector<int>::iterator pre = preNodeList.begin();
          pre != preNodeList.end(); pre++) {
        if (nNextNode == *pre) {
          ok = false;
          break;
        }
      }

      if (ok) {
        break;
      }

      next++;
      elm++;
    }

    assert(pElm != NULL);

    Segment new_segment;
    new_segment.nElm = pElm->elm_id();
    new_segment.nStartNode = nCurNode;
    new_segment.nEndNode = nNextNode;
    boundary_.push_back(new_segment);
    nCurNode = nNextNode;

    if (nNextNode == nEndNode)
      break;
  }

  CalcBoundarySegmentLength();

  for (BoundaryList_type::iterator it = boundary_.begin();
      it != boundary_.end(); it++) {
    printf("%d-(%d)->%d   %e %e\n", it->nStartNode, it->nElm,
           it->nEndNode, it->previousTotalLen, it->len);
  }
  printf("\n");
}

ElemNodeID GRID::GetNumNodesInElm(int e) const {
  return GetElm(e)->n_nodes();
}

void GRID::ResetNodeOwner() {
  for (int e = 0; e < n_elements(); e++) {
    ELEM* pElm = GetElm(e);
    for (ElemNodeID localNodeID = 0; localNodeID < pElm->n_nodes();
         localNodeID++) {
      LocalNodeID nodeID = pElm->node_id_array(localNodeID);
      node_array_[nodeID]->elem_id_ = e;
    }
  }
}

}  // namespace TALYFEMLIB
