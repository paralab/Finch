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
#include <talyfem/grid/elem.h>

#include <algorithm>  // for std::max

#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/node.h>

namespace TALYFEMLIB {

ELEM::ELEM()
    : surface_indicator_(),
      node_id_array_(NULL),
      elem_id_array_(NULL),
      face_id_array_(NULL),
      elm_id_(-1),
      n_nodes_(0) {
}

ELEM::~ELEM() {
  cleanup();
}

void ELEM::cleanup() {
  if (elem_id_array_) {
    free(elem_id_array_);
    elem_id_array_ = NULL;
  }
  if (face_id_array_) {
    free(face_id_array_);
    face_id_array_ = NULL;
  }
  if (node_id_array_) {
    free(node_id_array_);
    node_id_array_ = NULL;
    n_nodes_ = 0;
  }
}

const ZEROPTV& ELEM::get_node_loc(const GRID * grid, int node_idx) const {
  return grid->GetNode(node_id_array(node_idx))->location();
}

void ELEM::redim(ElemNodeID nodeno, const LocalNodeID* pNodeIDArray,
                 bool withNeighbor,
                 const PetscInt* pElemIDArray, const PetscInt* pFaceIDArray) {
  cleanup();

  set_n_nodes(nodeno);
  node_id_array_ = static_cast<LocalNodeID*>
  (malloc(sizeof(LocalNodeID) * n_nodes()));

  if (pNodeIDArray) {
    for (ElemNodeID i = 0; i < n_nodes(); i++) {
      node_id_array_[i] = pNodeIDArray[i];
    }
  }

  if (withNeighbor) {
    elem_id_array_ = static_cast<PetscInt*>
    (malloc(sizeof(PetscInt) * GetSurfaceCount()));
    face_id_array_ = static_cast<PetscInt*>
    (malloc(sizeof(PetscInt) * GetSurfaceCount()));
  } else {
    elem_id_array_ = NULL;
    face_id_array_ = NULL;
  }

  if (withNeighbor && pElemIDArray != NULL) {
    for (int i = 0; i < GetSurfaceCount(); i++) {
      elem_id_array_[i] = pElemIDArray[i];
    }
  }

  if (withNeighbor && pFaceIDArray != NULL) {
    for (int i = 0; i < GetSurfaceCount(); i++) {
      face_id_array_[i] = pFaceIDArray[i];
    }
  } else if (withNeighbor && pFaceIDArray == NULL) {
    for (int i = 0; i < GetSurfaceCount(); i++) {
      face_id_array_[i] = 0;
    }
  }
}

LocalNodeID ELEM::ElemToLocalNodeID(ElemNodeID idx) const {
  assert(idx >= 0 && idx < n_nodes());
  return node_id_array_[idx];
}

bool ELEM::Contains(LocalNodeID nodeID, LocalNodeID& preNodeID,
                    LocalNodeID& nextNodeID) const {
  // find which of our nodes has ID nodeID,
  // store the index of that node in A
  ElemNodeID A = 0;
  while (A < n_nodes()) {
    if (node_id_array_[A] == nodeID) {
      break;
    }
    A++;
  }

  if (A == n_nodes()) {  // we didn't find it
    return false;
  }
  // TODO: is this correct??? it is returning type ElemNodeID
  preNodeID = (A + n_nodes() - 1) % n_nodes();
  nextNodeID = (A + 1) % n_nodes();
  return true;
}

bool ELEM::Contains(LocalNodeID nodeID) const {
  for (ElemNodeID a = 0; a < n_nodes(); a++) {
    if (node_id_array_[a] == nodeID) {
      return true;
    }
  }
  return false;
}

void ELEM::PrintSurfaceIndicator(int elm_id_to_print) const {
  for (SurfaceList_type::const_iterator it = surface_indicator_.begin();
       it != surface_indicator_.end(); it++) {
    printf("Elm%d  Surface%d:", elm_id_to_print, it->surface_id());
    for (unsigned int i = 0; i < SurfaceIndicator::MAX_SURFACE_INDICATORS;
         i++) {
      if (it->has_indicator(i)) {
        printf("  %d", i);
      }
    }
    printf("\n");
  }
}

double ELEM::GetMeasure(const GRID* p_grid) const {
  throw NotImplementedException() << "GetMeasure not implemented for elmType: "
                                  << elmType();
}

bool ELEM::IsInnerPoint(const GRID* grid, const ZEROPTV& W) const {
  // For each surface (face) of this element, we test if the point is
  // "below" the the surface. If the point is "below" every surface,
  // the point must be inside the element.

  // This implementation assumes the element is convex.
  // May give false negatives if the polygon/polyhedron is concave.
  // http://math.stackexchange.com/a/7934

  /*for (ElemNodeID i = 0; i < n_nodes(); i++) {
    if ((grid->GetNode(ElemToLocalNodeID(i))->location() - W).norm() < 1e-16)
      return true;
  }*/

  const int n_surfs = GetSurfaceCount();
  const int row_len = GetSurfaceCheckArrayRowLength();
  const int* surf_arr = GetSurfaceCheckArray();

  for (int surf_idx = 0; surf_idx < n_surfs; surf_idx++) {
    const int surface_id = surf_arr[surf_idx * row_len];
    const ElemNodeID* nodes = surf_arr + surf_idx * row_len + 1;
    ZEROPTV normal = CalculateNormal(grid, surface_id);
    ZEROPTV vec = W - grid->GetNode(ElemToLocalNodeID(nodes[0]))->location();
    double dot = normal.innerProduct(vec);
    if (dot > 1e-14)  // epsilon to treat 0 to 1e-14 as "on plane"
      return false;
  }

  return true;
}

void ELEM::GenSurfaceIndicator(const GRID* pGrid,
                               const ZEROARRAY<int>& indicatorArray) {
  const int* check_array = GetSurfaceCheckArray();
  const int n_surfaces = GetSurfaceCount();
  const int nodes_per_surface = GetNodesPerSurface();
  const int row_len = GetSurfaceCheckArrayRowLength();
  surface_indicator_.clear();
  for (int i = 0; i < n_surfaces; i++) {
    const int idx = i * row_len;
    const int surfaceID = check_array[idx];
    SurfaceIndicator newIndicator(surfaceID);
    for (int k = 0; k < indicatorArray.size(); k++) {
      const int indicator = indicatorArray(k);
      int allInd = 1;
      for (int j = 1; j <= nodes_per_surface; j++) {
        ElemNodeID a = check_array[idx + j];
        LocalNodeID A = ElemToLocalNodeID(a);
        if (!pGrid->BoNode(A, indicator)) {
          allInd = 0;
          break;
        }
      }
      if (allInd) {
        newIndicator.add_indicator(indicator);
      }
    }
    if (newIndicator.num_indicators() > 0) {
      this->surface_indicator_.push_back(newIndicator);
    }
  }

  CalculateSurfaceNormals(pGrid);
}

const int* ELEM::GetNodesInSurface(int surface_id) const {
  const int* arr = GetSurfaceCheckArray();
  for (int i = 0; i < GetSurfaceCount(); i++) {
    if (*arr == surface_id) {
      return arr + 1;
    }
    arr += GetSurfaceCheckArrayRowLength();
  }

  throw TALYException() << "Unknown surface_id '" << surface_id << "' for "
                        << "element type '" << elmType() << "'";
}

void ELEM::CalculateSurfaceNormals(const GRID* grid) {
  for (SurfaceList_type::iterator it = surface_indicator_.begin();
       it != surface_indicator_.end(); it++) {
    it->set_normal(CalculateNormal(grid, it->surface_id()));
  }
}

ZEROPTV ELEM::CalculateCenter(const GRID* grid) const {
  ZEROPTV center(0, 0, 0, 0);
  for (ElemNodeID i = 0; i < n_nodes(); i++) {
    center += grid->GetNode(ElemToLocalNodeID(i))->location();
  }
  center *= 1.0 / n_nodes();
  return center;
}

double ELEM::CalculateRadius(const GRID* grid) const {
  const ZEROPTV center = CalculateCenter(grid);
  return CalculateRadius(grid, center);
}

double ELEM::CalculateRadius(const GRID* grid, const ZEROPTV& center) const {
  double max_radius = 0.0;
  for (ElemNodeID i = 0; i < n_nodes(); i++) {
    const NODE* node = grid->GetNode(ElemToLocalNodeID(i));
    const ZEROPTV p = node->location() - center;
#ifdef ENABLE_4D
    const double r = sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z() + p.t()*p.t());
#else
    const double r = sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z());
#endif
    max_radius = std::max(max_radius, r);
  }

  return max_radius;
}

}  // namespace TALYFEMLIB
