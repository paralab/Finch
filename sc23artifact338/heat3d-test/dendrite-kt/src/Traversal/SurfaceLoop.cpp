//
// Created by maksbh on 9/8/20.
//

#include <Traversal/SurfaceLoop.h>

SurfaceLoop::SurfaceLoop(DA *octDA, const std::vector<TREENODE> &treePart, const DomainExtents &domainExtents,
                         SubDomainBoundary *subDomainBoundary)
    : Traversal(octDA, treePart, domainExtents), eleOrder_(octDA->getElementOrder()),
      m_subdomainBoundary(subDomainBoundary) {

  init();

}

SurfaceLoop::SurfaceLoop(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo &v,
                         const DomainExtents &domainExtents,
                         SubDomainBoundary *subDomainBoundary)
    : Traversal(octDA, treePart, v, domainExtents), eleOrder_(octDA->getElementOrder()),
      m_subdomainBoundary(subDomainBoundary) {

  init();
  surfaceValues_ = new PetscScalar[v.ndof*octDA->getNumNodesPerElement()];
}

void SurfaceLoop::traverse() {
  Traversal::traverse();
}

void SurfaceLoop::traverseOperation(TALYFEMLIB::FEMElm &fe) {
  if (this->m_BoundaryOctant) {
    generateSurfaceFlags(this->m_coords);
    const int *surf_arr = m_elem->GetSurfaceCheckArray();
    const int surf_row_len = m_elem->GetSurfaceCheckArrayRowLength();

    for (int i = 0; i < 2 * DIM; i++) {
      const DENDRITE_UINT &boundaryType = isBoundaryFace_[i].boundaryType;
      int surf_id = surf_arr[i * surf_row_len];  // id determines which nodes on the surface
      if (boundaryType != -1) {
        TALYFEMLIB::SurfaceIndicator surf(surf_id);
        TALYFEMLIB::ZEROPTV  normal = m_elem->CalculateNormal(&m_grid, surf_id);
        if (boundaryType >= BoundaryTypes::WALL::MAX_WALL_TYPE_BOUNDARY) {
          surf.set_normal(normal * -1);
        } else {
          surf.set_normal(normal);
        }
        fe.refill_surface(m_elem,&surf,0);
        for(int n = 0; n < m_surfaceCoords.size();n++) {
          for(int d = 0; d < DIM; d++) {
            m_surfaceCoords[n].data()[d] = this->m_coords[surf_arr[i*surf_row_len + n + 1]*DIM + d];
          }
        }
        performSurfaceOperation(fe,m_surfaceCoords,isBoundaryFace_[i]);
      }
    }
  }
}

void SurfaceLoop::traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar * values) {
  if (this->m_BoundaryOctant) {
    generateSurfaceFlags(this->m_coords);
    const int *surf_arr = m_elem->GetSurfaceCheckArray();
    const int surf_row_len = m_elem->GetSurfaceCheckArrayRowLength();

    for (int i = 0; i < 2 * DIM; i++) {
      const DENDRITE_UINT &boundaryType = isBoundaryFace_[i].boundaryType;
      int surf_id = surf_arr[i * surf_row_len];  // id determines which nodes on the surface
      if (boundaryType != -1) {
        TALYFEMLIB::SurfaceIndicator surf(surf_id);
        TALYFEMLIB::ZEROPTV  normal = m_elem->CalculateNormal(&m_grid, surf_id);
        /// Here we are assuming that we are in the reference of a solid body.
        /// To be consistent with IBM framework, the normal direction should be pointing outward with respect to solid.
        surf.set_normal(normal * -1);
        fe.refill_surface(m_elem,&surf,0);
        for(int n = 0; n < m_surfaceCoords.size();n++) {
          for(int d = 0; d < DIM; d++) {
            m_surfaceCoords[n].data()[d] = this->m_coords[surf_arr[i*surf_row_len + n + 1]*DIM + d];
          }
          for(int dof = 0; dof < this->getNdof(); dof ++) {
            surfaceValues_[n*this->getNdof()+dof] = values[surf_arr[i*surf_row_len + n + 1]*getNdof() + dof];
          }

        }
        performSurfaceOperation(fe,m_surfaceCoords,isBoundaryFace_[i],values);
      }
    }
  }
}

void SurfaceLoop::init() {
  m_surfaceCoords.resize(this->m_elem->GetNodesPerSurface());

}

void SurfaceLoop::generateSurfaceFlags(double *physCoords) {
  static constexpr unsigned numPoints = 1u << DIM;

  if (this->m_octDA->getElementOrder() > 2) {
    throw std::runtime_error("Not implmeneted");
  }
#if (DIM == 3)
  static constexpr DENDRITE_UINT cornerMap[2][numPoints]{{0, 1, 2, 3, 4,  5,  6,  7},
                                                         {0, 2, 6, 8, 18, 20, 24, 26}};
  static constexpr DENDRITE_UINT numNodesPerFace = 4;
  static constexpr DENDRITE_UINT faceID[2 * DIM][4]
      {
          {0, 2, 4, 6}, // Left
          {1, 3, 5, 7}, // Right
          {0, 1, 4, 5}, // Bottom
          {2, 3, 6, 7}, // Top
          {0, 1, 2, 3}, // Back
          {4, 5, 6, 7},  // Front
      };
#elif(DIM == 2)
  static constexpr DENDRITE_UINT cornerMap[2][numPoints]{{0,1,2,3},{0,2,6,8}};
  static constexpr DENDRITE_UINT numNodesPerFace = 2;
  static constexpr DENDRITE_UINT faceID[2*DIM][2]
      {
          {2,0}, // Left
          {1,3}, // Right
          {0,1}, // Bottom
          {3,2}, // Top
      };

#else
  throw std::runtime_error("Not implemented\n");
#endif


  std::array<Point<DIM>, numPoints> arrayOfPoints;
  for (int i = 0; i < numPoints; i++) {
    arrayOfPoints[i] = Point<DIM>(&physCoords[cornerMap[eleOrder_ - 1][i] * DIM]);
  }

  std::array<DENDRITE_UINT, numPoints> ids{};
  std::array<std::bitset<MAX_BOUNDARY_TYPES>, numPoints> boundaryBits;
  for (int i = 0; i < numPoints; i++) {
    m_subdomainBoundary->generateBoundaryFlags(arrayOfPoints[i], ids[i]);
    boundaryBits[i] = m_subdomainBoundary->getBoundary();
  }
  BoundarySurface surface{static_cast<DENDRITE_UINT>(-1), static_cast<DENDRITE_UINT>(-1)};
  isBoundaryFace_.fill(surface);


  for (DENDRITE_UINT id = 0; id < 2 * DIM; id++) {
    for (DENDRITE_UINT i = 0; i < MAX_BOUNDARY_TYPES; i++) {

      bool check = boundaryBits[faceID[id][0]].test(i);
      DENDRITE_UINT objectID = ids[faceID[id][0]];

      for (DENDRITE_UINT face = 1; (face < numNodesPerFace) and (check == true); face++) {
        check = boundaryBits[faceID[id][face]].test(i) and (ids[faceID[id][face]] == objectID);
      }
      if (check == true) {
        isBoundaryFace_[id].boundaryType = i;
        isBoundaryFace_[id].id = objectID;
        break;
      }
    }
  }
}

SurfaceLoop::~SurfaceLoop() {
  if(this->m_traversalType == TRAVERSAL_TYPE::VALUES){
    delete [] surfaceValues_;
  }
}