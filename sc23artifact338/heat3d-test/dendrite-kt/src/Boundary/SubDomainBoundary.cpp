//
// Created by maksbh on 7/23/20.
//
#include <Boundary/SubDomainBoundary.h>

SubDomainBoundary::SubDomainBoundary(const SubDomain *subDomain, const DA *octDA, const DomainExtents &domainExtents)
        : m_subDomain(subDomain), DomainBoundary(octDA, domainExtents) {

}

void SubDomainBoundary::generateBoundaryFlags(const Point<DIM> &position, DENDRITE_UINT &id) {

    DomainBoundary::generateBoundaryFlags(position, id);
    if (isInsideCircle(position, id)) {
        this->m_boundary.set(BoundaryTypes::VOXEL::CIRCLE, true);
    }
    if (isInsideSphere(position, id)) {
        this->m_boundary.set(BoundaryTypes::VOXEL::SPHERE, true);
    }
    if (isInsideBox(position, id)) {
        this->m_boundary.set(BoundaryTypes::VOXEL::BOX, true);
    }
    if (isInsideGeometry(position, id)) {
        this->m_boundary.set(BoundaryTypes::VOXEL::GEOMETRY, true);
    }
}

void SubDomainBoundary::generateBoundaryFlags(const TALYFEMLIB::ZEROPTV &position, DENDRITE_UINT &id) {
    std::array<double, DIM> point;
    std::memcpy(point.data(), position.data(), sizeof(double) * DIM);
    generateBoundaryFlags(point, id);
}


bool SubDomainBoundary::isInsideCircle(const Point<DIM> &position, DENDRITE_UINT &id) const {
    DENDRITE_REAL coords[DIM];
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
        coords[dim] = position.x(dim);
    }
    for (DENDRITE_UINT i = 0; i < m_subDomain->getVoxelCircles().size(); i++) {
        const auto &circles = m_subDomain->getVoxelCircles()[i];
        if (circles.ifInside(coords)) {
            id = i;
            return true;
        }
    }
    return false;
}

bool SubDomainBoundary::isInsideSphere(const Point<DIM> &position, DENDRITE_UINT &id) const {
    DENDRITE_REAL coords[DIM];
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
        coords[dim] = position.x(dim);
    }
    for (DENDRITE_UINT i = 0; i < m_subDomain->getVoxelSpheres().size(); i++) {
        const auto &spheres = m_subDomain->getVoxelSpheres()[i];
        if (spheres.ifInside(coords)) {
            id = i;
            return true;
        }
    }
    return false;
}

bool SubDomainBoundary::isInsideBox(const Point<DIM> &position, DENDRITE_UINT &id) const {
    DENDRITE_REAL coords[DIM];
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
        coords[dim] = position.x(dim);
    }
    for (DENDRITE_UINT i = 0; i < m_subDomain->getVoxelBox().size(); i++) {
        const auto &box = m_subDomain->getVoxelBox()[i];
        if (box.ifInside(coords)) {
            id = i;
            return true;
        }
    }
    return false;
}

bool SubDomainBoundary::isInsideGeometry(const Point<DIM> &position, DENDRITE_UINT &id) const {
    DENDRITE_REAL coords[DIM];
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
        coords[dim] = position.x(dim);
    }
    for (DENDRITE_UINT i = 0; i < m_subDomain->getGeometry().size(); i++) {
        const auto &geom = m_subDomain->getGeometry()[i];
        if (geom->ifInside(coords)) {
            id = i;
            return true;
        }

    }
    return false;
}



