//
// Created by maksbh on 9/8/20.
//

#ifndef DENDRITEKT_SURFACELOOP_H
#define DENDRITEKT_SURFACELOOP_H

#include <Traversal/Traversal.h>
#include <Boundary/SubDomainBoundary.h>

class SurfaceLoop : public Traversal {

    std::vector<TALYFEMLIB::ZEROPTV> m_surfaceCoords;
    const DENDRITE_UINT eleOrder_;
    std::array<BoundarySurface, 2 * DIM> isBoundaryFace_;
    PetscScalar * surfaceValues_;
    void init();

protected:
    SubDomainBoundary *m_subdomainBoundary;


public:
    SurfaceLoop(DA *octDA, const std::vector<TREENODE> &treePart, const DomainExtents &domainExtents,
                SubDomainBoundary *subDomainBoundary);

    SurfaceLoop(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo &v, const DomainExtents &domainExtents,
                SubDomainBoundary *subDomainBoundary);

    void traverseOperation(TALYFEMLIB::FEMElm &fe) override;
    void traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar * values) override;

    void traverse();

    ~SurfaceLoop();



    virtual void
    performSurfaceOperation(TALYFEMLIB::FEMElm & surfaceFe, const std::vector<TALYFEMLIB::ZEROPTV> &surfaceCoords,
                       const BoundarySurface &boundarySurface, const PetscScalar *values) {
      throw TALYFEMLIB::TALYException() << "You need to overload this" << __func__ << "\n";
    }

    virtual void
    performSurfaceOperation(TALYFEMLIB::FEMElm surfaceFe, const std::vector<TALYFEMLIB::ZEROPTV> &surfaceCoords, const BoundarySurface &boundarySurface) {
      throw TALYFEMLIB::TALYException() << "You need to overload this" << __func__ << "\n";
    }

    void generateSurfaceFlags(double *physCoords);
};

#endif //DENDRITEKT_SURFACELOOP_H
