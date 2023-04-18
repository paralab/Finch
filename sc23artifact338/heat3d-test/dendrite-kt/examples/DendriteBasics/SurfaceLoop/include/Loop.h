//
// Created by maksbh on 10/15/20.
//

#ifndef DENDRITEKT_LOOP_H
#define DENDRITEKT_LOOP_H
#include <Traversal/SurfaceLoop.h>
class Loop: public SurfaceLoop{
    double length = 0;
public:
    Loop(DA * octDA,const std::vector<TREENODE> & treePart, const VecInfo & v, const DomainExtents & domain, SubDomainBoundary *subDomainBoundary);
    void performSurfaceOperation(TALYFEMLIB::FEMElm surfaceFe, const std::vector<TALYFEMLIB::ZEROPTV> &surfaceCoords,
                            const BoundarySurface &boundarySurface) override;
    void
    performSurfaceOperation(TALYFEMLIB::FEMElm & surfaceFe, const std::vector<TALYFEMLIB::ZEROPTV> &surfaceCoords,
                            const BoundarySurface &boundarySurface, const PetscScalar *values) override;

    void getTotalValue();
};

Loop::Loop(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo & v,const DomainExtents &domain,
           SubDomainBoundary *subDomainBoundary)
           :SurfaceLoop(octDA,treePart,v,domain,subDomainBoundary){
    this->traverse();
}

void
Loop::performSurfaceOperation(TALYFEMLIB::FEMElm surfaceFe, const std::vector<TALYFEMLIB::ZEROPTV> &surfaceCoords,
                                   const BoundarySurface &boundarySurface){

    if (boundarySurface.boundaryType == BoundaryTypes::X_MINUS) {
        while (surfaceFe.next_itg_pt()) {
            length += surfaceFe.detJxW();
        }
    }
}
void
Loop::performSurfaceOperation(TALYFEMLIB::FEMElm & surfaceFe, const std::vector<TALYFEMLIB::ZEROPTV> &surfaceCoords,
                        const BoundarySurface &boundarySurface, const PetscScalar *values){
    if (boundarySurface.boundaryType == BoundaryTypes::X_MINUS) {
        while (surfaceFe.next_itg_pt()) {
            length += surfaceFe.detJxW();
        }
    }
}
void Loop::getTotalValue() {
    double globalVal;
    MPI_Reduce(&length,&globalVal,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    std::cout << globalVal << "\n";
}
#endif //DENDRITEKT_LOOP_H
