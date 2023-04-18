//
// Created by maksbh on 8/27/20.
//

#include <Traversal/DomainBounds.h>

DomainBounds::DomainBounds(DA *octDA, const std::vector<TREENODE> & treePart,const DomainExtents &domain)
:Traversal(octDA,treePart,domain){
    std::memcpy(max_,domain.physicalDADomain.min.data(),sizeof(DENDRITE_REAL)*DIM);
    std::memcpy(min_,domain.physicalDADomain.max.data(),sizeof(DENDRITE_REAL)*DIM);
    this->traverse();
}

void DomainBounds::traverseOperation(TALYFEMLIB::FEMElm & fe){
    const DENDRITE_UINT npe = this->m_octDA->getNumNodesPerElement();
    if(this->m_BoundaryOctant){
        for(DENDRITE_UINT d = 0; d < DIM; d++){
            if(m_coords[(npe - 1)*DIM+d] > max_[d]){
                max_[d] = m_coords[(npe - 1)*DIM+d];
            }
            if(m_coords[d] < min_[d]){
                min_[d] = m_coords[d];
            }
        }
    }
}
void DomainBounds::updateDomain(DomainExtents &domain) const {
    DENDRITE_REAL globalMax[DIM];
    MPI_Allreduce(max_,globalMax,DIM,MPI_DOUBLE,MPI_MAX, MPI_COMM_WORLD);
    std::memcpy(domain.physicalDADomain.max.data(),globalMax,sizeof(DENDRITE_REAL)*DIM);
    DENDRITE_REAL globalMin[DIM];
    MPI_Allreduce(min_,globalMin,DIM,MPI_DOUBLE,MPI_MIN, MPI_COMM_WORLD);
    std::memcpy(domain.physicalDADomain.min.data(),globalMin,sizeof(DENDRITE_REAL)*DIM);

}

