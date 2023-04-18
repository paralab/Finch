//
// Created by maksbh on 7/23/20.
//

#include <Boundary/DomainBoundary.h>

DomainBoundary::DomainBoundary(const DA *octDA, const DomainExtents &domainExtents)
:m_octDA(octDA),m_domainExtents(domainExtents){

}

void DomainBoundary::generateBoundaryFlags(const Point<DIM> &position,  DENDRITE_UINT & id) {
  m_boundary.reset();
  if(FEQUALS(position.x(0),m_domainExtents.physicalDADomain.min[0])){
    m_boundary.set(BoundaryTypes::X_MINUS,true);
  }
  if(FEQUALS(position.x(0),m_domainExtents.physicalDADomain.max[0])){
    m_boundary.set(BoundaryTypes::X_PLUS,true);
  }
  if(FEQUALS(position.x(1),m_domainExtents.physicalDADomain.min[1])){
    m_boundary.set(BoundaryTypes::Y_MINUS,true);
  }
  if (FEQUALS(position.x(1),m_domainExtents.physicalDADomain.max[1])){
    m_boundary.set(BoundaryTypes::Y_PLUS,true);
  }
#if(DIM == 3)
  if (FEQUALS(position.x(2),m_domainExtents.physicalDADomain.min[2])){
    m_boundary.set(BoundaryTypes::Z_MINUS,true);
  }
  if (FEQUALS(position.x(2),m_domainExtents.physicalDADomain.max[2])){
    m_boundary.set(BoundaryTypes::Z_PLUS,true);
  }
#endif
  id = 0;
}

void DomainBoundary::generateBoundaryFlags(const TALYFEMLIB::ZEROPTV &position,  DENDRITE_UINT & id)  {
  std::array<double,DIM> point;
  std::memcpy(point.data(),position.data(),sizeof(double)*DIM);
  generateBoundaryFlags(Point<DIM>(point),id);
  id = 0;
}

bool DomainBoundary::checkBoundaryType(const DENDRITE_UINT & wallID) const {
  if(m_boundary.test(wallID)){
    return true;
  }
  else{
    return false;
  }
}