//
// Created by maksbh on 7/23/20.
//

#ifndef DENDRITEKT_SUBDOMAINBOUNDARY_H
#define DENDRITEKT_SUBDOMAINBOUNDARY_H

#include <Boundary/DomainBoundary.h>
#include <SubDA/SubDomain.h>

struct BoundarySurface{
    DENDRITE_UINT boundaryType;
    DENDRITE_UINT id;
};

class SubDomainBoundary:public DomainBoundary{

 protected:
  const SubDomain * m_subDomain;

  bool isInsideCircle(const Point<DIM> & position, DENDRITE_UINT & id) const;
  bool isInsideSphere(const Point<DIM> & position, DENDRITE_UINT & id) const;
  bool isInsideBox(const Point<DIM> & position, DENDRITE_UINT & id) const ;
  bool isInsideGeometry(const Point<DIM> & position, DENDRITE_UINT & id) const ;
 public:
  SubDomainBoundary(const SubDomain *subDomain,const DA * octDA, const DomainExtents & domainExtents);

  void generateBoundaryFlags(const Point<DIM> & position, DENDRITE_UINT & id) override;
  void generateBoundaryFlags(const TALYFEMLIB::ZEROPTV & position, DENDRITE_UINT & id) override;

};
#endif //DENDRITEKT_SUBDOMAINBOUNDARY_H
