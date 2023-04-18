//
// Created by maksbh on 7/23/20.
//

#ifndef DENDRITEKT_BOUNDARY_H
#define DENDRITEKT_BOUNDARY_H
#include <Boundary/BoundaryDataTypes.h>

class DomainBoundary{
 protected:
  const DA * m_octDA;
  const DomainExtents & m_domainExtents;
  std::bitset<MAX_BOUNDARY_TYPES> m_boundary;

 public:
  DomainBoundary(const DA * octDA, const DomainExtents & domainExtents);

  virtual void generateBoundaryFlags(const Point<DIM> & position, DENDRITE_UINT & id);

  virtual void generateBoundaryFlags(const TALYFEMLIB::ZEROPTV & position, DENDRITE_UINT & id);

  bool checkBoundaryType(const DENDRITE_UINT & wallID) const;

  inline const std::bitset<MAX_BOUNDARY_TYPES> & getBoundary() const{
    return m_boundary;
  }


};
#endif //DENDRITEKT_BOUNDARY_H
