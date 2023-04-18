//
// Created by maksbh on 7/19/20.
//

#ifndef DENDRITEKT_SUBDOMAIN_H
#define DENDRITEKT_SUBDOMAIN_H

#include <SubDA/Voxel.h>
#include <OctToPhysical.h>
class SubDomain: public VOXEL::Voxel{

  enum CreationStage : bool{
    INITIAL = true,
    REFINE = false
  };
  enum DIR{
    X = 0,
    Y = 1,
#if (DIM == 3)
    Z = 2
#endif
  };
  enum ReasonForIn:DENDRITE_UINT {
      DOMAIN_BOUNDARY = 0,
      CIRCLE = 1,
      SPHERE = 2,
      BOX = 3,
      GEOMETRY = 4,
      MAX = 5
  };
  OctToPhysical octToPhysical_;


  CreationStage creationStage_ = CreationStage::INITIAL;
 protected:



  const DomainExtents & m_domainExtents;

  bool m_carvedDir[DIM];

  ibm::Partition octantToRetain(const double * physCoords, const double * scalingFactor) const;

  void updatePhysicalDomain(DA * octDA,const std::vector<TREENODE> & treePart,DomainExtents & domainExtent) const;

  ibm::Partition classifyNodes(const DENDRITE_REAL * physCoords, ReasonForIn * reasonForIn, const DENDRITE_UINT & nodeID) const;


public:
  SubDomain(const DomainExtents & domainExtents, bool resumeFromCheckpoint = false);
  [[deprecated]]
  inline const DomainInfo & domain() const{
    return m_domainExtents.fullDADomain;
  }

  inline const DomainExtents & domainExtents() const{
    return m_domainExtents;
  }

  ibm::Partition functionToRetain(const double *octCoords, const double scale) const ;



  void finalize(DA * octDA,const std::vector<TREENODE> & treePart,DomainExtents & domainExtents);

};
#endif //DENDRITEKT_SUBDOMAIN_H
