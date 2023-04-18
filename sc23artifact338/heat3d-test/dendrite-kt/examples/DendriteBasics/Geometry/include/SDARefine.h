//
// Created by maksbh on 7/18/20.
//

#ifndef DENDRITEKT_SDAREFINE_H
#define DENDRITEKT_SDAREFINE_H

#include <Traversal/Refinement.h>
#include <SubDA/Voxel.h>
#include <Boundary/SubDomainBoundary.h>

class SDARefine : public Refinement {
  SubDomainBoundary * boundary_;
  const DENDRITE_UINT  maxLevel_;
 public:
  SDARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, SubDomainBoundary * domainBoundary, const DENDRITE_UINT maxLevel);
  ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;


};

SDARefine::SDARefine(DA *da, const std::vector<TREENODE> & treePart,const DomainExtents &domainInfo, SubDomainBoundary * domainBoundary, const DENDRITE_UINT maxLevel)
    : Refinement(da, treePart,domainInfo),boundary_(domainBoundary),maxLevel_(maxLevel){

  this->traverse();
}

ot::OCT_FLAGS::Refine
SDARefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) {
  ot::OCT_FLAGS::Refine refineFlag = ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
  DENDRITE_UINT id;
  if((this->m_BoundaryOctant) and (this->m_level < maxLevel_)){
    for(const auto & coord: coords){
       boundary_->generateBoundaryFlags(coord,id);
       bool isSTLBoundary = boundary_->checkBoundaryType(BoundaryTypes::GEOMETRY);
       if(isSTLBoundary){
         refineFlag = ot::OCT_FLAGS::Refine::OCT_REFINE;
       }
    }
  }
  return refineFlag;
}


#endif //DENDRITEKT_SDAREFINE_H
