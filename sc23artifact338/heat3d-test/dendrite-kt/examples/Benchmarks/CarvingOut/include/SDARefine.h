//
// Created by maksbh on 1/18/21.
//

#ifndef DENDRITEKT_SDAREFINE_H
#define DENDRITEKT_SDAREFINE_H
#include <Traversal/Refinement.h>
#include <Boundary/SubDomainBoundary.h>
class SDARefine: public Refinement{
  SubDomainBoundary * boundary_;
  const DENDRITE_UINT & maxLevel_;
public:
  SDARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, SubDomainBoundary * boundaryFlags, const DENDRITE_UINT & maxLevel);

  ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

};

SDARefine::SDARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo,SubDomainBoundary * boundary, const DENDRITE_UINT & maxLevel)
  :Refinement(da,treePart,domainInfo),boundary_(boundary),maxLevel_(maxLevel){
  this->initRefinement();
}

ot::OCT_FLAGS::Refine SDARefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords){
  DENDRITE_UINT  id;

  for(auto & coord:coords){
    boundary_->generateBoundaryFlags(coord,id);
    if((boundary_->checkBoundaryType(BoundaryTypes::VOXEL::GEOMETRY)) and (this->m_level < maxLevel_)){
      return ot::OCT_FLAGS::Refine::OCT_REFINE;
    }
  }
  return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
}
#endif //DENDRITEKT_SDAREFINE_H
