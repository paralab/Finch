//
// Created by maksbh on 6/18/21.
//

#ifndef DENDRITEKT_SDAREFINE_H
#define DENDRITEKT_SDAREFINE_H
#include <Traversal/Refinement.h>

class Refine: public Refinement{
  const DENDRITE_UINT boundaryLevel;
  const DENDRITE_UINT baseLevel;
  const bool allRefine_ = false;
public:
  Refine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, const DENDRITE_UINT blevel, const DENDRITE_UINT baseLevel, const bool allRefine = false);

  ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

};

Refine::Refine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo,const DENDRITE_UINT blevel, const DENDRITE_UINT baselevel, const bool allRefine)
  :Refinement(da,treePart,domainInfo),boundaryLevel(blevel),baseLevel(baselevel),allRefine_(allRefine){
  this->initRefinement();
}

ot::OCT_FLAGS::Refine Refine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords){

  if(allRefine_){
    return ot::OCT_FLAGS::Refine::OCT_REFINE;
  }
  if((this->m_level < boundaryLevel) and (this->m_BoundaryOctant)){
    return ot::OCT_FLAGS::Refine::OCT_REFINE;
  }
  return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
}

#endif //DENDRITEKT_SDAREFINE_H
