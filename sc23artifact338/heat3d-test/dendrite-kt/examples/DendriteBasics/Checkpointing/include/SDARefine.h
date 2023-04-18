//
// Created by maksbh on 11/8/20.
//

#ifndef DENDRITEKT_SDAREFINE_H
#define DENDRITEKT_SDAREFINE_H
#include <Traversal/Refinement.h>

class SDARefine: public Refinement{
  const DENDRITE_UINT maxLevel_;
public:
    SDARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, DENDRITE_UINT maxLevel);

    ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

};

SDARefine::SDARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, const DENDRITE_UINT maxLevel)
:Refinement(da,treePart,domainInfo),maxLevel_(maxLevel){
    this->initRefinement();
}

ot::OCT_FLAGS::Refine SDARefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords){
    if((this->m_BoundaryOctant) and (this->m_level < maxLevel_)){
        return ot::OCT_FLAGS::Refine::OCT_REFINE;
    }
    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
}
#endif //DENDRITEKT_SDAREFINE_H
