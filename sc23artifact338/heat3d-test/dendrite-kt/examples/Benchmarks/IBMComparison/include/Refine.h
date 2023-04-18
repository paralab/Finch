//
// Created by maksbh on 5/29/21.
//

#ifndef DENDRITEKT_REFINE_H
#define DENDRITEKT_REFINE_H

#include <Traversal/Refinement.h>

enum CaseType : int{
    IBM = 0,
    CARVED = 1
};
class Refine: public Refinement{
  SubDomainBoundary * boundary_;
  const DENDRITE_UINT & maxLevel_;
  const CaseType caseType;
  const GEOMETRY::STL *stl_;
  const bool final;
public:
  Refine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, SubDomainBoundary * boundaryFlags, const DENDRITE_UINT & maxLevel,const CaseType _caseType, const GEOMETRY::STL * stl, const bool ifFinal);

  ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

};

Refine::Refine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo,SubDomainBoundary * boundary, const DENDRITE_UINT & maxLevel,const CaseType _caseType,const GEOMETRY::STL *stl,const bool ifFinal)
  :Refinement(da,treePart,domainInfo),boundary_(boundary),maxLevel_(maxLevel),caseType(_caseType),stl_(stl),final(ifFinal){
  this->initRefinement();
}

ot::OCT_FLAGS::Refine Refine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords){
  DENDRITE_UINT  id;
  if(caseType == CARVED) {
    for (auto &coord:coords) {
      boundary_->generateBoundaryFlags(coord, id);
      if ((boundary_->checkBoundaryType(BoundaryTypes::VOXEL::GEOMETRY)) and (this->m_level < maxLevel_)) {
        return ot::OCT_FLAGS::Refine::OCT_REFINE;
      }
    }
    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
  } else{

    int numInside(0);
    for(const auto & coord:coords){
      if(stl_->ifInside(coord.data())){
        numInside++;
      }
    }
    if((numInside != 0) and (numInside != coords.size()) and (this->m_level < maxLevel_)){
      return ot::OCT_FLAGS::Refine::OCT_REFINE;
    }
    if(final) {
      if ((numInside == coords.size()) and (this->m_level < maxLevel_ - 1)) {
        return ot::OCT_FLAGS::Refine::OCT_REFINE;
      }
    }

    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
  }

}
#endif //DENDRITEKT_REFINE_H
