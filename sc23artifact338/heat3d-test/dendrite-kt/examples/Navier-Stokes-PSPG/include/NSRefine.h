//
// Created by maksbh on 5/22/20.
//

#ifndef DENDRITEKT_NSREFINE_H
#define DENDRITEKT_NSREFINE_H
#include <Traversal/Refinement.h>
#include <NSInputData.h>
class NSRefine : public Refinement {
  DENDRITE_UINT boundaryRefinementLevel_;
 public:
/**
 *
 * @param da The oldDA that you want to refine
 * @param numLocalElements Number of localElements
 * @param domainInfo physical domain information
 * @param inputData inputData
 */
  NSRefine(DA *da, const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, const size_t numLocalElements, const DENDRITE_UINT boundaryRefinementLevel);
  /**
   * @brief The refinement flags per element by element.
   * @param fe
   * @param coords vector of coords
   * @return the flags for each elements
   */
  ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;
};

NSRefine::NSRefine(DA *da, const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, const size_t numLocalElements, const DENDRITE_UINT boundaryRefinementLevel)
    : Refinement(da, treePart,domainInfo) {
  boundaryRefinementLevel_ = boundaryRefinementLevel;
  this->initRefinement();
}

ot::OCT_FLAGS::Refine NSRefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) {
  if(this->m_BoundaryOctant and (this->m_level < boundaryRefinementLevel_)){
    return ot::OCT_FLAGS::Refine::OCT_REFINE;
  }
  else{
    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
  }

}

#endif //DENDRITEKT_NSREFINE_H
