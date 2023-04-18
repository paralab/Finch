//
// Created by maksbh on 8/27/20.
//

#ifndef LINNS_PSPG_PNP_DOMAINBOUNDS_H
#define LINNS_PSPG_PNP_DOMAINBOUNDS_H
#include <Traversal/Traversal.h>
class DomainBounds : public Traversal{
    DENDRITE_REAL max_[DIM];
    DENDRITE_REAL min_[DIM];

public:
    DomainBounds(DA * octDA,const std::vector<TREENODE> & treePart,const DomainExtents & domain);
    virtual void traverseOperation(TALYFEMLIB::FEMElm & fe) override;
    void updateDomain(DomainExtents & domain) const;
};
#endif //LINNS_PSPG_PNP_DOMAINBOUNDS_H
