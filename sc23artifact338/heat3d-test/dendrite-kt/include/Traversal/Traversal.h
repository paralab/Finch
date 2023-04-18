//
// Created by maksbh on 5/9/20.
//

#ifndef DENDRITEKT_TRAVERSAL_H
#define DENDRITEKT_TRAVERSAL_H

#include <oda.h>
#include <DataTypes.h>
#include <sfcTreeLoop_matvec_io.h>
#include <talyfem/grid/femelm.h>
#include <TalyDendroSync.h>
#include <PETSc/VecInfo.h>
#include <OctToPhysical.h>

enum TRAVERSAL_TYPE {
    COORDS = 0,
    VALUES = 1
};

class Traversal {

    TalyDendroSync sync_;


    int *node_id_array_;
    DENDRITE_UINT m_uiDof;


    DENDRITE_UINT eleOrder_;
    DENDRITE_UINT nPe_;

    VecInfo v_;

    const int * relativeOrder_ = nullptr;
    bool traverseByCoords();

    void constructTalyGrid(const double *coords);

    bool traverseByValues();
    void init();

  inline int getRelativeOrder(int localElemNumber) const {
    int relativeOrder = 0;
    if(relativeOrder_) {
      relativeOrder = (relativeOrder_[localElemNumber]);
    }
    return relativeOrder;
  }
protected:
    TALYFEMLIB::ELEM *m_elem;
    ot::DA<DIM> *m_octDA;
    TRAVERSAL_TYPE m_traversalType;
    DENDRITE_REAL * m_coords;
    DENDRITE_UINT m_level;
    OctToPhysical m_octToPhysical;
    bool m_BoundaryOctant = false;
    const std::vector<TREENODE> &m_treePart;
    TALYFEMLIB::GRID m_grid;

public:
    Traversal(ot::DA<DIM> *octDA, const std::vector<TREENODE> &treePart, const DomainExtents &domain);

    Traversal(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo &v, const DomainExtents &domain);

    ~Traversal();

    void setRelativeOrderVec(const int * relativeOrder){
      relativeOrder_ = relativeOrder;
    }





    const DENDRITE_UINT getNdof() const;

    bool traverse();

    virtual void traverseOperation(TALYFEMLIB::FEMElm &fe) {
      throw TALYFEMLIB::TALYException() << " You need to override " << __func__ << "\n";
    }

    virtual void traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values) {
      throw TALYFEMLIB::TALYException() << " You need to override " << __func__ << "\n";
    }

};


#endif //DENDRITEKT_TRAVERSAL_H
