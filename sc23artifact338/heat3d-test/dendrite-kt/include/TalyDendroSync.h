//
// Created by maksbh on 9/18/19.
//

#ifndef DENDRITEKT_TALYDENDROSYNC_H
#define DENDRITEKT_TALYDENDROSYNC_H

#include "DataTypes.h"
#include <type_traits>
#include <DataTypes.h>
#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/gridfield.h>

#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/gridfield.h>
#include <talyfem/grid/elem.h>
#include <talyfem/grid/elem_types/elem3dhexahedral.h>
#include <talyfem/grid/elem_types/elem2dbox.h>

#ifdef ENABLE_4D
#include <talyfem/grid/elem_types/elem4dtesseract.h>
#endif

#include <talyfem/grid/femelm.h>
#include <talyfem/data_structures/zeroarray.h>
#include <talyfem/data_structures/zeromatrix.h>

class TalyDendroSync {
public:
    template<const DENDRITE_UINT order, typename std::enable_if<order == 1>::type * = nullptr>
    void syncCoords(const DENDRITE_REAL *coords, TALYFEMLIB::GRID *grid);

    template<const DENDRITE_UINT order, typename std::enable_if<order == 2>::type * = nullptr>
    void syncCoords(const DENDRITE_REAL *coords_dendro, TALYFEMLIB::GRID *grid);

    template<class NodeData>
    void syncValues(const DENDRITE_REAL *in_dendro, TALYFEMLIB::GridField<NodeData> *gf,
                    const DENDRITE_UINT ndof, const DENDRITE_UINT sdof);

    template<class NodeData>
    void syncValues(const DENDRITE_REAL *in_dendro, TALYFEMLIB::GridField<NodeData> *gf,
                    const DENDRITE_UINT ndof, const std::vector<DENDRITE_UINT> &dofID);

    template<const DENDRITE_UINT order, typename std::enable_if<order == 1>::type * = nullptr>
    static void convertBeToTalyOrder(const TALYFEMLIB::ZEROARRAY<DENDRITE_REAL> &beDendro,
                                     TALYFEMLIB::ZEROARRAY<DENDRITE_REAL> &beTaly,
                                     const DENDRITE_UINT ndof);

    template<const DENDRITE_UINT order, typename std::enable_if<order == 1>::type * = nullptr>
    static void convertAeToTalyOrder(const TALYFEMLIB::ZEROMATRIX<DENDRITE_REAL> &AeDendro,
                                     TALYFEMLIB::ZEROMATRIX<DENDRITE_REAL> &AeTaly,
                                     const DENDRITE_UINT ndof);


};

template<const DENDRITE_UINT order, typename std::enable_if<order == 1>::type *>
void TalyDendroSync::syncCoords(const DENDRITE_REAL *coords, TALYFEMLIB::GRID *grid) {
#if (DIM == 3)
    for (unsigned int i = 0; i < 8; i++) {
        grid->node_array_[i]->setCoor(coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]);
    }
#elif (DIM == 4)
    for (unsigned int i = 0; i < 16; i++) {
  //    std::cout << i << " " << coords[i * DIM] << " " <<  coords[i * DIM + 1] << " " << coords[i * DIM + 2] << " "<<  coords[i * DIM + 3] << "\n";
      grid->node_array_[i]->setCoor(coords[i * DIM], coords[i * DIM + 1], coords[i * DIM + 2], coords[i * DIM + 3]);
    }
#elif (DIM == 2)
    for (unsigned int i = 0; i < 4; i++) {
      grid->node_array_[i]->setCoor(coords[i * DIM], coords[i * DIM + 1]);
    }
#else

    throw TALYFEMLIB::TALYException() << "Not implemented\n";

#endif
}

template<const DENDRITE_UINT order, typename std::enable_if<order == 2>::type *>
void TalyDendroSync::syncCoords(const DENDRITE_REAL *coords, TALYFEMLIB::GRID *grid) {
#if (DIM == 3)
    for (unsigned int i = 0; i < 27; i++) {
        grid->node_array_[i]->setCoor(coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]);
    }
#elif(DIM == 4)
    for (unsigned int i = 0; i < 81; i++) {
      grid->node_array_[i]->setCoor(coords[i * DIM], coords[i * DIM + 1], coords[i * DIM + 2], coords[i*DIM + 3]);
    }
#elif (DIM == 2)
    for (unsigned int i = 0; i < 9; i++) {
      grid->node_array_[i]->setCoor(coords[i * DIM], coords[i * DIM + 1]);
    }
#else
    throw TALYFEMLIB::TALYException() << "Not implemented\n";
#endif

}

template<class NodeData>
void TalyDendroSync::syncValues(const DENDRITE_REAL *in_dendro, TALYFEMLIB::GridField<NodeData> *gf,
                                const DENDRITE_UINT ndof, const DENDRITE_UINT sdof) {
    for (unsigned int dof = 0; dof < ndof; dof++) {
        gf->NodeDataFromArray(&in_dendro[dof], sdof + dof, ndof);
    }
}

template<class NodeData>
void TalyDendroSync::syncValues(const DENDRITE_REAL *in_dendro, TALYFEMLIB::GridField<NodeData> *gf,
                                const DENDRITE_UINT ndof, const std::vector<DENDRITE_UINT> &dofID) {
    for (unsigned int dof = 0; dof < ndof; dof++) {
        gf->NodeDataFromArray(&in_dendro[dof], dofID[dof], ndof);
    }
}

template<const DENDRITE_UINT order, typename std::enable_if<order == 1>::type *>
void TalyDendroSync::convertBeToTalyOrder(const TALYFEMLIB::ZEROARRAY<DENDRITE_REAL> &beDendro,
                                          TALYFEMLIB::ZEROARRAY<DENDRITE_REAL> &beTaly, const DENDRITE_UINT ndof) {

#if (DIM == 2)
    static constexpr DENDRITE_UINT mapOrder[]{0,1,3,2};
    for(DENDRITE_UINT dof = 0; dof < ndof; dof++){
        for(DENDRITE_UINT a = 0; a < 4; a++){
           beTaly(a*ndof+dof) = beDendro(mapOrder[a]*ndof+dof);
        }
    }

#elif (DIM == 3)
    throw std::runtime_error("Not implemented\n");
#endif

}
template<const DENDRITE_UINT order, typename std::enable_if<order == 1>::type*>
void TalyDendroSync::convertAeToTalyOrder(const TALYFEMLIB::ZEROMATRIX<DENDRITE_REAL> &AeDendro,
                                 TALYFEMLIB::ZEROMATRIX<DENDRITE_REAL> &AeTaly,
                                 const DENDRITE_UINT ndof){
#if (DIM == 2)
    static constexpr DENDRITE_UINT mapOrder[]{0,1,3,2};
    for(DENDRITE_UINT dofa = 0; dofa < ndof; dofa++) {
        for (DENDRITE_UINT dofb = 0; dofb < ndof; dofb++) {
            for (DENDRITE_UINT a = 0; a < 4; a++) {
                for (DENDRITE_UINT b = 0; b < 4; b++) {
                    AeTaly(a * ndof + dofa,b * ndof + dofb) = AeDendro(mapOrder[a] * ndof + dofa,mapOrder[b] * ndof + dofb);
                }
            }
        }
    }
#elif (DIM == 3)
    throw std::runtime_error("Not implemented\n");
#endif
}

#endif //DENDRITEKT_TALYDENDROSYNC_H
