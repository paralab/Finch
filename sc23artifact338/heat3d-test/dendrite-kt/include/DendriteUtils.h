//
// Created by maksbh on 9/16/19.
//

#ifndef DENDRITEKT_DENDRITEUTILS_H
#define DENDRITEKT_DENDRITEUTILS_H

#include <iostream>
#include "DataTypes.h"
#include <vector>
#include <oda.h>
#include <talyfem/grid/femelm.h>

/**
 * @brief : Called at initialization.
 **/
void dendrite_init(int argc, char **argv);


/**
 * @brief creates a regular DA (Cartesian Grid).
 * @param level level of the regular DA
 * @param eleOrder element order 1/2/4.
 * @param max_depth max Depth of DA. default set to 25
 * @return regular DA.
 */
ot::DA<DIM> * createRegularDA(std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> & treePart,const DENDRITE_UINT level, const  DENDRITE_UINT eleOrder, const DENDRITE_UINT max_depth = 25);

/**
 *
 * @param[out] subDA SubDA
 * @param extent The extent in 2^{x,y,z}
 * @param level corresponding to max(x,y,z)
 * @param eleOrder element order
 * @param maxDepth maximum Depth
 *
 */
[[deprecated]]
void createRegularSubDA(ot::DA<DIM> & subDA,std::vector<TREENODE> & treenode,const std::array<unsigned int,DIM> & extent, const DENDRITE_UINT level, const DENDRITE_UINT eleOrder, const DENDRITE_UINT maxDepth = 25) ;



DA * createSubDA(DistTREE & distTREE,const PhysicalDomainDecider & funtion, const DENDRITE_UINT level, const DENDRITE_UINT eleOrder,const double sfcTol = 0.1, const DENDRITE_UINT maxDepth = 25);
/**
 * @brief finalize the program. Safely deletes the memory.
 * @param da octDA.
 *
 * */
void dendrite_finalize(ot::DA<DIM> * da);

void dendrite_finalize();
/**
 *
 * @brief Interpolates the value at integration point.
 * @param [in]fe : fe object
 * @param [in] ndof:  number of degrees of freedom corresponding to value
 * @param [in] value: value vector for the corresponding element in the order of (uvw,uvw,...)
 * @param [out] val_c: value FEM  at the integration point. Must be allocated outside the function
 * Usage:
 * while (fe.next_itg_pt()) {
 *  calcValueFEM(fe, ndof, val_taly, val_c);
 * }
 */

void calcValueFEM(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT ndof, const DENDRITE_REAL *value, DENDRITE_REAL *val_c) ;

void coordsToZeroptv(const DENDRITE_REAL *coords , std::vector<TALYFEMLIB::ZEROPTV >&node, const DENDRITE_UINT ele_order, bool isAllocated);

void calcValueDerivativeFEM(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT ndof, const DENDRITE_REAL *value, DENDRITE_REAL *val_c);
void petscPrintArray(const Vec & v);

void GetLocalPtv(TALYFEMLIB::FEMElm fe, const TALYFEMLIB::ZEROPTV &ptvg, TALYFEMLIB::ZEROPTV &ptvl);
#ifdef IBM
DENDRITE_REAL getNormalDistance(const TALYFEMLIB::FEMElm & fe,
                                const DENDRITE_REAL * gaussPointPosition,
                                const DENDRITE_REAL * gaussPointNormal,
                                const TALYFEMLIB::ZEROPTV & h,
                                const DENDRITE_REAL threshold = 1E-2);
#endif
#endif //DENDRITEKT_DENDRITEUTILS_H
