/**
 * @file:matvec.h
 * @author: Masado Ishii  --  UofU SoC,
 * @date: 2019-03-15
 * @brief: Variations of the TreeSort algorithm (tsort.h) for mesh-free matvec. Based on Dendro-5.0.
 */

#ifndef DENDRO_KT_MATVEC_H
#define DENDRO_KT_MATVEC_H

#include "treeNode.h"
#include "tsort.h"
#include "nsort.h"

namespace fem {

  template <typename T, typename da, unsigned int dim>
  struct SFC_Matvec
  {
    /**
     * @brief Finds the required buffer size at each level by bucketing points.
     * @note Also shuffles the parallel companion array. This is needed to get the shuffle map.
     */
    static ot::RankI countSubtreeSizes(ot::TNPoint<T,dim> *points, ot::RankI *companions,
        ot::RankI begin, ot::RankI end,
        ot::LevI subtreeRootLev,
        ot::RotI pRot,
        int order,
        std::vector<ot::RankI> &outSubtreeSizes);




  };

}//namespace fem

#include "matvecPreallocation.tcc"

#endif//DENDRO_KT_MATVEC_H
