//
// Created by maksbh on 12/4/18.
//


#ifndef DENDRITE2_0_VECINFO_H
#define DENDRITE2_0_VECINFO_H
#include <cstdio>
#include <DataTypes.h>
#include <petscvec.h>

enum VecPlaceholder {
    PLACEHOLDER_NONE = 0,
    PLACEHOLDER_GUESS  // "guess" solution given by SNES (FormFunction)
};

class VecInfo{
public:
    Vec v;
    const DENDRITE_REAL * val;
    DENDRITE_UINT ndof;
    DENDRITE_UINT nodeDataIndex;
    VecPlaceholder placeholder;
    VecInfo();
    void cleanup(){
        VecDestroy(&v);
        ndof = -1;
        nodeDataIndex = -1;
    }
    VecInfo(Vec v, DENDRITE_UINT _ndof, DENDRITE_UINT _nodeDataIndex, VecPlaceholder _placeholder = PLACEHOLDER_NONE);
    VecInfo(VecPlaceholder p, Vec & v, DENDRITE_UINT _ndof, DENDRITE_UINT _nodeDataIndex);
    VecInfo(VecPlaceholder p, DENDRITE_UINT _ndof, DENDRITE_UINT _nodeDataIndex);
//    void update_buffer(DENDRITE_REAL * new_val);

};







#endif //DENDRITE2_0_VECINFO_H
