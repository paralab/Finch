//
// Created by maksbh on 12/4/18.
//
#include <PETSc/VecInfo.h>

    VecInfo::VecInfo() {
        ndof = 0;
        v = NULL;
        nodeDataIndex = -1;
        placeholder = PLACEHOLDER_NONE;
    }

    VecInfo::VecInfo(Vec _v, DENDRITE_UINT _ndof, DENDRITE_UINT _nodeDataIndex, VecPlaceholder _placeholder) {
        ndof = _ndof;
        v = _v;
//        VecGetArrayRead(v, &val);
        nodeDataIndex = _nodeDataIndex;
        placeholder = _placeholder;
    }

    VecInfo::VecInfo(VecPlaceholder p, DENDRITE_UINT _ndof, DENDRITE_UINT _nodeDataIndex) {
        ndof = _ndof;
        v = NULL;
        val = NULL;
        nodeDataIndex = _nodeDataIndex;
        placeholder = p;
    }

    VecInfo::VecInfo(VecPlaceholder p, Vec & _v, DENDRITE_UINT _ndof, DENDRITE_UINT _nodeDataIndex) {

        ndof = _ndof;
        v = _v;
        nodeDataIndex = _nodeDataIndex;
        placeholder = p;
    }
