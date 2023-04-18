//
// Created by maksbh on 9/18/19.
//

#ifndef DENDRITEKT_TALYMESH_H
#define DENDRITEKT_TALYMESH_H
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

/**
 * @brief generates the mesh to be used in the Integrands.
 * @tparam NodeData Node Data class
 */

template <typename NodeData>
class TalyMesh {
 public:
  TALYFEMLIB::ELEM *taly_elem;
  TALYFEMLIB::GRID grid;
  TALYFEMLIB::GridField<NodeData> field;

  TalyMesh(int order);

  virtual ~TalyMesh();
};


template <typename NodeData>
TalyMesh<NodeData>::TalyMesh(int order){
#if (DIM == 3)
  if(order == 1) {
    int node_id_array[8];
    grid.redimArrays(8, 1);
    for (int i = 0; i < 8; i++) {
      grid.node_array_[i] = new TALYFEMLIB::NODE();
      node_id_array[i] = i;
    }

    taly_elem = new TALYFEMLIB::ELEM3dHexahedral();
    grid.elm_array_[0] = taly_elem;


    taly_elem->redim(8, node_id_array);

    field.redimGrid(&grid);
    field.redimNodeData();
  }
  else if(order == 2){
    int node_id_array[27];
    grid.redimArrays(27, 1);
    for (int i = 0; i < 27; i++) {
      grid.node_array_[i] = new TALYFEMLIB::NODE();
      node_id_array[i] = i;
    }

    taly_elem = new TALYFEMLIB::ELEM3dHexahedral();
    grid.elm_array_[0] = taly_elem;


    taly_elem->redim(27, node_id_array);

    field.redimGrid(&grid);
    field.redimNodeData();
  }
#elif (DIM == 4)

  if(order == 1) {
    int node_id_array[16];
    grid.redimArrays(16, 1);
    for (int i = 0; i < 16; i++) {
      grid.node_array_[i] = new TALYFEMLIB::NODE();
      node_id_array[i] = i;
    }

    taly_elem = new TALYFEMLIB::ELEM4dTesseract();
    grid.elm_array_[0] = taly_elem;


    taly_elem->redim(16, node_id_array);

    field.redimGrid(&grid);
    field.redimNodeData();
  }
  else if(order == 2){

    int node_id_array[81];
    grid.redimArrays(81, 1);
    for (int i = 0; i < 81; i++) {
      grid.node_array_[i] = new TALYFEMLIB::NODE();
      node_id_array[i] = i;
    }

    taly_elem = new TALYFEMLIB::ELEM4dTesseract();
    grid.elm_array_[0] = taly_elem;


    taly_elem->redim(81, node_id_array);

    field.redimGrid(&grid);

    field.redimNodeData();
  }

#elif (DIM == 2)
  if(order == 1) {
    static constexpr DENDRITE_UINT numNodes = 4;
    int node_id_array[numNodes];
    grid.redimArrays(numNodes, 1);
    for (int i = 0; i < numNodes; i++) {
      grid.node_array_[i] = new TALYFEMLIB::NODE();
      node_id_array[i] = i;
    }

    taly_elem = new TALYFEMLIB::ELEM2dBox();
    grid.elm_array_[0] = taly_elem;


    taly_elem->redim(numNodes, node_id_array);

    field.redimGrid(&grid);
    field.redimNodeData();
  }
  else if(order == 2){
    static constexpr DENDRITE_UINT numNodes = 9;
    int node_id_array[numNodes];
    grid.redimArrays(numNodes, 1);
    for (int i = 0; i < numNodes; i++) {
      grid.node_array_[i] = new TALYFEMLIB::NODE();
      node_id_array[i] = i;
    }

    taly_elem = new TALYFEMLIB::ELEM2dBox();
    grid.elm_array_[0] = taly_elem;


    taly_elem->redim(numNodes, node_id_array);

    field.redimGrid(&grid);

    field.redimNodeData();
  }
#else
  throw  TALYFEMLIB::TALYException() << "Check your dimension\n";
#endif
}

template <typename NodeData>
TalyMesh<NodeData>::~TalyMesh() {
}

#endif //DENDRITEKT_TALYMESH_H
