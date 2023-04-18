/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //
#ifndef GRID_GRID_TYPES_H_
#define GRID_GRID_TYPES_H_

#include <talyfem/grid/grid_types/grid1d.h>
#include <talyfem/grid/grid_types/grid2d.h>
#include <talyfem/grid/grid_types/gridbox2d.h>
#include <talyfem/grid/grid_types/gridbox3d.h>
#include <talyfem/grid/grid_types/gridsimp4d.h>

#ifdef ENABLE_4D

#include <talyfem/grid/grid_types/gridbox4d.h>

#endif

#include <talyfem/input_data/input_data.h>

#include <memory>


namespace TALYFEMLIB {

/**
 * Funtion to wrap the interface of creating different meshes.
 * @param[out] pGrid where to put the pointer to the new GRID
 * @param pIdata input data options
 * @throw FileIOException if a File IO-specific error occurs.
 * @throw TALYException if a general error occurs (catch this).
 */
inline void CreateGrid(GRID*&pGrid, const InputData* pIdata) {
  const int mesh_order = basis_get_mesh_order(pIdata->basisFunction);

  if (pIdata->ifBoxGrid) {
#ifdef ENABLE_4D
    if (pIdata->nsd == 4) {
      if (pIdata->ifTriElem) {
        pGrid = new GRIDSimp4D(mesh_order);
      } else {
        pGrid = new GRIDBox4D(mesh_order);
      }
    }
#endif
    if (pIdata->nsd == 3)
      pGrid = new GRIDBox3D(mesh_order);
    if (pIdata->nsd == 2) {
      if (pIdata->ifTriElem)
        pGrid = new GRID2D(mesh_order);
      else
        pGrid = new GRIDBox2D(mesh_order);
    }
    if (pIdata->nsd == 1)
      pGrid = new GRID1D(mesh_order);

    pGrid->CreateGrid(pIdata); // Generates the grid elements and nodes
  } else {
    pGrid = new GRID();
    pGrid->ReadGrid(pIdata);
  }
}

/**
 * CreateGrid but for shared_ptr.
 * @param[out] pGrid shared_ptr to put the GRID in
 * @param pIdata input data options
 */
inline void CreateGrid(std::shared_ptr<GRID> &pGrid, const InputData* pIdata) {
  GRID* grid_tmp = NULL;
  CreateGrid(grid_tmp, pIdata);
  pGrid.reset(grid_tmp);
}

/**
 * Frees memory for the given grid
 *
 * @param pGrid pointer to grid to delete.
 */
inline void DestroyGrid(GRID*&pGrid) {
  delete pGrid;
  pGrid = NULL;
}

}  // namespace TALYFEMLIB

#endif  // GRID_GRID_TYPES_H_
