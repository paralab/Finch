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
#ifndef GRID_GRID_TYPES_GRIDBOX2D_H_
#define GRID_GRID_TYPES_GRIDBOX2D_H_

#include <string>

#include <talyfem/data_structures/zeroarray.h>
#include <talyfem/grid/grid_types/grid.h>  // parent class


namespace TALYFEMLIB {

/**
 * Simple box type grid in 2D
 *
 * This class has a number of redim functions to create the grid in different
 * manners. Due to that, it does not share much in common with GridBox1D or
 * GridBox3D in the method of grid generation.
 *
 * Note that only the basic uniform box generation has been tested.
 *
 * TODO: details of this need to be explained more fully
 * TODO: the comments for this classes functions are not very clear
 */
class GRIDBox2D : public GRID {
 public:
  /**
   * Create a GRIDBox2D object
   *
   * @param basis_function_order order of basis functions to use
   */
  explicit GRIDBox2D(int basis_function_order = 1);

  ~GRIDBox2D() { }

  /**
   * Returns the node ID of the bottom line node with x index
   *
   * @param nodeXID x index of the node start from left to right
   * @return nodeID of desired node
   */
  int GetBottomNodeID(int nodeXID) const;

  /**
   * Returns the element id for element at (nXID, nYID)
   *
   * @param nXID the x count of the element
   * @param nYID the y count of the element
   * @return the element id
   */
  int GetElmID(int nXID, int nYID) const;

  /**
   * Returns the x index of element start from left to right
   *
   * @param elm_id id of element of interest
   * @return the x index of the element
   */
  int GetElmXID(int elm_id) const;

  /**
   * Returns the node id for node at index (i, j)
   *
   * @param i the x index of the node
   * @param j the y index of the node
   * @return the node id
   */
  int GetNodeID(int i, int j) const;

  /**
   * Returns the x index of node start from left to right
   *
   * @param node_id id of node of interest
   * @return the x index of the node
   */
  int GetNodeXID(int node_id) const;

  /**
   * Returns the node ID of the top line node with x index
   *
   * @param nodeXID x index of the node, starting from left
   * @return nodeID of the desired node in the top line
   */
  int GetTopNodeID(int nodeXID) const;

  /**
   * Gets the x & y index of given node
   *
   * @param node_id id of the node of interst
   * @param[out] I the x index of this node
   * @param[out] J the x index of this node
   */
  void GetXYID(int node_id, int& I, int& J) const;

  /**
   * Gets the y coordinate of node(i,j)
   *
   * @param i start from left to right
   * @param j start from bottom to top
   * @return the y coordinate of node at (i, j)
   */
  double GetY(int i, int j) const;

  /**
   * Generates a grid above anther grid
   *
   * @param grid the grid at the bottom
   * @param maxY the top line of the grid to be generated (flat)
   * @param nGridY the number of elements on each column
   */
  void redimAbove(const GRIDBox2D& grid, double maxY, int nGridY);

  /**
   * Generates a cosine shape grid
   *
   * @param amplitude the amplitude of the cosine line
   * @param wavelength the wavelength of the cosine line
   * @param noOfWavelength the number of wavelength the bottom line contains
   * @param DY the width of this grid
   * @param nGridX the number of elements on each row
   * @param nGridY the number of elements on each column
   */
  void redimCos(double amplitude, double wavelength, int noOfWavelength,
                double DY, int nGridX, int nGridY);

  /**
   * Generates a cosine shape grid (as redimCos but with topline as the base)
   *
   * @param amplitude the amplitude of the cosine line
   * @param wavelength the wavelength of the cosine line
   * @param noOfWavelength the number of wavelength the bottom line contains
   * @param thick the width of this grid
   * @param nGridX the number of elements on each row
   * @param nGridY the number of elements on each column
   */
  void redimCosMold(double amplitude, double wavelength, int noOfWavelength,
                    double thick, int nGridX, int nGridY);

  /**
   * Generates a cosine shape grid (the same as redimCos)
   *
   * @param amplitude the amplitude of the cosine line
   * @param wavelength the wavelength of the cosine line
   * @param noOfWavelength the number of wavelength the bottom line contains
   * @param thick the width of this grid
   * @param nGridX the number of elements on each row
   * @param nGridY the number of elements on each column
   */
  void redimCosShell(double amplitude, double wavelength, int noOfWavelength,
                     double thick, int nGridX, int nGridY);

  /**
   * Generates the grid
   *
   * @param DX the length of the grid
   * @param nGridX the number of element on each row
   * @param nGridY the number of element on each column
   * @param bottomline the y position of each node on the bottom line
   * @param topline the y position of each node on the top line
   */
  void redimGeneral(double DX, int nGridX, int nGridY,
                    const double bottomline[], const double topline[]);

//  void SetBndrIndicators(InputData* pIdata);

 private:
  void CreateElementsBasis1();
  void CreateElementsBasis2();
  void CreateElementsBasis3();
  void CreateNodes(const double *dimensions);

  std::string GridTypeName() const { return std::string("GRIDBOX2D"); }

  /**
   * Generates a simple rectangle grid
   *
   * @param dimensions lengths of grid in each dimension
   * @param n_elems number of elements in each direction
   */
  void redim(const double* dimensions, const int* n_elems);

  /**
   * Creates a 2D box domain using Domain Decomposition
   *
   * @param dimensions lengths of grid in each dimension
   * @param n_elems number of elements in each direction
   */
  void redimDD(const double* dimensions, const int* n_elems);

  ZEROARRAY<double> bottomline_;  ///< y position of each node on bottom line
  ZEROARRAY<double> topline_;  ///< y position of each node on top line

  double len_x_;  ///< length of the box in the x direction
};

}  // namespace TALYFEMLIB

#endif  // GRID_GRID_TYPES_GRIDBOX2D_H_
