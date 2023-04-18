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
#ifndef GRID_SEGMENT_H_
#define GRID_SEGMENT_H_


namespace TALYFEMLIB {

class GRID;
class ZEROPTV;

/**
 * Segements of the grid
 */
class Segment {
 public:
  Segment();
  ~Segment();

  /**
   * get the tangent of this segment
   *
   * @param pGrid the grid pointer
   * @param e the tangent direction, which is not normalized yet [out]
   */
  void getTangent(const GRID* pGrid, ZEROPTV& e) const;

  /**
   * get the nearest point to another point on this segment
   *
   * @param pGrid the grid pointer
   * @param W the point [in]
   * @param D the point on this segment nearest to W  [out]
   * @param ratio of PD/PQ range from 0 to 1 (P is the start point of this
   *              segment, Q is the end point) [out]
   * @param e the tangent direction of PQ (already normalized) [out]
   * @return the distance between W and D
   */
  double getNearest(const GRID* pGrid, const ZEROPTV& W, ZEROPTV& D,
                    double& ratio, ZEROPTV& e) const;

  /**
   * Calculate the length of this Segment, after this call you can use member
   * varible len
   */
  void calcLen(const GRID* pGrid);

  /**
   * Get the length of this Segment, member variable len will not be updated
   */
  double getLen(const GRID* pGrid) const;

  /**
   * Get the surface id of this segment in element ONLY SUPPORT 2D
   */
  int getSurfaceID(const GRID* pGrid) const;

  /**
   * Copy from another Segment
   */
  Segment& operator=(const Segment& A);

  /**
   * Get the local position of a point
   *
   * @param epsilon local position [out]
   * @param ratio PX/PQ  P, Q are end points X is your point
   * @param pGrid pointer of grid
   */
  void getLocalPosition(ZEROPTV& epsilon, double ratio,
                        const GRID* pGrid) const;

  int nStartNode;  ///< start node
  int nEndNode;    ///< end node
  int nElm;        ///< element ID
  double previousTotalLen;  ///< previous total length
  double len;  ///< length
};

}  // namespace TALYFEMLIB

#endif  // GRID_SEGMENT_H_
