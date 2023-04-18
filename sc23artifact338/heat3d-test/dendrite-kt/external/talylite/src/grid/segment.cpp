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
#include <talyfem/grid/segment.h>

#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/elem.h>
#include <talyfem/grid/zeroptv.h>
#include <talyfem/grid/node.h>


namespace TALYFEMLIB {

Segment::Segment() { }

Segment::~Segment() { }

void Segment::getTangent(const GRID* pGrid, ZEROPTV& e) const {
  const NODE* pNode1 = pGrid->GetNode(nStartNode);
  const NODE* pNode2 = pGrid->GetNode(nEndNode);
  e.x() = pNode2->x() - pNode1->x();
  e.y() = pNode2->y() - pNode1->y();
  e.z() = pNode2->z() - pNode1->z();
}

double Segment::getNearest(const GRID* pGrid, const ZEROPTV& W,
                           ZEROPTV& D, double& ratio,
                           ZEROPTV& e) const {
  const NODE* pNode1 = pGrid->GetNode(nStartNode);
  const NODE* pNode2 = pGrid->GetNode(nEndNode);
  ZEROPTV P, Q;
  P.x() = pNode1->x();
  P.y() = pNode1->y();
  Q.x() = pNode2->x();
  Q.y() = pNode2->y();

  ZEROPTV PQ, PW;
  PQ.x() = Q.x() - P.x();
  PQ.y() = Q.y() - P.y();
  PW.x() = W.x() - P.x();
  PW.y() = W.y() - P.y();
  double distancePQ = PQ.Normalize();
  double distancePD = PW.innerProduct(PQ);

  if (distancePD < 0) {
    D = P;
  } else if (distancePD > distancePQ) {
    D = Q;
  } else {
    D.x() = P.x() + distancePD * PQ.x();
    D.y() = P.y() + distancePD * PQ.y();
  }

  ratio = distancePD / distancePQ;
  e = PQ;
  return W.distanceTo(D);
}

void Segment::calcLen(const GRID* pGrid) {
  len = getLen(pGrid);
}

double Segment::getLen(const GRID* pGrid) const {
  const NODE* pNode1 = pGrid->GetNode(nStartNode);
  const NODE* pNode2 = pGrid->GetNode(nEndNode);
  ZEROPTV P, Q;
  P.x() = pNode1->x();
  P.y() = pNode1->y();
  Q.x() = pNode2->x();
  Q.y() = pNode2->y();
  return P.distanceTo(Q);
}

int Segment::getSurfaceID(const GRID* pGrid) const {
  const ELEM* pElm = pGrid->GetElm(nElm);
  switch (pElm->elmType()) {
    case kElem2dBox:
      // {-1, 1, 4},  start 4 ->surfaceID 1
      // {+1, 2, 3},  start 2 ->surfaceID 2
      // {-2, 1, 2},  start 1 ->surfaceID 3
      // {+2, 3, 4}   start 3 ->surfaceID 4
    {
      if (nStartNode == pElm->ElemToLocalNodeID(3)) {
        return 1;
      }
      if (nStartNode == pElm->ElemToLocalNodeID(1)) {
        return 2;
      }
      if (nStartNode == pElm->ElemToLocalNodeID(0)) {
        return 3;
      }
      if (nStartNode == pElm->ElemToLocalNodeID(2)) {
        return 4;
      }
    }
      break;
    case kElem2dTriangle:
      // {-1, 1, 3},  start 3 -> surfaceID 1
      // {+1, 2, 3},  start 2 -> surfaceID 2
      // {-2, 1, 2}   start 1 -> surfaceID 3
    {
      if (nStartNode == pElm->ElemToLocalNodeID(2)) {
        return 1;
      }
      if (nStartNode == pElm->ElemToLocalNodeID(1)) {
        return 2;
      }
      if (nStartNode == pElm->ElemToLocalNodeID(0)) {
        return 3;
      }
    }
      break;
    default:
      throw TALYException() << "Segment::getSurfaceID does not support "
          << "elmType " << pElm->elmType();
  }
  return 1;
}

Segment& Segment::operator=(const Segment& A) {
  nStartNode = A.nStartNode;
  nEndNode = A.nEndNode;
  nElm = A.nElm;
  previousTotalLen = A.previousTotalLen;
  len = A.len;

  return *this;
}

void Segment::getLocalPosition(ZEROPTV& epsilon, double ratio,
                               const GRID* pGrid) const {
  const ELEM* pElm = pGrid->GetElm(nElm);

  switch (pElm->elmType()) {
    case kElem2dBox:
      if (nStartNode == pElm->ElemToLocalNodeID(0)) {
        epsilon.x() = -1 + ratio * 2;
        epsilon.y() = -1;
      } else if (nStartNode == pElm->ElemToLocalNodeID(1)) {
        epsilon.x() = +1;
        epsilon.y() = -1 + ratio * 2;
      } else if (nStartNode == pElm->ElemToLocalNodeID(2)) {
        epsilon.x() = +1 - ratio * 2;
        epsilon.y() = +1;
      } else if (nStartNode == pElm->ElemToLocalNodeID(3)) {
        epsilon.x() = -1;
        epsilon.y() = +1 - ratio * 2;
      }
      epsilon.z() = 0;
      break;
    case kElem2dTriangle:
      if (nStartNode == pElm->ElemToLocalNodeID(0)) {
        epsilon.x() = -1 + ratio * 2;
        epsilon.y() = -1;
      } else if (nStartNode == pElm->ElemToLocalNodeID(1)) {
        epsilon.x() = +1;
        epsilon.y() = -1 + ratio * 2;
      } else if (nStartNode == pElm->ElemToLocalNodeID(2)) {
        epsilon.x() = -1;
        epsilon.y() = +1 - ratio * 2;
      }
      epsilon.z() = 0;
      break;
    default:
      printf("Segment::getLocalPosition only support ElmT3n2D & ElmB4n2D\n");
      exit(0);
  }
}

}  // namespace TALYFEMLIB
