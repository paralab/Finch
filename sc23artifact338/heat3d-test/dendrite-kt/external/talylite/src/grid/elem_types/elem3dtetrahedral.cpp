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
#include <talyfem/grid/elem_types/elem3dtetrahedral.h>
#include <talyfem/grid/node.h>
#include <talyfem/grid/grid_types/grid.h>


namespace TALYFEMLIB {

const int* ELEM3dTetrahedral::GetSurfaceCheckArray() const {
  /*
   * (from http://gmsh.info/doc/texinfo/gmsh.html)
   *                         Y
   *                       .
   *                     ,/
   *                    /
   *                 2
   *               ,/|`\
   *             ,/  |  `\
   *           ,/    '.   `\
   *         ,/       |     `\
   *       ,/         |       `\
   *      0-----------'.--------1 --> X
   *       `\.         |      ,/
   *          `\.      |    ,/
   *             `\.   '. ,/
   *                `\. |/
   *                   `3
   *                      `\.
   *                         ` Z
   */

  static const int T4n3DCheckArray[] = {
      // {surface_id, node_id1, node_id2, node_id3 }
      1, 1, 2, 3,    // that one weird face (hypotenuse?)
      2, 0, 3, 2,    // -x (YZ plane)
      3, 0, 1, 3,    // -y (XZ plane)
      4, 1, 0, 2     // -z (XY plane)
  };

  /*
        Tetrahedron10:

                2
              ,/|`\
            ,/  |  `\
          ,6    '.   `5
        ,/       8     `\
      ,/         |       `\
     0--------4--'.--------1
      `\.         |      ,/
         `\.      |    ,9
            `7.   '. ,/
               `\. |/
                  `3
   */
  static const int T10n3DCheckArray[] = {
      // {surface_id, node_id1, node_id2, node_id3, ... }
      1, 1, 2, 3, 5, 8, 9,  // hypotenuse
      2, 0, 3, 2, 7, 8, 6,  // YZ plane (normal = X- direction)
      3, 0, 1, 3, 4, 9, 7,  // XZ plane (normal = Y- direction)
      4, 1, 0, 2, 4, 6, 5   // XY plane (normal = Z- direction)
  };

  /*
    Cubic tetrahedron node order from the Gmsh source code:
    1. Vertices (isoparametric coords): (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)
    2. Edges (filled from first local node ID to second): {0, 1}, {1, 2}, {2, 0}, {3, 0}, {3, 2}, {3, 1}
    3. Faces (at center): {0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {3, 1, 2}
   */
  static const int T20n3DCheckArray[(10 + 1) * 4] = {
      // {surface_id, node_id1, node_id2, node_id3, ... }
      1, 1, 2, 3, 6, 7, 13, 12, 14, 15, 19,  // hypotenuse
      2, 0, 3, 2, 11, 10, 12, 13, 8, 9, 18,  // X-
      3, 0, 1, 3, 4, 5, 15, 14, 10, 11, 17,  // Y-
      4, 1, 0, 2, 5, 4, 9, 8, 7, 6, 16       // Z-
  };

  switch (n_nodes()) {
    case 4: return T4n3DCheckArray;
    case 10: return T10n3DCheckArray;
    case 20: return T20n3DCheckArray;
    default: throw NotImplementedException();
  }
}

int ELEM3dTetrahedral::GetNodesPerSurface() const {
  switch (n_nodes()) {
    case 4: return 3;
    case 10: return 6;
    case 20: return 10;
    default: throw NotImplementedException();
  }
}

ZEROPTV ELEM3dTetrahedral::CalculateNormal(const GRID* grid, int srf_id) const {
  const int* nodes = GetNodesInSurface(srf_id);
  const ZEROPTV& p1 = grid->GetNode(ElemToLocalNodeID(nodes[0]))->location();
  const ZEROPTV& p2 = grid->GetNode(ElemToLocalNodeID(nodes[1]))->location();
  const ZEROPTV& p3 = grid->GetNode(ElemToLocalNodeID(nodes[2]))->location();
  ZEROPTV normal;
  normal.crossProduct(p2 - p1, p3 - p1);
  normal.Normalize();
  return normal;
}

// calculate the volume using the Cayley-Menger determinant
double ELEM3dTetrahedral::GetMeasure(const GRID* p_grid) const {
  // get the nodes of the element
  const NODE* node0 = p_grid->node_array_[node_id_array(0)];
  const NODE* node1 = p_grid->node_array_[node_id_array(1)];
  const NODE* node2 = p_grid->node_array_[node_id_array(2)];
  const NODE* node3 = p_grid->node_array_[node_id_array(3)];
  // calculate lengths of all edges lenij is the distance between i and j
  // variables a-j are the squares of the lengths as needed by the determinant
  // difference between two node locations is a ZEROPTV value on which the
  // norm() function is called for the length
  const double len01 = (node0->location() - node1->location()).norm();
  const double a = len01 * len01;
  const double len02 = (node0->location() - node2->location()).norm();
  const double b = len02 * len02;
  const double len03 = (node0->location() - node3->location()).norm();
  const double c = len03 * len03;
  const double len12 = (node1->location() - node2->location()).norm();
  const double d = len12 * len12;
  const double len13 = (node1->location() - node3->location()).norm();
  const double e = len13 * len13;
  const double len23 = (node2->location() - node3->location()).norm();
  const double f = len23 * len23;
  const double volume_squared = ((a+b+c+d+e+f) * (a*f + b*e + c*d) -
      2.0*a*f*(a+f) - 2*b*e*(b+e) - 2.0*c*d*(c+d) -
      0.5*(a+f)*(b+e)*(c+d) +
      0.5*(a-f)*(b-e)*(c-d)) / 144.0;
  return sqrt(volume_squared);
}

kBasisFunction ELEM3dTetrahedral::basis_function() const {
  switch (n_nodes()) {
    case 4: return BASIS_LINEAR;
    case 10: return BASIS_QUADRATIC;
    case 20: return BASIS_CUBIC;
    default: throw NotImplementedException();
  }
}

double ELEM3dTetrahedral::Volume(const GRID * grid) const {
  auto v1 = get_node_loc(grid, 1) - get_node_loc(grid, 0);
  auto v2 = get_node_loc(grid, 2) - get_node_loc(grid, 0);
  auto v3 = get_node_loc(grid, 3) - get_node_loc(grid, 0);

  ZEROPTV v1xv2;
  v1xv2.crossProduct(v1, v2);

  // Scalar triple product
  double stp = v3.innerProduct(v1xv2);

  return 1.0 / 6.0 * fabs(stp);
}

double ELEM3dTetrahedral::Angle(const GRID * grid, QualityMetricType type) const {
  static const int num_angles = 12;

  double min = 3.14159;
  double max = 0;
  double sum = 0;

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = j + 1; k < 4; k++) {
        if ((i == j) || (i == k) || (j == k)) {
          continue;
        }

        auto v1 = get_node_loc(grid, j) - get_node_loc(grid, i);
        auto v2 = get_node_loc(grid, k) - get_node_loc(grid, i);
        double angle = fabs(v1.angleTo(v2));

        if (angle < min) {
          min = angle;
        }

        if (angle > max) {
          max = angle;
        }

        sum += angle;
      }
    }
  }

  switch (type) {
    case kMin:
      return min;
    case kMax:
      return max;
    default:
      break;
  }

  return sum / (double)num_angles;
}

double ELEM3dTetrahedral::FaceArea(const GRID * grid, QualityMetricType type) const {
  throw NotImplementedException();
}

double ELEM3dTetrahedral::AspectRatio(const GRID * grid) const {

  // double ideal_circum_r =
  //   sqrt(6.0) / 4.0 * cbrt(
  //     6.0 * sqrt(2.0) / 2.0 * volume);

  // double ideal_circum_vol =
  //     4.0 / 3.0 * 3.14159 * ideal_circum_r * ideal_circum_r * ideal_circum_r;

  throw NotImplementedException();
}

}  // namespace TALYFEMLIB
