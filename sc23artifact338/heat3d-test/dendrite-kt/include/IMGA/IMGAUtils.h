////
//// Created by maksbh on 7/2/20.
////
#pragma once

#include <Geometry/Geometry.h>
#include <DataTypes.h>

struct GaussPoint {
  DENDRITE_REAL pos[DIM];
  DENDRITE_REAL normal[DIM];
  DENDRITE_UINT geomID;
  DENDRITE_UINT elemID; // Only Active when the In-Out Test is not Analytical
  DENDRITE_REAL elemArea; // Area of triangle (in 3D)  or line length (in 2D) after splitting

  bool isInsideDomain(const DomainExtents & domainExtents){
    for(int d = 0; d < DIM; d++){
      if((pos[d] <= domainExtents.physicalDADomain.min[d]) or (pos[d] >= domainExtents.physicalDADomain.max[d])){
        // Its safe to have equality sign especially at the max limit to ignore the effect of +ve SubDA boundaries
        return false;
      }
    }
    return true;
  }
  static MPI_Datatype dataType() {
    static bool first = true;
    static MPI_Datatype _datatype;
    if (first) {
      first = false;
      MPI_Type_contiguous(sizeof(GaussPoint), MPI_BYTE, &_datatype);
      MPI_Type_commit(&_datatype);
    }
    return _datatype;
  }
};

struct GeomRefinement {
  DENDRITE_UINT octantLevel = 0; /// Octant refinement level near the geometry
  DENDRITE_UINT maxSplitIteration = 1000; /// maximum Number of times the split needs to be carried out
  DENDRITE_REAL ratioArea = 0.25; /// Ratio of triangle required/ Octant Area ratio
};


#if(DIM == 3)
/**
 * @brief computes the triangle with Heron's formula
 * @param triCoords Coordinates of triangle (3 X 3)
 * @return Area of triangle
 */
static DENDRITE_REAL computeTriangleArea(const DENDRITE_REAL triCoords[][3]) {
  /// calculate the normal and size of the surface
  TALYFEMLIB::ZEROPTV pts_0{triCoords[0][0], triCoords[0][1], triCoords[0][2]};
  TALYFEMLIB::ZEROPTV pts_1{triCoords[1][0], triCoords[1][1], triCoords[1][2]};
  TALYFEMLIB::ZEROPTV pts_2{triCoords[2][0], triCoords[2][1], triCoords[2][2]};
  TALYFEMLIB::ZEROPTV side_1, side_2;
  for (int dim = 0; dim < DIM; dim++) {
    side_1(dim) = pts_1(dim) - pts_0(dim);
    side_2(dim) = pts_2(dim) - pts_0(dim);
  }

  /// this normal is either outside or inside, need to be pointed inside afterwards
  TALYFEMLIB::ZEROPTV normal;
  normal.crossProduct(side_1, side_2);
  const DENDRITE_REAL area = 0.5 * normal.SafeNormalize();  /// note: this line also normalizes normal!
  return area;
}
/**
 * @brief Splits each triangles into 4. Used to fill more Gauss points.
 * @param [out] splitTriangles the splitted triangles
 * @param [in] coords original coordinates of  triangle to be splitted
 * @param numSplit number of times the split has to be carried out
 *
 */

static void performTriangleSplitting(std::vector<GEOMETRY::Triangles> &splitTriangles,
                                     const DENDRITE_REAL coords[][3],
                                     const DENDRITE_UINT numSplit) {
  int currentSplit = 0;
  const DENDRITE_UINT numTriangles = (1u << (2 * numSplit));
  splitTriangles.resize(numTriangles);
  std::vector<GEOMETRY::Triangles> _splitTriangles(numTriangles);
  std::memcpy(splitTriangles[0].triangleCoord, coords, sizeof(DENDRITE_REAL) * 9);
  while (currentSplit < numSplit) {
    const DENDRITE_UINT triangleSplitLevel = (1u << (2 * currentSplit)); // Triangle at this split level
    for (DENDRITE_UINT i = 0; i < triangleSplitLevel; i++) {
      DENDRITE_REAL midEdges[3][DIM]; // 3 mid - point edges of triangles.

      midEdges[0][0] = 0.5 * (splitTriangles[i].triangleCoord[0][0] + splitTriangles[i].triangleCoord[1][0]);
      midEdges[0][1] = 0.5 * (splitTriangles[i].triangleCoord[0][1] + splitTriangles[i].triangleCoord[1][1]);
      midEdges[0][2] = 0.5 * (splitTriangles[i].triangleCoord[0][2] + splitTriangles[i].triangleCoord[1][2]);

      midEdges[1][0] = 0.5 * (splitTriangles[i].triangleCoord[1][0] + splitTriangles[i].triangleCoord[2][0]);
      midEdges[1][1] = 0.5 * (splitTriangles[i].triangleCoord[1][1] + splitTriangles[i].triangleCoord[2][1]);
      midEdges[1][2] = 0.5 * (splitTriangles[i].triangleCoord[1][2] + splitTriangles[i].triangleCoord[2][2]);

      midEdges[2][0] = 0.5 * (splitTriangles[i].triangleCoord[2][0] + splitTriangles[i].triangleCoord[0][0]);
      midEdges[2][1] = 0.5 * (splitTriangles[i].triangleCoord[2][1] + splitTriangles[i].triangleCoord[0][1]);
      midEdges[2][2] = 0.5 * (splitTriangles[i].triangleCoord[2][2] + splitTriangles[i].triangleCoord[0][2]);


      // Triangle 1 (a, m0, m2)
      std::memcpy(&_splitTriangles[4 * i + 0].triangleCoord[0][0], &splitTriangles[i].triangleCoord[0][0], sizeof(DENDRITE_REAL) * DIM);
      std::memcpy(&_splitTriangles[4 * i + 0].triangleCoord[1][0], &midEdges[0][0], sizeof(DENDRITE_REAL) * DIM);
      std::memcpy(&_splitTriangles[4 * i + 0].triangleCoord[2][0], &midEdges[2][0], sizeof(DENDRITE_REAL) * DIM);

      // Triangle 2 (m0, b, m1)
      std::memcpy(&_splitTriangles[4 * i + 1].triangleCoord[0][0], &midEdges[0][0], sizeof(DENDRITE_REAL) * DIM);
      std::memcpy(&_splitTriangles[4 * i + 1].triangleCoord[1][0], &splitTriangles[i].triangleCoord[1][0], sizeof(DENDRITE_REAL) * DIM);
      std::memcpy(&_splitTriangles[4 * i + 1].triangleCoord[2][0], &midEdges[1][0], sizeof(DENDRITE_REAL) * DIM);

      // Triangle 3 (m1, c, m2)
      std::memcpy(&_splitTriangles[4 * i + 2].triangleCoord[0][0], &midEdges[1][0], sizeof(DENDRITE_REAL) * DIM);
      std::memcpy(&_splitTriangles[4 * i + 2].triangleCoord[1][0], &splitTriangles[i].triangleCoord[2][0], sizeof(DENDRITE_REAL) * DIM);
      std::memcpy(&_splitTriangles[4 * i + 2].triangleCoord[2][0], &midEdges[2][0], sizeof(DENDRITE_REAL) * DIM);

      // Triangle 4 (m0,m1,m2)
      std::memcpy(&_splitTriangles[4 * i + 3].triangleCoord[0][0], &midEdges[0][0], sizeof(DENDRITE_REAL) * DIM);
      std::memcpy(&_splitTriangles[4 * i + 3].triangleCoord[1][0], &midEdges[1][0], sizeof(DENDRITE_REAL) * DIM);
      std::memcpy(&_splitTriangles[4 * i + 3].triangleCoord[2][0], &midEdges[2][0], sizeof(DENDRITE_REAL) * DIM);

    }
    std::swap(splitTriangles, _splitTriangles);
    currentSplit++;
  }
}

#endif

#if(DIM == 2)
/**
 * @brief computes the line element length
 * @param lineCoords Coordinates of line element (2 X 2)
 * @return length of line element
 */
static DENDRITE_REAL computeLineLength(const DENDRITE_REAL lineCoords[][2]) {
  static_assert(DIM == 2, "computeLineLength should only be called with DIM == 2");
  /// calculate the normal and size of the surface
  TALYFEMLIB::ZEROPTV pts_0{lineCoords[0][0], lineCoords[0][1], 0.0};
  TALYFEMLIB::ZEROPTV pts_1{lineCoords[1][0], lineCoords[1][1], 0.0};
  TALYFEMLIB::ZEROPTV line = pts_1 - pts_0;

  return line.norm();

}

static void performLineSplitting(std::vector<GEOMETRY::Lines> &splitLines,
                                     const DENDRITE_REAL coords[][2],
                                     const DENDRITE_UINT numSplit) {
  int currentSplit = 0;
  const DENDRITE_UINT numLines = (1u << (numSplit));
  splitLines.resize(numLines);
  std::vector<GEOMETRY::Lines> _splitLines(numLines);
  std::memcpy(splitLines[0].lineCoord, coords, sizeof(DENDRITE_REAL) * 4);
  while (currentSplit < numSplit) {
    const DENDRITE_UINT lineSplitLevel = (1u << (currentSplit)); // Triangle at this split level
    for (DENDRITE_UINT i = 0; i < lineSplitLevel; i++) {
      DENDRITE_REAL midNode[1][DIM]; // 1 mid - point of line.

      midNode[0][0] = 0.5 * (splitLines[i].lineCoord[0][0] + splitLines[i].lineCoord[1][0]);
      midNode[0][1] = 0.5 * (splitLines[i].lineCoord[0][1] + splitLines[i].lineCoord[1][1]);


      // line1 1 (a, m0)
      std::memcpy(&_splitLines[2 * i + 0].lineCoord[0][0], &splitLines[i].lineCoord[0][0], sizeof(DENDRITE_REAL) * DIM);
      std::memcpy(&_splitLines[2 * i + 0].lineCoord[1][0], &midNode[0][0], sizeof(DENDRITE_REAL) * DIM);

      // line2 2 (m0, b)
      std::memcpy(&_splitLines[2 * i + 1].lineCoord[0][0], &midNode[0][0], sizeof(DENDRITE_REAL) * DIM);
      std::memcpy(&_splitLines[2 * i + 1].lineCoord[1][0], &splitLines[i].lineCoord[1][0], sizeof(DENDRITE_REAL) * DIM);

    }
    std::swap(splitLines, _splitLines);
    currentSplit++;
  }

}

static void gpPositionLines(const DENDRITE_REAL lineCoord[][2], DENDRITE_REAL *_pos, int gp_no){
  DENDRITE_REAL x1 = lineCoord[0][0];
  DENDRITE_REAL y1 = lineCoord[0][1];
  DENDRITE_REAL x2 = lineCoord[1][0];
  DENDRITE_REAL y2 = lineCoord[1][1];
  if (gp_no == 0) {
    _pos[0] = x1 + (x2 - x1) / 2 * (1 - sqrt(3) / 3);
    _pos[1] = y1 + (y2 - y1) / 2 * (1 - sqrt(3) / 3);
  } else {
    _pos[0] = x1 + (x2 - x1) / 2 * (1 + sqrt(3) / 3);
    _pos[1] = y1 + (y2 - y1) / 2 * (1 + sqrt(3) / 3);
  }
}
#endif


/**
 *
 * @param domainExtents domainExtents
 * @param levelOctant level of refinement near the object
 * @return the minimum of the area
 */
static DENDRITE_REAL computeOctantArea(const DomainExtents &domainExtents, const DENDRITE_UINT &levelOctant) {
  DENDRITE_REAL sideLength[DIM];
  for (int d = 0; d < DIM; d++) {
    sideLength[d] = domainExtents.fullDADomain.max[d] - domainExtents.fullDADomain.min[d];
  }

#if (DIM == 2)
  DENDRITE_REAL minLength = *std::min_element(sideLength, sideLength+DIM);
  return minLength/((1u << levelOctant)*1.0);
#endif
#if (DIM == 3)
  DENDRITE_REAL area[DIM];
  for (int d = 0; d < DIM; d++) {
    area[d] = sideLength[d] * sideLength[(d + 1) % DIM];
  }
  DENDRITE_REAL minArea = *std::min_element(area, area + DIM);
  return minArea / ((1u << levelOctant * 2) * 1.0);
#endif

}

static void TESTGaussPoints(const std::string &filename, const std::vector<GaussPoint> &gaussPoints_test, const std::string& mode = "a") {
  FILE *fp = fopen(filename.c_str(), mode.c_str());
  assert(fp != nullptr);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  fprintf(fp, "gp_x,gp_y,gp_nx,gp_ny,geomID,lineID,lineLen,rank\n");
  for (const auto &gp : gaussPoints_test) {
#if (DIM == 2)
    fprintf(fp, "%.10e,%.10e,"
                "%.10e,%.10e,"
                "%1d,%1d,%10e,%1d\n",
            gp.pos[0], gp.pos[1],
            gp.normal[0], gp.normal[1],
            gp.geomID, gp.elemID, gp.elemArea, rank);
#endif
#if (DIM == 3)
    fprintf(fp, "%.10e,%.10e,%.10e,"
                "%.10e,%.10e,%.10e,"
                "%1d,%1d,%10e,%1d\n",
            gp.pos[0], gp.pos[1], gp.pos[2],
            gp.normal[0], gp.normal[1], gp.normal[2],
            gp.geomID, gp.elemID, gp.elemArea, rank);
#endif
  }
  fclose(fp);
}


static void printGaussPointsToFile(const std::string &filename, const std::vector<NodeAndValues<DENDRITE_REAL>> &gaussPoints) {
  FILE *fp = fopen(filename.c_str(),"w");
  assert(fp != nullptr);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  fprintf(fp, "gp_x gp_y gp_nx gp_ny geomID lineID lineLen rank\n");
  for (const auto &gp : gaussPoints) {
#if (DIM == 2)
    fprintf(fp, "%.10e %.10e "
                "%.10e %.10e "
                "%1d %1d %10e %1d\n",
            gp.location[0], gp.location[1],
            gp.normal[0], gp.normal[1],
            gp.geomID, gp.elemID, gp.elemArea, rank);
#endif
#if (DIM == 3)
    fprintf(fp, "%.10e,%.10e,%.10e,"
                "%.10e,%.10e,%.10e,"
                "%1d,%1d,%10e,%1d\n",
            gp.location[0], gp.location[1], gp.location[2],
            gp.normal[0], gp.normal[1], gp.normal[2],
            gp.geomID, gp.elemID, gp.elemArea, rank);
#endif
  }
  fclose(fp);
}

