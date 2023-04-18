//
// Created by maksbh on 8/12/20.
//

#include <Geometry/STL.h>
#include <talyfem/utils/utils.h>

namespace GEOMETRY {

STL::STL(const std::string &fileName, const InOutTest inOutTest)
    : inOutTest_(inOutTest) {

  int nProc = TALYFEMLIB::GetMPISize();
  int rank = TALYFEMLIB::GetMPIRank();
  // Modified from: https://github.com/dillonhuff/stl_parser/blob/master/parse_stl.cpp
  std::ifstream stlFile(fileName.c_str(), std::ios::in | std::ios::binary);
  if (!stlFile) {
    throw std::runtime_error("ERROR: COULD NOT READ STL FILE\n");
  }

  char header_info[80] = "";
  char n_triangles[4];
  stlFile.read(header_info, 80);
  stlFile.read(n_triangles, 4);

  DENDRITE_UINT *r = (DENDRITE_UINT *) n_triangles;
  DENDRITE_UINT numTriangles = *r;
  if (numTriangles > 1000000) {
    throw std::runtime_error("ERROR: STL can't be this large (>1000000), something is wrong!");
  }
  DENDRITE_UINT numTrianglesToRead = ceil(numTriangles*1.0/(nProc*1.0));
  DENDRITE_UINT numTrianglesStart = (rank*numTrianglesToRead >= numTriangles) ? numTriangles : (rank*numTrianglesToRead);
  DENDRITE_UINT numTrianglesEnd = ((numTrianglesStart + numTrianglesToRead) > numTriangles) ? numTriangles : (numTrianglesStart + numTrianglesToRead);

  numTrianglesToRead = numTrianglesEnd - numTrianglesStart;
  m_start = 0;
  m_end = numTrianglesToRead;
  m_triangles.resize(numTrianglesToRead);

  static constexpr DENDRITE_UINT eachTriangleSize = 50; // normals = 12, 3 coords = 12*3, dummy = 2
  stlFile.seekg(eachTriangleSize*numTrianglesStart, std::ios::cur);

  for (DENDRITE_UINT i = 0; i < numTrianglesToRead; i++) {
    for (DENDRITE_UINT d = 0; d < 3; d++) {
      // this normal is override in flipNormal(), but do not comment out this line as it causes reading STL to fail.
      m_triangles[i].normal[d] = static_cast<DENDRITE_REAL>(parse_float(stlFile));
    }

    for (DENDRITE_UINT d = 0; d < 3; d++) {
      m_triangles[i].triangleCoord[d][0] = static_cast<DENDRITE_REAL>(parse_float(stlFile));
      m_triangles[i].triangleCoord[d][1] = static_cast<DENDRITE_REAL>(parse_float(stlFile));
      m_triangles[i].triangleCoord[d][2] = static_cast<DENDRITE_REAL>(parse_float(stlFile));
    }
    char dummy[2];
    stlFile.read(dummy, 2);
  }
  stlFile.close();

  // TODO : For parallel ray tracing we need DA information.
  if (inOutTest_==GEOMETRY::InOutTest::RAY_TRACING) {
    m_start = numTrianglesStart;
    m_end = numTrianglesEnd;
    MPI_Comm comm = MPI_COMM_WORLD;
    std::vector<int> eachProcData(nProc);
    MPI_Allgather(&numTrianglesToRead, 1, MPI_INT, eachProcData.data(), 1, MPI_INT, comm);
#ifndef NDEBUG
    int totalNumTriangles = std::accumulate(eachProcData.begin(), eachProcData.end(), 0);
    assert(totalNumTriangles==numTriangles);
#endif
    std::vector<int> displacement(nProc, 0);
    for (int i = 1; i < displacement.size(); i++) {
      displacement[i] = displacement[i - 1] + eachProcData[i - 1];
    }

    std::vector<Triangles> totalTriangles(numTriangles);
    MPI_Allgatherv(m_triangles.data(), m_triangles.size(), GEOMETRY::Triangles::dataType(),
                   totalTriangles.data(), eachProcData.data(), displacement.data(), GEOMETRY::Triangles::dataType(), comm);
    MPI_Barrier(comm);
    std::swap(m_triangles, totalTriangles);
    m_raytracer = new RayTracer((const RayTracer::TriangleData *) m_triangles.data(), m_triangles.size());
  }
}


void STL::initSphere(DENDRITE_REAL radius, DENDRITE_REAL *center) {
  assert(inOutTest_==GEOMETRY::InOutTest::SPHERE);
  assert(not(initAnalyticalSurface_));
  sphere_.radius = radius;
  std::memcpy(sphere_.center, center, sizeof(DENDRITE_REAL)*3);
  initAnalyticalSurface_ = true;
}

void STL::correctNormals() {
  for (int elemID = 0; elemID < m_triangles.size(); elemID++) {
    std::vector<TALYFEMLIB::ZEROPTV> pts(3);
    for (int i = 0; i < 3; i++) {
      for (int dim = 0; dim < DIM; dim++) {
        pts[i](dim) = m_triangles[elemID].triangleCoord[i][dim];
      }
    }
    /// calculate the normal, this normal used for calc normal_flip
    TALYFEMLIB::ZEROPTV side_1, side_2;
    for (int dim = 0; dim < DIM; dim++) {
      side_1(dim) = pts[1](dim) - pts[0](dim);
      side_2(dim) = pts[2](dim) - pts[0](dim);
    }
    TALYFEMLIB::ZEROPTV normal;
    normal.crossProduct(side_1, side_2);
    normal.SafeNormalize();
    TALYFEMLIB::ZEROPTV triangle_centroid = (pts[0] + pts[1] + pts[2])*(1.0/3.0);

    /// create a new point using pts[0] of the triangle and its normal, the factor is chosen with no apparent reason
    TALYFEMLIB::ZEROPTV vec = normal*((side_1.norm() + side_2.norm())*0.005);
    TALYFEMLIB::ZEROPTV test_coor1 = triangle_centroid + vec;  // in the direction of the normal
    TALYFEMLIB::ZEROPTV test_coor2 = triangle_centroid - vec;  // opposite the normal

    InsideType::Type res = m_raytracer->ifInside(test_coor1.data());
    assert (res!=InsideType::ON_SURFACE && res!=InsideType::UNKNOWN);
    bool inside1 = (res==InsideType::INSIDE);


    // verify that at least one direction of the normal is inside the geometry
    // if this assert fails, the geometry near this point is very very thin,
    // and we can't tell which direction the normal should point by this test
    // if this happens, you can try a smaller constant for vec...
    res = m_raytracer->ifInside(test_coor2.data());
    assert (res!=InsideType::ON_SURFACE && res!=InsideType::UNKNOWN);
    bool inside2 = (res==InsideType::INSIDE);
    if (inside1==inside2) {
      // second check with raw ray-tracer
      res = m_raytracer->ifInside(test_coor1.data(), true);
      bool inside3 = (res==InsideType::INSIDE);

      res = m_raytracer->ifInside(test_coor2.data(), true);
      bool inside4 = (res==InsideType::INSIDE);

      if (inside3==inside4) {
        PrintStatus("even with fall back:", test_coor1, test_coor2);
//        throw TALYFEMLIB::TALYException() << "Surface normal flip confused for ??";
      } else {
        inside1 = inside3;
        inside2 = inside4;
      }
    }

    /// normals are expected to point inside the geometry, so flip if inside1 is false
    /// gp_normal_flip is multiplied by normal, so 1 means "don't flip"
    int gp_normal_flip = (inside1 ? 1 : -1);

    for (int dim = 0; dim < DIM; dim++) {
      m_triangles[elemID].normal[dim] = normal(dim) * gp_normal_flip;
    }
    if (gp_normal_flip > 0) {
      // change order of nodes from 0, 1, 2 to 2, 1, 0
      DENDRITE_REAL tempnode[3];
      for (int dim = 0; dim < DIM; dim++) {
        tempnode[dim] = m_triangles[elemID].triangleCoord[0][dim];
        m_triangles[elemID].triangleCoord[0][dim] = m_triangles[elemID].triangleCoord[2][dim];
        m_triangles[elemID].triangleCoord[2][dim] = tempnode[dim];
      }
    }
  }
}

void STL::print() const {
  for (const auto &m_triangle : m_triangles) {
    std::cout << "----------------------- Triangles -------------------------------------\n";
    std::cout << "Normal = " << m_triangle.normal[0] << " " << m_triangle.normal[1] << " " << m_triangle.normal[2] << " \n";
    std::cout << "Coord1 = " << m_triangle.triangleCoord[0][0] << " " << m_triangle.triangleCoord[0][1] << " "
              << m_triangle.triangleCoord[0][2] << " \n";
    std::cout << "Coord2 = " << m_triangle.triangleCoord[1][0] << " " << m_triangle.triangleCoord[1][1] << " "
              << m_triangle.triangleCoord[1][2] << " \n";
    std::cout << "Coord3 = " << m_triangle.triangleCoord[2][0] << " " << m_triangle.triangleCoord[2][1] << " "
              << m_triangle.triangleCoord[2][2] << " \n";

  }
}

void STL::printToFile(const std::string &fPrefix) const {
  char fname[256];
  sprintf(fname, "%s_%d_%d.vtu", fPrefix.c_str(), TALYFEMLIB::GetMPIRank(), TALYFEMLIB::GetMPISize());
  std::ofstream fout(fname);
  for (const auto &m_triangle : m_triangles) {
    fout << "----------------------- Triangles -------------------------------------\n";
    fout << "Normal = " << m_triangle.normal[0] << " " << m_triangle.normal[1] << " " << m_triangle.normal[2] << " \n";
    fout << "Coord1 = " << m_triangle.triangleCoord[0][0] << " " << m_triangle.triangleCoord[0][1] << " "
         << m_triangle.triangleCoord[0][2] << " \n";
    fout << "Coord2 = " << m_triangle.triangleCoord[1][0] << " " << m_triangle.triangleCoord[1][1] << " "
         << m_triangle.triangleCoord[1][2] << " \n";
    fout << "Coord3 = " << m_triangle.triangleCoord[2][0] << " " << m_triangle.triangleCoord[2][1] << " "
         << m_triangle.triangleCoord[2][2] << " \n";

  }
  fout.close();

}

STL::~STL() {
  if (inOutTest_==GEOMETRY::InOutTest::RAY_TRACING) {
    delete m_raytracer;
  }
}

bool STL::ifInside(const DENDRITE_REAL *position) const {
  switch (inOutTest_) {
    case GEOMETRY::InOutTest::RAY_TRACING: {
      return (m_raytracer->ifInside(position)==InsideType::Type::INSIDE or m_raytracer->ifInside(position)==InsideType::Type::ON_SURFACE);
      break;
    }
    case GEOMETRY::InOutTest::SPHERE: {
      assert(initAnalyticalSurface_);
      auto &center = sphere_.center;
      auto &radius = sphere_.radius;

      DENDRITE_REAL dist[3];
      dist[0] = center[0] - position[0];
      dist[1] = center[1] - position[1];
      dist[2] = center[2] - position[2];

      DENDRITE_REAL distFromCenter = dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];
      if (distFromCenter < radius*radius) {
        return true;
      } else {
        return false;
      }
      break;
    }
    default:throw std::logic_error("In Out test can not be of this type");
      break;
  }
}

}