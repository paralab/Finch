//
// Created by boshun on 10/1/20.
//

#include <Geometry/MSH.h>
#include <talyfem/utils/utils.h>

namespace GEOMETRY {

MSH::MSH(const std::string &fileName, const InOutTest2D inOutTest)
    : inOutTest_(inOutTest) {

  int nProc = TALYFEMLIB::GetMPISize();
  int rank = TALYFEMLIB::GetMPIRank();

  std::ifstream gmsh_line2d_file(fileName.c_str(), std::ios::in);
  if (!gmsh_line2d_file) {
    throw std::runtime_error("ERROR: COULD NOT READ GMSH MESH FILE");
    assert(false);
  }

  LINE2D::Line2DData info("");
  std::string line;
  std::vector<std::vector<int>> connectivity;
  std::vector<Vector2d> line_nodes;

  while (std::getline(gmsh_line2d_file, line)) {
    std::istringstream iss(line);
    std::vector<int> int_value_temp;
    int val_temp;
    if (iss.str() == "$Nodes") {
      std::getline(gmsh_line2d_file, line);
      iss.clear();
      iss.str(line);
      // read header
      while ((iss >> val_temp)) {
        int_value_temp.push_back(val_temp);
      }
      int_value_temp.clear();
      while (std::getline(gmsh_line2d_file, line)) {
        // read end of nodes
        iss.clear();
        iss.str(line);
        if (iss.str() == "$EndNodes") {
          break;
        }
        // read nodes in each section
        while ((iss >> val_temp)) {
          int_value_temp.push_back(val_temp);
        }
        int no_nodes_per_section = int_value_temp.back();
        int_value_temp.clear();
        for (int i = 0; i < no_nodes_per_section; i++) {
          std::getline(gmsh_line2d_file, line);
        }
        for (int i = 0; i < no_nodes_per_section; i++) {
          std::getline(gmsh_line2d_file, line);
          iss.clear();
          iss.str(line);
          double val_coor;
          std::vector<double> node;
          while ((iss >> val_coor)) {
            node.push_back(val_coor);
          }
          assert(node.size() == 3 and fabs(node[2]) < 1e-15);
          line_nodes.emplace_back(node[0], node[1]);
        }
      }
    }

    if (iss.str() == "$Elements") {
      std::getline(gmsh_line2d_file, line);
      iss.clear();
      iss.str(line);
      // read header
      while ((iss >> val_temp)) {
        int_value_temp.push_back(val_temp);
      }
      int_value_temp.clear();
      while (std::getline(gmsh_line2d_file, line)) {
        // read entity
        iss.clear();
        iss.str(line);
        if (iss.str() == "$EndElements") {
          break;
        }
        while ((iss >> val_temp)) {
          int_value_temp.push_back(val_temp);
        }
        // line segments
        if (int_value_temp.size() >= 4 and int_value_temp.at(2) == 1) {
          for (int i = 0; i < int_value_temp[3]; i++) {
            std::vector<int> nodes;
            std::getline(gmsh_line2d_file, line);
            iss.clear();
            iss.str(line);
            while ((iss >> val_temp)) {
              nodes.push_back(val_temp);
            }
            connectivity.push_back(nodes);
          }
        } else {
          for (int i = 0; i < int_value_temp[3]; i++) {
            std::getline(gmsh_line2d_file, line);
          }
        }
        int_value_temp.clear();
      }
    }
  }
  info.lines.reserve(connectivity.size());
  for (const auto &r: connectivity) {
    auto v1 = line_nodes.at(r[1] - 1);
    auto v2 = line_nodes.at(r[2] - 1);
    info.lines.emplace_back(LINE2D::Line(v1, v2));
    assert(r.size() == 3);
  }
  DENDRITE_UINT numLines = info.lines.size();
  m_lines.resize(numLines);


  for (DENDRITE_UINT i = 0; i < numLines; i++) {
    for (DENDRITE_UINT d = 0; d < 2; d++) {
      m_lines[i].normal[d] = 0.0;
    }
    for (DENDRITE_UINT d = 0; d < 2; d++) {
      m_lines[i].lineCoord[d][0] = static_cast<float>(info.lines[i].v[d].data[0]);
      m_lines[i].lineCoord[d][1] = static_cast<float>(info.lines[i].v[d].data[1]);
    }
  }



  // TODO : no parallel read for 2D
  if (inOutTest_ == GEOMETRY::InOutTest2D::RAY_TRACING_2D) {
    m_raytracer = new RayTracer2D((const RayTracer2D::Line2DData *) m_lines.data(), m_lines.size());
  }
  /// Since every processor reads MSH file
  const DENDRITE_UINT perProcLine = ceil(m_lines.size()*1.0/nProc);
  m_start = rank*(perProcLine);
  m_end = (rank + 1)*(perProcLine);

  m_start = (rank*perProcLine >= m_lines.size()) ? m_lines.size() : (rank*(perProcLine));
  m_end = (m_end > m_lines.size() ? m_lines.size() : m_end);
  gmsh_line2d_file.close();
  correctNormals();

}
MSH::MSH(std::vector<std::array<DENDRITE_REAL,DIM>> & points, const int numSize, const InOutTest2D inOutTest)
:inOutTest_(inOutTest){
#ifndef DNDEBUG
  DENDRITE_REAL distX = fabs(points[0][0] - points[(numSize - 1)][0]);
  DENDRITE_REAL distY = fabs(points[0][1] - points[(numSize - 1)][1]);
  if(not(FEQUALS(distX,0.0) and FEQUALS(distY,0.0))){
    std::cout << " The surface is not a closed loop " << distX << " " << distY << "\n";
  }
#endif
  m_lines.resize(numSize - 1);
  for(int i = 0; i < numSize - 1; i++){
    m_lines[i].lineCoord[0][0] = points[i][0];
    m_lines[i].lineCoord[0][1] = points[i][1];
    m_lines[i].lineCoord[1][0] = points[i+1][0];
    m_lines[i].lineCoord[1][1] = points[i+1][1];
  }
  std::vector<std::array<DENDRITE_REAL,DIM>>  _points;
  std::swap(_points,points);
  int rank = TALYFEMLIB::GetMPIRank();
  int nProc = TALYFEMLIB::GetMPISize();

  /// Since every processor reads MSH file
  const DENDRITE_UINT perProcLine = ceil(m_lines.size()*1.0/nProc);
  m_start = rank*(perProcLine);
  m_end = (rank + 1)*(perProcLine);

  m_start = (rank*perProcLine >= m_lines.size()) ? m_lines.size() : (rank*(perProcLine));
  m_end = (m_end > m_lines.size() ? m_lines.size() : m_end);
}

void MSH::correctNormals() {
  switch (inOutTest_) {
    case RAY_TRACING_2D:
      correctNormalsRayTracing();
      break;
    case CIRCLE:
      assert((initAnalyticalSurface_));
      correctNormalsCircle();
      break;
    case RECTANGLE:
      throw std::runtime_error("Not implemented");
      break;
    default:
      throw std::runtime_error("Wrong type of inOut test");
  }
}

void MSH::initCircle(const DENDRITE_REAL radius, const DENDRITE_REAL *center) {
  assert(inOutTest_ == GEOMETRY::InOutTest2D::CIRCLE);
  assert(not(initAnalyticalSurface_));
  circle_.radius = radius;
  std::memcpy(circle_.center, center, sizeof(DENDRITE_REAL) * 2);
  initAnalyticalSurface_ = true;
}

void MSH::print() const {
  for (const auto & line : m_lines) {
    std::cout << "----------------------- Triangles -------------------------------------\n";
    std::cout << "Normal = " << line.normal[0] << " " << line.normal[1] << " \n";
    std::cout << "Coord1 = " << line.lineCoord[0][0] << " " << line.lineCoord[0][1] << " \n";
    std::cout << "Coord2 = " << line.lineCoord[1][0] << " " << line.lineCoord[1][1] << " \n";

  }
}

void MSH::printToFile(const std::string &fPrefix) const {
  char fname[256];
  sprintf(fname, "%s_%d_%d.vtu", fPrefix.c_str(), TALYFEMLIB::GetMPIRank(), TALYFEMLIB::GetMPISize());
  std::ofstream fout(fname);
  for (const auto & line : m_lines) {
    fout << "----------------------- Lines -------------------------------------\n";
    fout << "Normal = " << line.normal[0] << " " << line.normal[1] << " \n";
    fout << "Coord1 = " << line.lineCoord[0][0] << " " << line.lineCoord[0][1] << " \n";
    fout << "Coord2 = " << line.lineCoord[1][0] << " " << line.lineCoord[1][1] << " \n";
  }
  fout.close();
}

void MSH::printToCSV(const std::string &fPrefix) const {
  if (TALYFEMLIB::GetMPIRank() == 0) {
    char fname[256];
    sprintf(fname, "%s.csv", fPrefix.c_str());
    std::ofstream fout(fname);
    fout << "x0,y0,x1,y1,nx,ny\n";
    for (const auto &line : m_lines) {
      fout << line.lineCoord[0][0] << "," << line.lineCoord[0][1] << ",";
      fout << line.lineCoord[1][0] << "," << line.lineCoord[1][1] << ",";
      fout << line.normal[0] << "," << line.normal[1] << "\n";
    }
    fout.close();
  }
}

MSH::~MSH() {
  if (inOutTest_ == GEOMETRY::InOutTest2D::RAY_TRACING_2D) {
    delete m_raytracer;
  }
}

void MSH::correctNormalsCircle(){
  /// The normals are pointed towards the center of circle.
  for (auto & m_line : m_lines) {
    DENDRITE_REAL centerOfLine[DIM];
    centerOfLine[0] = 0.5*(m_line.lineCoord[0][0] + m_line.lineCoord[1][0]);
    centerOfLine[1] = 0.5*(m_line.lineCoord[0][1] + m_line.lineCoord[1][1]);

    TALYFEMLIB::ZEROPTV normal;
    normal.x() = -centerOfLine[0] + circle_.center[0];
    normal.y() = -centerOfLine[1] + circle_.center[1];
    normal.z() = 0.0;

    normal.SafeNormalize();
    for(int d = 0; d < DIM; d++){
      m_line.normal[d] = normal[d];
    }
  }

}
void MSH::correctNormalsRayTracing() {
  for (int elemID = 0; elemID < m_lines.size(); elemID++) {
    std::vector<TALYFEMLIB::ZEROPTV> pts(2);
    for (int i = 0; i < 2; i++) {
      for (int dim = 0; dim < DIM; dim++) {
        pts[i](dim) = m_lines[elemID].lineCoord[i][dim];
      }
    }
    /// calculate the normal, this normal used for calc normal_flip
    TALYFEMLIB::ZEROPTV side_1;
    for (int dim = 0; dim < DIM; dim++) {
      side_1(dim) = pts[1](dim) - pts[0](dim);
    }
    TALYFEMLIB::ZEROPTV normal;
    normal(0) = side_1(1);
    normal(1) = -side_1(0);
    normal.SafeNormalize();
    TALYFEMLIB::ZEROPTV line_centroid = (pts[0] + pts[1]) * (1.0 / 2.0);

    /// create a new point using pts[0] of the triangle and its normal, the factor is chosen with no apparent reason
    TALYFEMLIB::ZEROPTV vec = normal * ((side_1.norm()) * 0.005);
    TALYFEMLIB::ZEROPTV test_coor1 = line_centroid + vec;  // in the direction of the normal
    TALYFEMLIB::ZEROPTV test_coor2 = line_centroid - vec;  // opposite the normal

    InsideType2D::Type res = m_raytracer->ifInside(test_coor1.data());
    assert (res != InsideType2D::ON_SURFACE && res != InsideType2D::UNKNOWN);
    bool inside1 = (res == InsideType2D::INSIDE);


    // verify that at least one direction of the normal is inside the geometry
    // if this assert fails, the geometry near this point is very very thin,
    // and we can't tell which direction the normal should point by this test
    // if this happens, you can try a smaller constant for vec...
    res = m_raytracer->ifInside(test_coor2.data());
    assert (res != InsideType2D::ON_SURFACE && res != InsideType2D::UNKNOWN);
    bool inside2 = (res == InsideType2D::INSIDE);
    if (inside1 == inside2) {
      // second check with raw ray-tracer
      res = m_raytracer->ifInside(test_coor1.data(), true);
      bool inside3 = (res == InsideType2D::INSIDE);

      res = m_raytracer->ifInside(test_coor2.data(), true);
      bool inside4 = (res == InsideType2D::INSIDE);

      if (inside3 == inside4) {
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
//    std::cout << elemID << " " << gp_normal_flip << std::endl;
    for (int dim = 0; dim < DIM; dim++) {
      m_lines[elemID].normal[dim] = normal(dim) * gp_normal_flip;
    }
    if (gp_normal_flip > 0) {
      // change order of nodes from 0, 1 to 1, 0
      DENDRITE_REAL tempnode[2];
      for (int dim = 0; dim < DIM; dim++) {
        tempnode[dim] = m_lines[elemID].lineCoord[0][dim];
        m_lines[elemID].lineCoord[0][dim] = m_lines[elemID].lineCoord[1][dim];
        m_lines[elemID].lineCoord[1][dim] = tempnode[dim];
      }
    }
  }
}

bool MSH::ifInside(const DENDRITE_REAL *position) const {
  switch (inOutTest_) {
    case GEOMETRY::InOutTest2D::RAY_TRACING_2D: {
      return (m_raytracer->ifInside(position) == InsideType2D::Type::INSIDE);
      break;
    }
    case GEOMETRY::InOutTest2D::CIRCLE: {
      assert(initAnalyticalSurface_);
      auto &center = circle_.center;
      auto &radius = circle_.radius;

      DENDRITE_REAL dist[2];
      dist[0] = center[0] - position[0];
      dist[1] = center[1] - position[1];

      DENDRITE_REAL distFromCenter = dist[0] * dist[0] + dist[1] * dist[1];
      if (distFromCenter < radius * radius) {
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
