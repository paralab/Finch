//
// Created by boshun on 10/1/20.
//

#ifndef DENDRITEKT_MSH_H
#define DENDRITEKT_MSH_H

#include <DataTypes.h>
#include <raytracer/ray_tracer_2d.h>
namespace GEOMETRY {
struct Lines {
  DENDRITE_REAL lineCoord[2][2];
  DENDRITE_REAL normal[2];

  static MPI_Datatype dataType() {
    static bool first = true;
    static MPI_Datatype _datatype;
    if (first) {
      first = false;
      MPI_Type_contiguous(sizeof(Lines), MPI_BYTE, &_datatype);
      MPI_Type_commit(&_datatype);
    }
    return _datatype;
  }
};

enum InOutTest2D : DENDRITE_UINT {
  RAY_TRACING_2D = 0,
  CIRCLE = 1,
  RECTANGLE = 2,
};

class MSH {
 private:

  struct {
    DENDRITE_REAL center[2];
    DENDRITE_REAL radius;
  } circle_;

  bool initAnalyticalSurface_ = false;

  const InOutTest2D inOutTest_;
 protected:
  std::vector<Lines> m_lines;
  RayTracer2D *m_raytracer = NULL;

  DENDRITE_UINT m_start = 0;
  DENDRITE_UINT m_end = 0;

  void correctNormalsRayTracing();
  void correctNormalsCircle();
 public:
  /**
   * @brief: MSH constructor
   *
   * @param fileName MSH filename
   * @param inOutTest Type of IN - OUT test
   */
  MSH(const std::string &fileName, InOutTest2D inOutTest = InOutTest2D::RAY_TRACING_2D);

  /**
   *
   * @param points : Points in 2D Must be of size (numSize X 2). Must be an ordered set.
   * The first point must match with last point to get a closed loop.
   * After the loop the points will be deleted to prevent memory
   * @param numSize : number of points
   * @param inOutTest Type of IN-OUT test
   */
  MSH(std::vector<std::array<DENDRITE_REAL,DIM>> & points, const int numSize, InOutTest2D inOutTest = InOutTest2D::RAY_TRACING_2D);

  void initCircle(const DENDRITE_REAL radius, const DENDRITE_REAL *center);

  inline const std::vector<Lines> &getLines() const {
    return m_lines;
  }

  void updateGeoLocationTranslate(const TALYFEMLIB::ZEROPTV &LinearDisplacementUpdate) {
    for (int line_id = 0; line_id < m_lines.size(); line_id++) {
      for (int no_node = 0; no_node < 2; no_node++) {
        for (int dim = 0; dim < 2; dim++) {
          m_lines[line_id].lineCoord[no_node][dim] += LinearDisplacementUpdate(dim);
        }
      }
    }
    updateRayTracer();
  }

  void updateGeoLocationRotation(const std::vector<TALYFEMLIB::ZEROPTV> &RotMat,
                                 const TALYFEMLIB::ZEROPTV &center,
                                 const TALYFEMLIB::ZEROPTV &init_translate) {
    for (int line_id = 0; line_id < m_lines.size(); line_id++) {
      for (int no_node = 0; no_node < 2; no_node++) {
        TALYFEMLIB::ZEROPTV newLocation;
        for (int dim = 0; dim < 3; dim++) {
          for (int k = 0; k < 2; k++) {
            newLocation(dim) += (m_lines[line_id].lineCoord[no_node][k] + init_translate[k] - center(k)) * RotMat[dim](k);
          }
        }
        for (int dim = 0; dim < 2; dim++) {
          m_lines[line_id].lineCoord[no_node][dim] = center[dim] + newLocation[dim] - init_translate[dim];
        }
      }
    }
    updateRayTracer();
  }

  void updateRayTracer() {
    delete m_raytracer;
    m_raytracer = new RayTracer2D((const RayTracer2D::Line2DData *) m_lines.data(), m_lines.size());
  }

  inline const InOutTest2D getInOutTest() const {
    return inOutTest_;
  }

  inline const RayTracer2D *getRayTracer() const {
    return m_raytracer;
  };
  void print() const;

  void printToFile(const std::string &fPrefix) const;

  void printToCSV(const std::string &fPrefix) const;

  ~MSH();

  bool ifInside(const DENDRITE_REAL *position) const;

  inline const DENDRITE_UINT getStart() const {
    return m_start;
  }

  inline const DENDRITE_UINT getEnd() const {
    return m_end;
  }

  void correctNormals();

};

}

#endif //NSHTIBM_KT_DENDRITE_KT_INCLUDE_GEOMETRY_MSH_H_
