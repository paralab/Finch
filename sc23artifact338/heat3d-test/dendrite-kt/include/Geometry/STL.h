//
// Created by maksbh on 8/12/20.
//

#ifndef DENDRITEKT_STL_H
#define DENDRITEKT_STL_H

#include <DataTypes.h>
#include <raytracer/ray_tracer.h>
namespace GEOMETRY{
struct Triangles{
    DENDRITE_REAL triangleCoord[3][3];
    DENDRITE_REAL normal[3];

    static MPI_Datatype dataType(){
      static bool         first = true;
      static MPI_Datatype _datatype;
      if (first){
        first = false;
        MPI_Type_contiguous(sizeof(Triangles), MPI_BYTE, &_datatype);
        MPI_Type_commit(&_datatype);
      }
      return _datatype;
  }
};

enum InOutTest:DENDRITE_UINT {
  RAY_TRACING = 0,
  SPHERE = 1,
};


class STL{
 private:

  struct{
    DENDRITE_REAL center[3];
    DENDRITE_REAL radius;
  }sphere_;

  bool initAnalyticalSurface_ = false;

  static inline float parse_float(std::ifstream& s){
    char f_buf[sizeof(float)];
    s.read(f_buf, 4);
    float* fptr = (float*) f_buf;
    return *fptr;
  }

  const InOutTest inOutTest_;
 protected:
  std::vector<Triangles> m_triangles;
  RayTracer *m_raytracer = NULL;

  DENDRITE_UINT m_start = 0;
  DENDRITE_UINT m_end = 0;

 public:
  /**
   * @brief: STL constructor
   *
   * @param fileName STL filename
   * @param inOutTest Type of IN - OUT test
   */
  STL(const std::string & fileName, InOutTest inOutTest = InOutTest::RAY_TRACING);

  void initSphere(DENDRITE_REAL radius, DENDRITE_REAL * center);

  inline const std::vector<Triangles> & getTriangles() const {
    return m_triangles;
  }

  void updateGeoLocationTranslate(const TALYFEMLIB::ZEROPTV &LinearDisplacementUpdate) {
    for (int tri_id = 0; tri_id < m_triangles.size(); tri_id++) {
      for (int no_node = 0; no_node < 3; no_node++) {
        for (int dim = 0; dim < 3; dim++) {
          m_triangles[tri_id].triangleCoord[no_node][dim] += LinearDisplacementUpdate(dim);
        }
      }
    }
    updateRayTracer();
  }

  void updateGeoLocationRotation(const std::vector<TALYFEMLIB::ZEROPTV> &RotMat,
                                 const TALYFEMLIB::ZEROPTV &center,
                                 const TALYFEMLIB::ZEROPTV &init_translate) {
    for (int tri_id = 0; tri_id < m_triangles.size(); tri_id++) {
      for (int no_node = 0; no_node < 3; no_node++) {
        TALYFEMLIB::ZEROPTV newLocation;
        for (int dim = 0; dim < 3; dim++) {
          for (int k = 0; k < 3; k++) {
            newLocation(dim) += (m_triangles[tri_id].triangleCoord[no_node][k] + init_translate[k] - center(k))*RotMat[dim](k);
          }
        }
        for (int dim = 0; dim < 3; dim++) {
          m_triangles[tri_id].triangleCoord[no_node][dim] = newLocation[dim] + center[dim] - init_translate[dim];
        }
      }
    }
    updateRayTracer();
  }

  void updateRayTracer() {
    delete m_raytracer;
    m_raytracer = new RayTracer((const RayTracer::TriangleData *) m_triangles.data(), m_triangles.size());
  }

  inline const InOutTest getInOutTest() const {
    return inOutTest_;
  }

  inline const RayTracer * getRayTracer() const{
    return m_raytracer;
  };
  void print() const;

  void printToFile(const std::string& fPrefix) const ;

  ~STL();

  bool ifInside(const DENDRITE_REAL * position) const;

  inline DENDRITE_UINT getStart() const{
    return m_start;
  }

  inline DENDRITE_UINT getEnd() const{
      return m_end;
  }

  void correctNormals();
};

}
#endif //DENDRITEKT_STL_H
