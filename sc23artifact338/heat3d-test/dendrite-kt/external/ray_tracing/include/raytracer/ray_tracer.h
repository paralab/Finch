//
// Created by boshun on 1/11/18.
//
#pragma once
#include <string>
#include <cmath>
#include <chrono>
#include "stl_reader.h"

using namespace std::chrono;
/** Define namespace for return type of ifInside
 * */
namespace InsideType {
/**
 * Type < 0   -> error
 * Type == 0  -> outside
 * Type > 0   -> inside or on surface
 */
enum Type {
  UNKNOWN = -1,
  OUTSIDE = 0,
  INSIDE = 1,
  ON_SURFACE = 2,
};
}

/**
 * return type of line intersection
 */
namespace IntersectionType {
enum Type {
  INTERSECT = 0,
  PARALLEL = 1,
  ON_EDGE = 2,
  ON_SURFACE = 3,
  NO_INTERSECTION = 4,
  ON_STARTPOINT = 5,
};
}

namespace BinningType {
enum Type {
  NOBIN = -1,
  X_1D = 0,
  Y_1D = 1,
  Z_1D = 2,
  XY_2D = 3,
  YZ_2D = 4,
  ZX_2D = 5,
};
}

class RayTracer {
 protected:
  inline static void Vec_Cross(double* dest, const double* v1, const double* v2) {
    dest[0] = v1[1] * v2[2] - v1[2] * v2[1];
    dest[1] = v1[2] * v2[0] - v1[0] * v2[2];
    dest[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }

  inline static double Vec_Dot(const double* v1, const double* v2) {
    return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
  }

  inline static void Vec_Sub(double *dest, const double *v1, const double *v2) {
    dest[0] = v1[0] - v2[0];
    dest[1] = v1[1] - v2[1];
    dest[2] = v1[2] - v2[2];
  }

  /**
   * Determine if a ray intersects with a sphere.
   * @param orig source point of ray
   * @param dir direction of ray (orig + dir = endpoint)
   * @param center center of sphere
   * @param radius radius of sphere
   * @return
   */
  static bool OutsideBoundingSphere(const double orig[3], const double dir[3],
                                    const double center[3], double radius);

  /**
   * Determine if the triangle is outside the BBOX constructed by the ray
   * @param orig
   * @param dir
   * @param center
   * @param radius
   * @return
   */
  static bool OutsideRayBBOX(const double orig[3], const double dir[3], const double center[3], double radius);

  /**
   * Determine if a line intersects with a triangle surface.
   * @param orig source point of ray
   * @param dir scaled direction of ray (orig + dir = endpoint)
   * @param vert0 vertex 0 of triangle
   * @param vert1 vertex 1 of triangle
   * @param vert2 vertex 2 of triangle
   * @return type of intersection
   */
  static IntersectionType::Type ifLineIntersectTriangle(const double orig[3],
                              const double dir[3],
                              const double vert0[3],
                              const double vert1[3],
                              const double vert2[3]);

  /**
   * Determine if a point lies within an axis-aligned bounding box.
   * @param ray_end point to test
   * @param bounds bounding box
   * @return true if outside, false if inside
   */
  static bool ifOutsideBoundingBox(const double* point, const Box3d& bounds);

 public:
#ifdef BENCHMARK
  duration<double> d_bounding_sphere;
  duration<double> d_bounding_box;
  duration<double> d_set_triangle;
  std::pair<duration<double>, unsigned long> d_in_out;
  std::pair<duration<double>, unsigned long> d_out_bounding_box;
  /// for each point tested, we have a breakdown of performance, time/no. of lines
  std::vector<std::pair<duration<double>, unsigned long>> d_in_bounding_boxs;
  std::vector<std::pair<duration<double>, unsigned long>> d_out_raybbs;
  std::vector<std::pair<duration<double>, unsigned long>> d_out_triangle_spheres;
  std::vector<std::pair<duration<double>, unsigned long>> d_line_tris;
#endif

  struct TriangleData {
    Vector3d v[3];
    Vector3d normal;
    Sphere boundingSphere() const{
      return Sphere::fromTriangle(this->v);
    };
  };
  BinningType::Type binning_type = BinningType::Type::NOBIN;
  std::vector<std::pair<std::vector<double>, std::vector<unsigned int>>> bins_1d;
  std::vector<std::vector<std::pair<std::vector<double>, std::vector<unsigned int>>>> bins_2d;
//  void setTriangles(std::vector<stl::Triangle>& lines);
//  void setTriangles(std::vector<TriangleData>&& lines);
  void setBinning();
  void binning_nobin();
  void binning_1d(int no_bins_1d);
  void binning_2d(int no_bins_dir1, int no_bins_dir2);

  RayTracer(const TriangleData * triangleData,const size_t numTriangles);

  ~RayTracer();

#ifdef BENCHMARK
  void printTimer(const std::string &printTimer = "") const;
  void cleanUpTimer();
  InsideType::Type ifInside(const double point[3], bool forceNoBin = false);
#else
  InsideType::Type ifInside(const double point[3], bool forceNoBin = false) const;
#endif
  // note: ray_end must be outside the geometry
#ifdef BENCHMARK
  InsideType::Type ifInside(const double ray_start[3], const double ray_end[3], int binno = -1);
#else
  InsideType::Type ifInside(const double ray_start[3], const double ray_end[3], int binno = -1) const;
#endif

#ifdef BENCHMARK
  InsideType::Type ifInsideBin2D(const double ray_start[3], const double ray_end[3], const unsigned int *binno);
#else
  InsideType::Type ifInsideBin2D(const double ray_start[3], const double ray_end[3], const unsigned int *binno) const;
#endif

  /**
   * Generates a random ray endpoint guaranteed to traverse bounding box (?)
   * @param pick affects which direction the ray is aimed
   * @param[out] result double[3] the result ray will be written to
   */
  void generateRay(int pick, double* result);

  void generateRay1D(const double *rayend, double* result, int &binno) const;

  void generateRay2D(const double *rayend, double* result, unsigned int binno[2]) const;

  inline const TriangleData* triangles() const {
    return triangles_;
  }

  inline const Box3d& bounds() const {
    return bounds_;
  }

  inline const size_t getNumTriangles() const{
    return numTriangles_;
  }

 protected:
  const TriangleData * triangles_;
  const size_t numTriangles_;

  Box3d bounds_;
  Vector3d ray_starts_[8];

  static Box3d calcBoundingBox(const TriangleData * triangles, size_t numTriangles);
};
