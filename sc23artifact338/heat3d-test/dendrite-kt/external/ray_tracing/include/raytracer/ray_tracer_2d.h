//
// Created by boshun on 7/13/20.
//
#pragma once
#include <string>
#include <cmath>
#include <chrono>
#include "gmsh_2d_reader.h"

using namespace std::chrono;
/** Define namespace for return type of ifInside
 * */
namespace InsideType2D {
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
namespace IntersectionType2D {
enum Type {
  INTERSECT = 0,
  PARALLEL = 1,
  ON_NODE = 2,
  ON_LINE = 3,
  NO_INTERSECTION = 4,
  ON_STARTPOINT = 5,
};
}

namespace BinningType2D {
enum Type {
  NOBIN = -1,
  X_1D = 0,
  Y_1D = 1,
};
}

class RayTracer2D {
 protected:
  inline static double Vec_Dot(const double *v1, const double *v2) {
    return (v1[0] * v2[0] + v1[1] * v2[1]);
  }

  inline static void Vec_Sub(double *dest, const double *v1, const double *v2) {
    dest[0] = v1[0] - v2[0];
    dest[1] = v1[1] - v2[1];
  }

  /**
   * Determine if a ray intersects with a sphere.
   * @param orig source point of ray
   * @param dir direction of ray (orig + dir = endpoint)
   * @param center center of sphere
   * @param radius radius of sphere
   * @return
   */
  static bool OutsideBoundingCircle(const double orig[2], const double dir[2],
                                    const double center[2], double radius);

  /**
   * Determine if the line is outside the BBOX constructed by the ray
   * @param orig
   * @param dir
   * @param center
   * @param radius
   * @return
   */
  static bool OutsideRayBBOX(const double orig[2], const double dir[2], const double center[2], double radius);

  /**
   * Determine if a line intersects with a line.
   * @param orig source point of ray
   * @param dir scaled direction of ray (orig + dir = endpoint)
   * @param vert0 vertex 0 of line
   * @param vert1 vertex 1 of line
   * @return type of intersection
   */
  static IntersectionType2D::Type ifLineIntersectLine(const double orig[2],
                                                      const double dest[2],
                                                      const double vert0[2],
                                                      const double vert1[2]);

  /**
   * Determine if a point lies within an axis-aligned bounding box.
   * @param ray_end point to test
   * @param bounds bounding box
   * @return true if outside, false if inside
   */
  static bool ifOutsideBoundingBox(const double *point, const Box2d &bounds);

 public:
#ifdef BENCHMARK
  duration<double> d_bounding_circle;
  duration<double> d_bounding_box;
  duration<double> d_set_line;
  std::pair<duration<double>, unsigned long> d_in_out;
  std::pair<duration<double>, unsigned long> d_out_bounding_box;
  /// for each point tested, we have a breakdown of performance, time/no. of lines
  std::vector<std::pair<duration<double>, unsigned long>> d_in_bounding_boxs;
  std::vector<std::pair<duration<double>, unsigned long>> d_out_raybbs;
  std::vector<std::pair<duration<double>, unsigned long>> d_out_line_circles;
  std::vector<std::pair<duration<double>, unsigned long>> d_line_lines;
#endif

  struct Line2DData {
    Vector2d v[2];
    Vector2d normal;
    Circle boundingCircle() const{
      return Circle::fromLine(this->v);
    };
  };
  BinningType2D::Type binning_type = BinningType2D::Type::NOBIN;
  std::vector<std::pair<std::vector<double>, std::vector<unsigned int>>> bins_1d;
//  void setLines(std::vector<LINE2D::Line> &lines);
//  void setLines(std::vector<Line2DData> &&lines);
  void setBinning();
  void binning_nobin();
  void binning_1d(int no_bins_1d);

  RayTracer2D(const Line2DData * lineData,const size_t numLines);

  ~RayTracer2D();
#ifdef BENCHMARK
  void printTimer(const std::string &printTimer = "") const;
  void cleanUpTimer();
  InsideType2D::Type ifInside(const double point[2], bool forceNoBin = false);
#else
  InsideType2D::Type ifInside(const double point[2], bool forceNoBin = false) const;
#endif
  // note: ray_end must be outside the geometry
#ifdef BENCHMARK
  InsideType2D::Type ifInside(const double ray_start[2], const double ray_end[2], int binno = -1);
#else
  InsideType2D::Type ifInside(const double ray_start[2], const double ray_end[2], int binno = -1) const;
#endif

  /**
   * Generates a random ray endpoint guaranteed to traverse bounding box (?)
   * @param pick affects which direction the ray is aimed
   * @param[out] result double[3] the result ray will be written to
   */
  void generateRay(int pick, double *result);

  void generateRay1D(const double *rayend, double *result, int &binno) const;

  inline const Line2DData* lines() const {
    return lines_;
  }

  inline const Box2d &bounds() const {
    return bounds_;
  }

  inline const size_t getNumLines() const{
    return numLines_;
  }

 protected:
  const Line2DData * lines_;
  const size_t numLines_;

  Box2d bounds_;
  Vector2d ray_starts_[4];

  static Box2d calcBoundingBox(const Line2DData *lines, size_t numLines);
};
