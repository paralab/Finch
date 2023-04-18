//
// Created by lofquist on 1/18/18.
//

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include "raytracer/ray_tracer.h"

bool RayTracer::OutsideBoundingSphere(const double *orig, const double *dir, const double *center, const double radius) {
  double w[3], b, Pb[3], distance_square;
  Vec_Sub(w, center, orig);
  b = Vec_Dot(w, dir) / Vec_Dot(dir, dir);
  Pb[0] = orig[0] + b * dir[0];
  Pb[1] = orig[1] + b * dir[1];
  Pb[2] = orig[2] + b * dir[2];
  distance_square = (Pb[0] - center[0]) * (Pb[0] - center[0]) +
      (Pb[1] - center[1]) * (Pb[1] - center[1]) +
      (Pb[2] - center[2]) * (Pb[2] - center[2]);
  return (distance_square > radius * radius);
}

bool RayTracer::OutsideRayBBOX(const double *orig, const double *dir, const double *center, const double radius) {

  double lb[3], ub[3];
  lb[0] = std::min(orig[0], orig[0] + dir[0]) - radius;
  lb[1] = std::min(orig[1], orig[1] + dir[1]) - radius;
  lb[2] = std::min(orig[2], orig[2] + dir[2]) - radius;
  ub[0] = std::max(orig[0], orig[0] + dir[0]) + radius;
  ub[1] = std::max(orig[1], orig[1] + dir[1]) + radius;
  ub[2] = std::max(orig[2], orig[2] + dir[2]) + radius;
  return center[0] < lb[0] || center[1] < lb[1] || center[2] < lb[2]
      || center[0] > ub[0] || center[1] > ub[1] || center[2] > ub[2];
}

IntersectionType::Type
RayTracer::ifLineIntersectTriangle(const double *orig, const double *dir, const double *vert0, const double *vert1,
                                   const double *vert2) {
  double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3], ray_length_scale = 1.0;
  Vec_Sub(edge1, vert1, vert0);
  Vec_Sub(edge2, vert2, vert0);
  Vec_Cross(pvec, dir, edge2);
  double det = Vec_Dot(edge1, pvec);
  const double epsilon_edge = 1e-7;
  const double epsilon_surface = 1e-7;
  if (det < 1e-10 && det > -1e-10) {
    /// Ray parallel to the triangle
    return IntersectionType::PARALLEL;
  } else {
    const double inv_det = 1.0 / det;
    Vec_Sub(tvec, orig, vert0);
    const double u = Vec_Dot(tvec, pvec) * inv_det;
    Vec_Cross(qvec, tvec, edge1);
    const double v = Vec_Dot(dir, qvec) * inv_det;
    const double t = Vec_Dot(edge2, qvec) * inv_det;
    /// Judge if the intersection point is inside the triangle
    /// Outside the triangle && outside reach
    if (u < -epsilon_edge || v < -epsilon_edge ||
        u > 1 + epsilon_edge || v > 1 + epsilon_edge || (u + v) > 1 + epsilon_edge || t > 1 + epsilon_surface / ray_length_scale || t < 0) {
      return IntersectionType::NO_INTERSECTION;
    } else {
      /// On the triangle edge (possible outcome are on edge or on surface)
      if (u < epsilon_edge || v < epsilon_edge || (u + v) > 1 - epsilon_edge) {
        /// Judge t (t already satisfies condition)
        if (t > 1 - epsilon_surface / ray_length_scale) {
          return IntersectionType::ON_SURFACE;
        } else {
          return IntersectionType::ON_EDGE;
        }
      } else {
        /// Inside triangle (possible outcome are on surface or intersect)
        if (t > 1 - epsilon_surface / ray_length_scale) {
          return IntersectionType::ON_SURFACE;
        } else {
          return IntersectionType::INTERSECT;
        }
      }
    }
    /*    /// Judge if the intersection Vector3d is inside the Triangle, works slow
        if (u < -epsilon_edge || v < -epsilon_edge ||
            u > 1 + epsilon_edge || v > 1 + epsilon_edge || (u + v) > 1 + epsilon_edge || t > 1 + epsilon_edge / ray_length) {
          return IntersectionType::NO_INTERSECTION;
        } else {
          if (u < 1 - epsilon_edge && v < 1 - epsilon_edge && u > epsilon_edge && v > epsilon_edge && u + v < 1 - epsilon_edge) {
            /// Judge t
            if (t < 1 - epsilon_edge / ray_length && t > epsilon_edge / ray_length) {
              return IntersectionType::INTERSECT;
            } else {
              if (t < epsilon_edge / ray_length || t > 1 - epsilon_edge / ray_length) {
                return IntersectionType::ON_SURFACE;
              }
              return IntersectionType::NO_INTERSECTION;
            }
          } else {
            return IntersectionType::ON_EDGE;
          }
        }*/
  }
}

bool RayTracer::ifOutsideBoundingBox(const double *ray_end, const Box3d &bounds) {
  return (ray_end[0] < bounds.min[0]) || (ray_end[0] > bounds.max[0])
      || (ray_end[1] < bounds.min[1]) || (ray_end[1] > bounds.max[1])
      || (ray_end[2] < bounds.min[2]) || (ray_end[2] > bounds.max[2]);
}

#ifdef BENCHMARK
InsideType::Type RayTracer::ifInside(const double *ray_start, const double *ray_end, const int binno) {
#else
InsideType::Type RayTracer::ifInside(const double *ray_start, const double *ray_end, const int binno) const {
#endif
#ifdef BENCHMARK
  high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
#endif
  /// Ignore points that are outside the bounding box for STL file
  if (ifOutsideBoundingBox(ray_end, bounds_)) {
#ifdef BENCHMARK
    high_resolution_clock::time_point t_outbb = std::chrono::high_resolution_clock::now();
    d_out_bounding_box.first += duration_cast<duration<double>>(t_outbb - t0);
    d_out_bounding_box.second++;
#endif
    return InsideType::OUTSIDE;
  }

  int noOfIntersections = 0;
  double dir[3];
  Vec_Sub(dir, ray_end, ray_start);
#ifdef BENCHMARK
  std::pair<duration<double>, unsigned long> d_in_bounding_box;
  std::pair<duration<double>, unsigned long> d_out_raybb;
  std::pair<duration<double>, unsigned long> d_out_triangle_sphere;
  std::pair<duration<double>, unsigned long> d_line_tri;
  t0 = std::chrono::high_resolution_clock::now();
  duration<double> duration_line_tri = std::chrono::milliseconds::zero();
#endif
  std::vector<unsigned int> triIDS;
  if (binno < 0) {
    triIDS.clear();
    triIDS.resize(numTriangles_);
    std::iota(std::begin(triIDS), std::end(triIDS), 0);
  } else {
    triIDS = bins_1d[binno].second;
  }
  for (const auto &triID: triIDS) {
    const auto &tri = triangles_[triID];

#ifdef BENCHMARK
    high_resolution_clock::time_point t_out_s = std::chrono::high_resolution_clock::now();
#endif


    ///  Determine if the ray intersects with the bounding sphere or a triangle, skipped if not
    if (OutsideRayBBOX(ray_start, dir, tri.boundingSphere().center, tri.boundingSphere().radius)) {
#ifdef BENCHMARK
      d_in_bounding_box.second++;
      high_resolution_clock::time_point t_out_raybb_e = std::chrono::high_resolution_clock::now();
      d_out_raybb.second++;
      d_out_raybb.first += duration_cast<duration<double>>(t_out_raybb_e - t_out_s);
#endif
      continue;
    }
#ifdef BENCHMARK
    t_out_s = std::chrono::high_resolution_clock::now();
#endif
    ///  Determine if the ray intersects with the bounding sphere or a triangle, skipped if not
    if (OutsideBoundingSphere(ray_start, dir, tri.boundingSphere().center, tri.boundingSphere().radius)) {
#ifdef BENCHMARK
      d_in_bounding_box.second++;
      high_resolution_clock::time_point t_out_sphere_e = std::chrono::high_resolution_clock::now();
      d_out_triangle_sphere.second++;
      d_out_triangle_sphere.first += duration_cast<duration<double>>(t_out_sphere_e - t_out_s);
#endif
      continue;
    }

#ifdef BENCHMARK
    d_in_bounding_box.second++;
    d_line_tri.second++;
#endif
    high_resolution_clock::time_point t_s_line_tri = std::chrono::high_resolution_clock::now();
    IntersectionType::Type result_from_ifIntersect =
        ifLineIntersectTriangle(ray_start,
                                dir,
                                tri.v[0].data,
                                tri.v[1].data,
                                tri.v[2].data);
    /// If ON_EDGE or ON_SURFACE or ON_STARTPOINT happens, then this position is vague, jump out of the loop and choose another ray
    if (result_from_ifIntersect == IntersectionType::ON_EDGE ||
        result_from_ifIntersect == IntersectionType::ON_STARTPOINT) {
#ifdef BENCHMARK
      high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
      d_in_bounding_box.first += duration_cast<duration<double>>(t1 - t0);
      d_out_raybbs.emplace_back(d_out_raybb);
      d_out_triangle_spheres.emplace_back(d_out_triangle_sphere);
      d_in_bounding_boxs.emplace_back(d_in_bounding_box);
      high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
      d_line_tri.first += duration_cast<duration<double>>(t2 - t_s_line_tri);
      d_line_tris.emplace_back(d_line_tri);
#endif
      return InsideType::UNKNOWN;
    }
    if (result_from_ifIntersect == IntersectionType::ON_SURFACE) {
#ifdef BENCHMARK
      high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
      d_in_bounding_box.first += duration_cast<duration<double>>(t1 - t0);
      d_out_raybbs.emplace_back(d_out_raybb);
      d_out_triangle_spheres.emplace_back(d_out_triangle_sphere);
      d_in_bounding_boxs.emplace_back(d_in_bounding_box);
      high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
      d_line_tri.first += duration_cast<duration<double>>(t2 - t_s_line_tri);
      d_line_tris.emplace_back(d_line_tri);
#endif
      return InsideType::ON_SURFACE;
    }
    if (result_from_ifIntersect == IntersectionType::INTERSECT) {
      noOfIntersections++;
    }
#ifdef BENCHMARK
    high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    duration_line_tri += duration_cast<duration<double>>(t2 - t_s_line_tri);
#endif
  }
#ifdef BENCHMARK
  high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  d_line_tri.first += duration_line_tri;
  d_in_bounding_box.first += duration_cast<duration<double>>(t1 - t0);
  d_out_raybbs.emplace_back(d_out_raybb);
  d_out_triangle_spheres.emplace_back(d_out_triangle_sphere);
  d_in_bounding_boxs.emplace_back(d_in_bounding_box);
  d_line_tris.emplace_back(d_line_tri);
#endif
  if ((noOfIntersections % 2) != 0) {
    return InsideType::INSIDE;
  } else {
    return InsideType::OUTSIDE;
  }
}

#ifdef BENCHMARK
InsideType::Type RayTracer::ifInsideBin2D(const double *ray_start, const double *ray_end, const unsigned int *binno) {
#else
InsideType::Type RayTracer::ifInsideBin2D(const double *ray_start, const double *ray_end, const unsigned int *binno) const {
#endif
#ifdef BENCHMARK
  high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
#endif
  /// Ignore points that are outside the bounding box for STL file
  if (ifOutsideBoundingBox(ray_end, bounds_)) {
#ifdef BENCHMARK
    high_resolution_clock::time_point t_outbb = std::chrono::high_resolution_clock::now();
    d_out_bounding_box.first += duration_cast<duration<double>>(t_outbb - t0);
    d_out_bounding_box.second++;
#endif
    return InsideType::OUTSIDE;
  }

  int noOfIntersections = 0;
  double dir[3];
  Vec_Sub(dir, ray_end, ray_start);
#ifdef BENCHMARK
  std::pair<duration<double>, unsigned long> d_in_bounding_box;
  std::pair<duration<double>, unsigned long> d_out_raybb;
  std::pair<duration<double>, unsigned long> d_out_triangle_sphere;
  std::pair<duration<double>, unsigned long> d_line_tri;
  t0 = std::chrono::high_resolution_clock::now();
  duration<double> duration_line_tri = std::chrono::milliseconds::zero();
#endif

  for (const auto &triID: bins_2d[binno[0]][binno[1]].second) {
    const auto &tri = triangles_[triID];

#ifdef BENCHMARK
    high_resolution_clock::time_point t_out_s = std::chrono::high_resolution_clock::now();
#endif


    ///  Determine if the ray intersects with the bounding sphere or a triangle, skipped if not
    if (OutsideRayBBOX(ray_start, dir, tri.boundingSphere().center, tri.boundingSphere().radius)) {
#ifdef BENCHMARK
      d_in_bounding_box.second++;
      high_resolution_clock::time_point t_out_raybb_e = std::chrono::high_resolution_clock::now();
      d_out_raybb.second++;
      d_out_raybb.first += duration_cast<duration<double>>(t_out_raybb_e - t_out_s);
#endif
      continue;
    }
#ifdef BENCHMARK
    t_out_s = std::chrono::high_resolution_clock::now();
#endif
    ///  Determine if the ray intersects with the bounding sphere or a triangle, skipped if not
    if (OutsideBoundingSphere(ray_start, dir, tri.boundingSphere().center, tri.boundingSphere().radius)) {
#ifdef BENCHMARK
      d_in_bounding_box.second++;
      high_resolution_clock::time_point t_out_sphere_e = std::chrono::high_resolution_clock::now();
      d_out_triangle_sphere.second++;
      d_out_triangle_sphere.first += duration_cast<duration<double>>(t_out_sphere_e - t_out_s);
#endif
      continue;
    }

#ifdef BENCHMARK
    d_in_bounding_box.second++;
    d_line_tri.second++;
#endif
    high_resolution_clock::time_point t_s_line_tri = std::chrono::high_resolution_clock::now();
    IntersectionType::Type result_from_ifIntersect =
        ifLineIntersectTriangle(ray_start,
                                dir,
                                tri.v[0].data,
                                tri.v[1].data,
                                tri.v[2].data);
    /// If ON_EDGE or ON_SURFACE or ON_STARTPOINT happens, then this position is vague, jump out of the loop and choose another ray
    if (result_from_ifIntersect == IntersectionType::ON_EDGE ||
        result_from_ifIntersect == IntersectionType::ON_STARTPOINT) {
#ifdef BENCHMARK
      high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
      d_in_bounding_box.first += duration_cast<duration<double>>(t1 - t0);
      d_out_raybbs.emplace_back(d_out_raybb);
      d_out_triangle_spheres.emplace_back(d_out_triangle_sphere);
      d_in_bounding_boxs.emplace_back(d_in_bounding_box);
      high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
      d_line_tri.first += duration_cast<duration<double>>(t2 - t_s_line_tri);
      d_line_tris.emplace_back(d_line_tri);
#endif
      return InsideType::UNKNOWN;
    }
    if (result_from_ifIntersect == IntersectionType::ON_SURFACE) {
#ifdef BENCHMARK
      high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
      d_in_bounding_box.first += duration_cast<duration<double>>(t1 - t0);
      d_out_raybbs.emplace_back(d_out_raybb);
      d_out_triangle_spheres.emplace_back(d_out_triangle_sphere);
      d_in_bounding_boxs.emplace_back(d_in_bounding_box);
      high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
      d_line_tri.first += duration_cast<duration<double>>(t2 - t_s_line_tri);
      d_line_tris.emplace_back(d_line_tri);
#endif
      return InsideType::ON_SURFACE;
    }
    if (result_from_ifIntersect == IntersectionType::INTERSECT) {
      noOfIntersections++;
    }
#ifdef BENCHMARK
    high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    duration_line_tri += duration_cast<duration<double>>(t2 - t_s_line_tri);
#endif
  }
#ifdef BENCHMARK
  high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  d_line_tri.first += duration_line_tri;
  d_in_bounding_box.first += duration_cast<duration<double>>(t1 - t0);
  d_out_raybbs.emplace_back(d_out_raybb);
  d_out_triangle_spheres.emplace_back(d_out_triangle_sphere);
  d_in_bounding_boxs.emplace_back(d_in_bounding_box);
  d_line_tris.emplace_back(d_line_tri);
#endif
  if ((noOfIntersections % 2) != 0) {
    return InsideType::INSIDE;
  } else {
    return InsideType::OUTSIDE;
  }
}

void RayTracer::generateRay(int pick, double *result) {
  const auto &boundingBox = bounds_;

  double random_perturb[3] = {
      1 + double(rand()) / RAND_MAX,
      1 + double(rand()) / RAND_MAX,
      1 + double(rand()) / RAND_MAX
  };
  double x_length = boundingBox.max[0] - boundingBox.min[0];
  double y_length = boundingBox.max[1] - boundingBox.min[1];
  double z_length = boundingBox.max[2] - boundingBox.min[2];

  switch (pick) {
    case 0:result[0] = (boundingBox.min[0] - x_length * random_perturb[0]);
      result[1] = (boundingBox.min[1] - y_length * random_perturb[1]);
      result[2] = (boundingBox.min[2] - z_length * random_perturb[2]);
      break;
    case 1:result[0] = (boundingBox.max[0] + x_length * random_perturb[0]);
      result[1] = (boundingBox.min[1] - y_length * random_perturb[1]);
      result[2] = (boundingBox.min[2] - z_length * random_perturb[2]);
      break;
    case 2:result[0] = (boundingBox.min[0] - x_length * random_perturb[0]);
      result[1] = (boundingBox.max[1] + y_length * random_perturb[1]);
      result[2] = (boundingBox.min[2] - z_length * random_perturb[2]);
      break;
    case 3:result[0] = (boundingBox.max[0] + x_length * random_perturb[0]);
      result[1] = (boundingBox.max[1] + y_length * random_perturb[1]);
      result[2] = (boundingBox.min[2] - z_length * random_perturb[2]);
      break;
    case 4:result[0] = (boundingBox.min[0] - x_length * random_perturb[0]);
      result[1] = (boundingBox.min[1] - y_length * random_perturb[1]);
      result[2] = (boundingBox.max[2] + z_length * random_perturb[2]);
      break;
    case 5:result[0] = (boundingBox.max[0] + x_length * random_perturb[0]);
      result[1] = (boundingBox.min[1] - y_length * random_perturb[1]);
      result[2] = (boundingBox.max[2] + z_length * random_perturb[2]);
      break;
    case 6:result[0] = (boundingBox.min[0] - x_length * random_perturb[0]);
      result[1] = (boundingBox.max[1] + y_length * random_perturb[1]);
      result[2] = (boundingBox.max[2] + z_length * random_perturb[2]);
      break;
    case 7:
    default:result[0] = (boundingBox.max[0] + x_length * random_perturb[0]);
      result[1] = (boundingBox.max[1] + y_length * random_perturb[1]);
      result[2] = (boundingBox.max[2] + z_length * random_perturb[2]);
      break;
  }
}

void RayTracer::generateRay1D(const double *rayend, double *result, int &binno) const {
  const auto &boundingBox = bounds_;

  result[0] = (rayend[0]);
  result[1] = (rayend[1]);
  result[2] = (rayend[2]);
  int ray_dir = (binning_type + 1) % 3;
  double bb_ray_dir_length = boundingBox.max[ray_dir] - boundingBox.min[ray_dir];
  result[ray_dir] = (boundingBox.min[ray_dir] - bb_ray_dir_length);

  const int no_bins_1d = bins_1d.size();
  double bin_width = (bounds_.max[binning_type] - bounds_.min[binning_type]) / no_bins_1d;
  binno = static_cast<int>((rayend[binning_type] - bounds_.min[binning_type]) / bin_width);
}

void RayTracer::generateRay2D(const double *rayend, double *result, unsigned int binno[2]) const {
  const auto &boundingBox = bounds_;

  result[0] = (rayend[0]);
  result[1] = (rayend[1]);
  result[2] = (rayend[2]);
  int ray_dir = (binning_type + 2) % 3;
  int bin_dir_1 = (binning_type) % 3;
  int bin_dir_2 = (binning_type + 1) % 3;
  double bb_ray_dir_length = boundingBox.max[ray_dir] - boundingBox.min[ray_dir];
  result[ray_dir] = (boundingBox.min[ray_dir] - bb_ray_dir_length);

  const int no_bins_dir1 = bins_2d.size();
  const int no_bins_dir2 = bins_2d[0].size();
  double bin_width_dir1 = (bounds_.max[bin_dir_1] - bounds_.min[bin_dir_1]) / no_bins_dir1;
  double bin_width_dir2 = (bounds_.max[bin_dir_2] - bounds_.min[bin_dir_2]) / no_bins_dir2;
  binno[0] = static_cast<unsigned int>((rayend[bin_dir_1] - bounds_.min[bin_dir_1]) / bin_width_dir1);
  binno[1] = static_cast<unsigned int>((rayend[bin_dir_2] - bounds_.min[bin_dir_2]) / bin_width_dir2);
}

Box3d RayTracer::calcBoundingBox(const RayTracer::TriangleData *triangles, const size_t numTriangles) {
  Box3d box;
  for (int d = 0; d < 3; d++) {
    box.min[d] = triangles[0].v[0].data[d];
    box.max[d] = box.min[d];
  }

  for (unsigned int i = 0; i < numTriangles; i++) {
    for (const auto &vert: triangles[i].v) {
      for (int dim = 0; dim < 3; dim++) {
        box.min[dim] = std::min(box.min[dim], vert.data[dim]);
        box.max[dim] = std::max(box.max[dim], vert.data[dim]);
      }
    }
  }

  double eps = 1e-6;
  for (int d = 0; d < 3; d++) {
    box.min[d] -= eps;
    box.max[d] += eps;
  }

  return box;
}

#ifdef BENCHMARK
InsideType::Type RayTracer::ifInside(const double *point, bool forceNoBin) {
#else
InsideType::Type RayTracer::ifInside(const double *point, bool forceNoBin) const {
#endif

#ifdef BENCHMARK
  high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
#endif
  if ((binning_type == BinningType::Type::X_1D ||
      binning_type == BinningType::Type::Y_1D ||
      binning_type == BinningType::Type::Z_1D) and (not forceNoBin)) {
    double ray_start[3];
    int binno;
    generateRay1D(point, ray_start, binno);
    InsideType::Type res = ifInside(ray_start, point, binno);
    if (res == InsideType::UNKNOWN) {
      /// falling back to original ray-tracing
      for (int i = 0; i < 8; i++) {
        res = ifInside(ray_starts_[i].data, point);
        if (res == InsideType::UNKNOWN) {
          std::cout << "Changing ray to no: " << i << std::endl;
          continue;
        }
#ifdef BENCHMARK
        high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        d_in_out.first += duration_cast<duration<double>>(t2 - t0);
        d_in_out.second += (i + 1);
#endif
        return res;
      }
#ifdef BENCHMARK
      high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
      d_in_out.first += duration_cast<duration<double>>(t2 - t0);
      d_in_out.second += 1;
#endif
      return InsideType::UNKNOWN;
    }
#ifdef BENCHMARK
    high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    d_in_out.first += duration_cast<duration<double>>(t2 - t0);
    d_in_out.second++;
#endif
    return res;
  } else if (binning_type == BinningType::Type::NOBIN or forceNoBin) {
    for (int i = 0; i < 8; i++) {
      InsideType::Type res = ifInside(ray_starts_[i].data, point);
      if (res == InsideType::UNKNOWN) {
        std::cout << "Changing ray to no: " << i << std::endl;
        continue;
      }
#ifdef BENCHMARK
      high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
      d_in_out.first += duration_cast<duration<double>>(t1 - t0);
      d_in_out.second += (i + 1);
#endif
      return res;
    }
#ifdef BENCHMARK
    high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    d_in_out.first += duration_cast<duration<double>>(t1 - t0);
    d_in_out.second += 1;
#endif
    return InsideType::UNKNOWN;
  } else if (not forceNoBin) {
    double ray_start[3];
    unsigned int binno[2];
    generateRay2D(point, ray_start, binno);
    InsideType::Type res = ifInsideBin2D(ray_start, point, binno);
    if (res == InsideType::UNKNOWN) {
      /// falling back to original ray-tracing
      for (int i = 0; i < 8; i++) {
        res = ifInside(ray_starts_[i].data, point);
        if (res == InsideType::UNKNOWN) {
          std::cout << "Changing ray to no: " << i << std::endl;
          continue;
        }
#ifdef BENCHMARK
        high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        d_in_out.first += duration_cast<duration<double>>(t2 - t0);
        d_in_out.second += (i + 1);
#endif
        return res;
      }
#ifdef BENCHMARK
      high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
      d_in_out.first += duration_cast<duration<double>>(t2 - t0);
      d_in_out.second += 1;
#endif
      return InsideType::UNKNOWN;
    }
#ifdef BENCHMARK
    high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    d_in_out.first += duration_cast<duration<double>>(t2 - t0);
    d_in_out.second++;
#endif
    return res;
  }
}

RayTracer::RayTracer(const TriangleData *triangleData, const size_t numTriangles)
    : triangles_(triangleData), numTriangles_(numTriangles) {
  bounds_ = calcBoundingBox(triangles_, numTriangles);
  for (int i = 0; i < 8; i++) {
    generateRay(i, ray_starts_[i].data);
  }
  setBinning();
}

RayTracer::~RayTracer() {
}


//void RayTracer::setTriangles(std::vector<stl::Triangle> &lines) {
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
//#endif
//  triangles_ = new TriangleData[lines.size()];
//  numTriangles_ = lines.size();
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//#endif
//  for (unsigned int i = 0; i < lines.size(); i++) {
//    const auto &from = lines[i];
//    auto &to = triangles_[i];
//    for (int vert = 0; vert < 3; vert++)
//      to.v[vert] = from.v[vert];
////    to.boundingSphere = Sphere::fromTriangle(to.v);
//  }
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
//  d_bounding_sphere = duration_cast<duration<double>>(t2 - t1);
//  high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
//#endif
//  bounds_ = calcBoundingBox(triangles_,numTriangles);
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
//  d_bounding_box = duration_cast<duration<double>>(t4 - t3);
//#endif
//  for (int i = 0; i < 8; i++) {
//    generateRay(i, ray_starts_[i].data);
//  }
//  setBinning();
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t5 = std::chrono::high_resolution_clock::now();
//  d_set_triangle = duration_cast<duration<double>>(t5 - t0);
//#endif
//
//}

//void RayTracer::setTriangles(std::vector<RayTracer::TriangleData> &&lines) {
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
//#endif
//  triangles_ = lines;
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//#endif
////  for (unsigned int i = 0; i < lines.size(); i++) {
////    triangles_[i].boundingSphere = Sphere::fromTriangle(triangles_[i].v);
////  }
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
//  d_bounding_sphere = duration_cast<duration<double>>(t2 - t1);
//  high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
//#endif
//  bounds_ = calcBoundingBox(triangles_);
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
//  d_bounding_box = duration_cast<duration<double>>(t4 - t3);
//#endif
//  for (int i = 0; i < 8; i++) {
//    generateRay(i, ray_starts_[i].data);
//  }
//  setBinning();
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t5 = std::chrono::high_resolution_clock::now();
//  d_set_triangle = duration_cast<duration<double>>(t5 - t0);
//#endif
//}

void RayTracer::setBinning() {
  const size_t no_tris = numTriangles_;
  double x_range = bounds_.max[0] - bounds_.min[0];
  double y_range = bounds_.max[1] - bounds_.min[1];
  double z_range = bounds_.max[2] - bounds_.min[2];

  binning_type = BinningType::Type::YZ_2D;
  if (x_range <= y_range and x_range <= z_range and binning_type == BinningType::Type::NOBIN) {
    binning_type = BinningType::Type::YZ_2D;
  }
  if (y_range <= x_range and y_range <= z_range and binning_type == BinningType::Type::NOBIN) {
    binning_type = BinningType::Type::ZX_2D;
  }
  if (z_range <= y_range and z_range <= x_range and binning_type == BinningType::Type::NOBIN) {
    binning_type = BinningType::Type::XY_2D;
  }
  binning_2d(int(sqrt(no_tris / 10)), int(sqrt(no_tris / 10)));
//  binning_nobin();
}

void RayTracer::binning_nobin() {
  bins_1d.resize(1);
  bins_1d[0].first.resize(2);
  // lower and upper bound with overlap of 2 * eps
  bins_1d[0].first[0] = bounds_.min[0] - 0.01;
  bins_1d[0].first[1] = bounds_.max[0] + 0.01;
  bins_1d[0].second.resize(numTriangles_);
  std::iota(std::begin(bins_1d[0].second), std::end(bins_1d[0].second), 0);
}

void RayTracer::binning_1d(int no_bins_1d = 1) {
  // classify lines to bins
  int bin_dir_1d = -1;
  if (binning_type == BinningType::Type::X_1D ||
      binning_type == BinningType::Type::Y_1D ||
      binning_type == BinningType::Type::Z_1D) {
    bin_dir_1d = binning_type;
  }
  const double bin_width = (bounds_.max[bin_dir_1d] - bounds_.min[bin_dir_1d]) / no_bins_1d;
  const double eps_percent = 0.01;
  const double eps = eps_percent * bin_width;
  bins_1d.resize(no_bins_1d);
  for (unsigned int i = 0; i < no_bins_1d; i++) {
    bins_1d[i].first.resize(2);
    // lower and upper bound with overlap of 2 * eps
    bins_1d[i].first[0] = bounds_.min[bin_dir_1d] + bin_width * i - eps;
    bins_1d[i].first[1] = bounds_.min[bin_dir_1d] + bin_width * (i + 1) + eps;
    bins_1d[i].second.reserve(numTriangles_ / no_bins_1d * 1.5);
  }

  for (unsigned int triID = 0; triID < numTriangles_; triID++) {
    const auto &tri = triangles_[triID];
    double c0 = tri.v[0].data[bin_dir_1d], c1 = tri.v[1].data[bin_dir_1d], c2 = tri.v[2].data[bin_dir_1d];
    double cmin = std::min(std::min(c0, c1), c2);
    double cmax = std::max(std::max(c0, c1), c2);

    auto cbin_min = static_cast<unsigned int>((cmin - eps - bounds_.min[bin_dir_1d]) / bin_width);
    auto cbin_max = static_cast<unsigned int>((cmax + eps - bounds_.min[bin_dir_1d]) / bin_width);
    if (cbin_max == no_bins_1d) {
      cbin_max--;
    }
    assert(cmin > bins_1d[cbin_min].first[0]);
    assert(cmin < bins_1d[cbin_min].first[1]);

    assert(cmax > bins_1d[cbin_max].first[0]);
    assert(cmax < bins_1d[cbin_max].first[1]);

    assert(cbin_min <= cbin_max);
    for (unsigned int bin_no = cbin_min; bin_no <= cbin_max; bin_no++) {
      bins_1d[bin_no].second.emplace_back(triID);
    }
  }
}

void RayTracer::binning_2d(int no_bins_dir_1 = 1, int no_bins_dir_2 = 1) {
  // classify lines to bins
  int bin_dir_1 = -1, bin_dir_2 = -1;
  if (binning_type == BinningType::Type::XY_2D) {
    bin_dir_1 = 0;
    bin_dir_2 = 1;
  } else if (binning_type == BinningType::Type::YZ_2D) {
    bin_dir_1 = 1;
    bin_dir_2 = 2;
  } else if (binning_type == BinningType::Type::ZX_2D) {
    bin_dir_1 = 2;
    bin_dir_2 = 0;
  }
  const double bin_width_dir1 = (bounds_.max[bin_dir_1] - bounds_.min[bin_dir_1]) / no_bins_dir_1;
  const double bin_width_dir2 = (bounds_.max[bin_dir_2] - bounds_.min[bin_dir_2]) / no_bins_dir_2;
  const double eps_percent = 0.01;
  const double eps_dir1 = eps_percent * bin_width_dir1;
  const double eps_dir2 = eps_percent * bin_width_dir2;
  bins_2d.clear();
  bins_2d.resize(no_bins_dir_1);
  for (auto &bins_2d_dir2 : bins_2d) {
    bins_2d_dir2.resize(no_bins_dir_2);
  }

  for (unsigned int i = 0; i < no_bins_dir_1; i++) {
    for (unsigned int j = 0; j < no_bins_dir_2; j++) {
      bins_2d[i][j].first.resize(4); // x-, x+, y-, y+
      // lower and upper bound with overlap of 2 * eps
      bins_2d[i][j].first[0] = bounds_.min[bin_dir_1] + bin_width_dir1 * i - eps_dir1;
      bins_2d[i][j].first[1] = bounds_.min[bin_dir_1] + bin_width_dir1 * (i + 1) + eps_dir1;
      bins_2d[i][j].first[2] = bounds_.min[bin_dir_2] + bin_width_dir2 * j - eps_dir2;
      bins_2d[i][j].first[3] = bounds_.min[bin_dir_2] + bin_width_dir2 * (j + 1) + eps_dir2;
      bins_2d[i][j].second.reserve(numTriangles_ / no_bins_dir_1 / no_bins_dir_2 * 4);
    }
  }

  for (unsigned int triID = 0; triID < numTriangles_; triID++) {
    const auto &tri = triangles_[triID];
    double c0_dir1 = tri.v[0].data[bin_dir_1], c1_dir1 = tri.v[1].data[bin_dir_1], c2_dir1 = tri.v[2].data[bin_dir_1];
    double cmin_dir1 = std::min(std::min(c0_dir1, c1_dir1), c2_dir1);
    double cmax_dir1 = std::max(std::max(c0_dir1, c1_dir1), c2_dir1);

    double c0_dir2 = tri.v[0].data[bin_dir_2], c1_dir2 = tri.v[1].data[bin_dir_2], c2_dir2 = tri.v[2].data[bin_dir_2];
    double cmin_dir2 = std::min(std::min(c0_dir2, c1_dir2), c2_dir2);
    double cmax_dir2 = std::max(std::max(c0_dir2, c1_dir2), c2_dir2);

    auto cbin_min_dir1 = static_cast<unsigned int>((cmin_dir1 - eps_dir1 - bounds_.min[bin_dir_1]) / bin_width_dir1);
    auto cbin_max_dir1 = static_cast<unsigned int>((cmax_dir1 + eps_dir1 - bounds_.min[bin_dir_1]) / bin_width_dir1);
    auto cbin_min_dir2 = static_cast<unsigned int>((cmin_dir2 - eps_dir2 - bounds_.min[bin_dir_2]) / bin_width_dir2);
    auto cbin_max_dir2 = static_cast<unsigned int>((cmax_dir2 + eps_dir2 - bounds_.min[bin_dir_2]) / bin_width_dir2);

    if (cbin_max_dir1 == no_bins_dir_1) {
      cbin_max_dir1--;
    }
    if (cbin_max_dir2 == no_bins_dir_2) {
      cbin_max_dir2--;
    }

    for (unsigned int i = cbin_min_dir1; i <= cbin_max_dir1; i++) {
      for (unsigned int j = cbin_min_dir2; j <= cbin_max_dir2; j++) {
        bins_2d[i][j].second.emplace_back(triID);
      }
    }
  }

  // generate a typical bin-1d as backup
  std::pair<std::vector<double>, std::vector<unsigned int>> bin0;
  bins_1d.push_back(bin0);
  bins_1d[0].second.resize(numTriangles_);
  std::iota(std::begin(bins_1d[0].second), std::end(bins_1d[0].second), 0);

}

#ifdef BENCHMARK
void RayTracer::printTimer(const std::string &filename) const {
  std::ofstream logfile;
  if (not(filename.empty())) {
    logfile.open((filename + ".log").c_str(), std::ios::out);
    int sum = 0;
    if (binning_type == BinningType::Type::X_1D ||
        binning_type == BinningType::Type::Y_1D ||
        binning_type == BinningType::Type::Z_1D ||
        binning_type == BinningType::Type::NOBIN) {
      logfile << "bins: ";
      for (const auto &bin: bins_1d) {
        sum += bin.second.size();
        logfile << bin.second.size() << " ";
      }
      logfile << "\nTotal triangles: " << getNumTriangles() << ", sum of triangles in " << bins_1d.size() << " bins: " << sum << "\n";
    } else {
      logfile << "bins2D: ";
      for (const auto &bin_dir1: bins_2d) {
        for (const auto &bin: bin_dir1) {
          sum += bin.second.size();
          logfile << bin.second.size() << " ";
        }
      }
      logfile << "\nTotal triangles: " << getNumTriangles()
              << ", sum of triangles in " << bins_2d.size() * bins_2d[0].size() << " bins: " << sum << "\n";
    }

    logfile << "Generate ray-tracer: " << d_set_triangle.count() << "\n";
    logfile << "  - calc bounding spheres: " << d_bounding_sphere.count() << "\n";
    logfile << "  - calc bounding box: " << d_bounding_box.count() << "\n";
    logfile << "  - overhead: " << d_set_triangle.count() - d_bounding_sphere.count() - d_bounding_box.count() << "\n";
    logfile << "Testing points: " << d_in_out.second
            << " --- " << d_in_out.first.count()
            << " --- " << d_in_out.second / d_in_out.first.count() << " points/sec.\n";
    logfile << "  - out bbox: " << d_out_bounding_box.second <<
            " --- " << d_out_bounding_box.first.count() <<
            " --- " << d_out_bounding_box.second / d_out_bounding_box.first.count() << " points/sec.\n";
  }
  std::cout << "Generate ray-tracer: " << d_set_triangle.count() << "\n";
  std::cout << "  - calc bounding spheres: " << d_bounding_sphere.count() << "\n";
  std::cout << "  - calc bounding box: " << d_bounding_box.count() << "\n";
  std::cout << "  - overhead: " << d_set_triangle.count() - d_bounding_sphere.count() - d_bounding_box.count() << "\n";
  std::cout << "Testing points: " << d_in_out.second
            << " --- " << d_in_out.first.count()
            << " --- " << d_in_out.second / d_in_out.first.count() << " points/sec.\n";
  std::cout << "  - out bbox: " << d_out_bounding_box.second <<
            " --- " << d_out_bounding_box.first.count() <<
            " --- " << d_out_bounding_box.second / d_out_bounding_box.first.count() << " points/sec.\n";
  assert(d_in_bounding_boxs.size() == d_out_triangle_spheres.size());
  assert(d_line_tris.size() == d_out_triangle_spheres.size());
  assert(d_out_raybbs.size() == d_out_triangle_spheres.size());
  double inbbox_time = 0.0, out_raybb_time = 0.0, out_sphere_time = 0.0, out_line_tris_time = 0.0;
  unsigned long inbbox_points = d_in_bounding_boxs.size();
  unsigned long inbbox_points_ave_out_raybb = 0;
  unsigned long inbbox_points_ave_out_triangle_sphere = 0;
  unsigned long inbbox_points_ave_line_tris = 0;
  for (unsigned int i = 0; i < d_in_bounding_boxs.size(); i++) {
    assert(d_in_bounding_boxs[i].second - d_out_triangle_spheres[i].second - d_line_tris[i].second - d_out_raybbs[i].second == 0);
    inbbox_time += d_in_bounding_boxs[i].first.count();
    inbbox_points_ave_out_raybb += d_out_raybbs[i].second;
    inbbox_points_ave_out_triangle_sphere += d_out_triangle_spheres[i].second;
    inbbox_points_ave_line_tris += d_line_tris[i].second;
    out_raybb_time += d_out_raybbs[i].first.count();
    out_sphere_time += d_out_triangle_spheres[i].first.count();
    out_line_tris_time += d_line_tris[i].first.count();
  }
  if (not(filename.empty())) {
    logfile << "  - in bbox: " << inbbox_points
            << " --- " << inbbox_time
            << " --- " << inbbox_points / inbbox_time << " points/sec.\n";
    logfile << "    - outside raybb: " << inbbox_points_ave_out_raybb
            << " --- " << out_raybb_time
            << " --- " << inbbox_points_ave_out_raybb / out_raybb_time << " rays/sec"
            << " --- " << inbbox_points_ave_out_raybb / (out_raybb_time * inbbox_points) << " points/sec.\n";
    logfile << "    - outside bbsphere: " << inbbox_points_ave_out_triangle_sphere
            << " --- " << out_sphere_time
            << " --- " << inbbox_points_ave_out_triangle_sphere / out_sphere_time << " rays/sec"
            << " --- " << inbbox_points_ave_out_triangle_sphere / (out_sphere_time * inbbox_points) << " points/sec.\n";
    logfile << "    - line-triangle: " << inbbox_points_ave_line_tris
            << " --- " << out_line_tris_time
            << " --- " << inbbox_points_ave_line_tris / out_line_tris_time << " rays/sec"
            << " --- " << inbbox_points_ave_line_tris / (out_line_tris_time * inbbox_points) << " points/sec.\n";
  }
  logfile.close();
  std::cout << "  - in bbox: " << inbbox_points
            << " --- " << inbbox_time
            << " --- " << inbbox_points / inbbox_time << " points/sec.\n";
  std::cout << "    - outside raybb: " << inbbox_points_ave_out_raybb
            << " --- " << out_raybb_time
            << " --- " << inbbox_points_ave_out_raybb / out_raybb_time << " rays/sec"
            << " --- " << inbbox_points_ave_out_raybb / (out_raybb_time * inbbox_points) << " points/sec.\n";
  std::cout << "    - outside bbsphere: " << inbbox_points_ave_out_triangle_sphere
            << " --- " << out_sphere_time
            << " --- " << inbbox_points_ave_out_triangle_sphere / out_sphere_time << " rays/sec"
            << " --- " << inbbox_points_ave_out_triangle_sphere / (out_sphere_time * inbbox_points) << " points/sec.\n";
  std::cout << "    - line-triangle: " << inbbox_points_ave_line_tris
            << " --- " << out_line_tris_time
            << " --- " << inbbox_points_ave_line_tris / out_line_tris_time << " rays/sec"
            << " --- " << inbbox_points_ave_line_tris / (out_line_tris_time * inbbox_points) << " points/sec.\n";
}

void RayTracer::cleanUpTimer() {
  d_bounding_sphere = std::chrono::milliseconds::zero();
  d_bounding_box = std::chrono::milliseconds::zero();
  d_set_triangle = std::chrono::milliseconds::zero();
  d_in_out.first = std::chrono::milliseconds::zero();
  d_in_out.second = 0;
  d_out_bounding_box.first = std::chrono::milliseconds::zero();
  d_out_bounding_box.second = 0;
  /// for each point tested, we have a breakdown of performance, time/no. of lines
  d_in_bounding_boxs.clear();
  d_out_raybbs.clear();
  d_out_triangle_spheres.clear();
  d_line_tris.clear();
}
#endif

