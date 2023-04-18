//
// Created by boshun on 7/13/20.
//

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include "raytracer/ray_tracer_2d.h"

bool RayTracer2D::OutsideBoundingCircle(const double *orig, const double *dir, const double *center, const double radius) {
  double w[2], b, Pb[2], distance_square;
  Vec_Sub(w, center, orig);
  b = Vec_Dot(w, dir) / Vec_Dot(dir, dir);
  Pb[0] = orig[0] + b * dir[0];
  Pb[1] = orig[1] + b * dir[1];
  distance_square = (Pb[0] - center[0]) * (Pb[0] - center[0]) +
      (Pb[1] - center[1]) * (Pb[1] - center[1]);
  return (distance_square > radius * radius);
}

bool RayTracer2D::OutsideRayBBOX(const double *orig, const double *dir, const double *center, const double radius) {

  double lb[2], ub[2];
  lb[0] = std::min(orig[0], orig[0] + dir[0]) - radius;
  lb[1] = std::min(orig[1], orig[1] + dir[1]) - radius;
  ub[0] = std::max(orig[0], orig[0] + dir[0]) + radius;
  ub[1] = std::max(orig[1], orig[1] + dir[1]) + radius;
  return center[0] < lb[0] || center[1] < lb[1] || center[0] > ub[0] || center[1] > ub[1];
}

IntersectionType2D::Type
RayTracer2D::ifLineIntersectLine(const double *orig, const double *dest, const double *vert0, const double *vert1) {
  // Line AB represented as a1x + b1y = c1
  double a1 = dest[1] - orig[1];
  double b1 = orig[0] - dest[0];
  double c1 = a1 * (orig[0]) + b1 * (orig[1]);

  // Line CD represented as a2x + b2y = c2
  double a2 = vert0[1] - vert1[1];
  double b2 = vert1[0] - vert0[0];
  double c2 = a2 * (vert1[0]) + b2 * (vert1[1]);

  double det = a1 * b2 - a2 * b1;
  const double epsilon_node = 1e-7;
  const double epsilon_edge = 1e-7;
  if (det < 1e-10 && det > -1e-10) {
    /// Ray parallel to the line
    return IntersectionType2D::PARALLEL;
  } else {
    double x = (b2 * c1 - b1 * c2) / det;
    double y = (a1 * c2 - a2 * c1) / det;

    /// Judge if the intersection point is inside the line
    /// Outside the line && outside reach
    if (x > std::max(orig[0], dest[0]) + epsilon_edge ||
        y > std::max(orig[1], dest[1]) + epsilon_edge ||
        x > std::max(vert0[0], vert1[0]) + epsilon_node ||
        y > std::max(vert0[1], vert1[1]) + epsilon_node ||
        x < std::min(orig[0], dest[0]) - epsilon_edge ||
        y < std::min(orig[1], dest[1]) - epsilon_edge ||
        x < std::min(vert0[0], vert1[0]) - epsilon_node ||
        y < std::min(vert0[1], vert1[1]) - epsilon_node) {
      return IntersectionType2D::NO_INTERSECTION;
    } else {
      /// On line edge (possible outcome are on edge or on surface)
      if ((fabs(x - vert0[0]) < epsilon_node && fabs(y - vert0[1]) < epsilon_node) ||
          (fabs(x - vert1[0]) < epsilon_node && fabs(y - vert1[1]) < epsilon_node)) {
        /// Judge t (t already satisfies condition)
        if (fabs(x - dest[0]) < epsilon_edge &&
            fabs(x - dest[1]) < epsilon_edge) {
          return IntersectionType2D::ON_LINE;
        } else {
          return IntersectionType2D::ON_NODE;
        }
      } else {
        /// Inside line (possible outcome are on surface or intersect)
        if (fabs(x - dest[0]) < epsilon_edge &&
            fabs(x - dest[1]) < epsilon_edge) {
          return IntersectionType2D::ON_LINE;
        } else {
          return IntersectionType2D::INTERSECT;
        }
      }
    }
  }

}

bool RayTracer2D::ifOutsideBoundingBox(const double *ray_end, const Box2d &bounds) {
  return (ray_end[0] < bounds.min[0]) || (ray_end[0] > bounds.max[0])
      || (ray_end[1] < bounds.min[1]) || (ray_end[1] > bounds.max[1]);
}

#ifdef BENCHMARK
InsideType2D::Type RayTracer2D::ifInside(const double *ray_start, const double *ray_end, const int binno) {
#else
  InsideType2D::Type RayTracer2D::ifInside(const double *ray_start, const double *ray_end, const int binno) const {
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
    return InsideType2D::OUTSIDE;
  }

  int noOfIntersections = 0;
  double dir[2];
  Vec_Sub(dir, ray_end, ray_start);
#ifdef BENCHMARK
  std::pair<duration<double>, unsigned long> d_in_bounding_box;
  std::pair<duration<double>, unsigned long> d_out_raybb;
  std::pair<duration<double>, unsigned long> d_out_line_sphere;
  std::pair<duration<double>, unsigned long> d_line_line;
  t0 = std::chrono::high_resolution_clock::now();
  duration<double> duration_line_line = std::chrono::milliseconds::zero();
#endif

  std::vector<unsigned int> lineIDs;
  if (binno < 0) {
    lineIDs.clear();
    lineIDs.resize(numLines_);
    std::iota(std::begin(lineIDs), std::end(lineIDs), 0);
  } else {
    lineIDs = bins_1d[binno].second;
  }

  for (const auto &lineID: lineIDs) {
    const auto &line = lines_[lineID];

#ifdef BENCHMARK
    high_resolution_clock::time_point t_out_s = std::chrono::high_resolution_clock::now();
#endif


    ///  Determine if the ray intersects with the bounding sphere or a line, skipped if not
    if (OutsideRayBBOX(ray_start, dir, line.boundingCircle().center, line.boundingCircle().radius)) {
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
    ///  Determine if the ray intersects with the bounding sphere or a line, skipped if not
    if (OutsideBoundingCircle(ray_start, dir, line.boundingCircle().center, line.boundingCircle().radius)) {
#ifdef BENCHMARK
      d_in_bounding_box.second++;
      high_resolution_clock::time_point t_out_sphere_e = std::chrono::high_resolution_clock::now();
      d_out_line_sphere.second++;
      d_out_line_sphere.first += duration_cast<duration<double>>(t_out_sphere_e - t_out_s);
#endif
      continue;
    }

#ifdef BENCHMARK
    d_in_bounding_box.second++;
    d_line_line.second++;
#endif
    high_resolution_clock::time_point t_s_line_line = std::chrono::high_resolution_clock::now();
//    IntersectionType::Type result_from_ifIntersect = IntersectionType::NO_INTERSECTION;
    IntersectionType2D::Type result_from_ifIntersect =
        ifLineIntersectLine(ray_start, ray_end, line.v[0].data, line.v[1].data);
//    std::cout << lineID << " " << ray_start[0] << " " << ray_start[1] << " " << ray_end[0] << " " << ray_end[1] << " "
//              << line.v[0].data[0] << " " << line.v[0].data[1] << " " << line.v[1].data[0] << " " << line.v[1].data[1] << " "
//              << result_from_ifIntersect << "\n";
    /// If ON_EDGE or ON_SURFACE or ON_STARTPOINT happens, then this position is vague, jump out of the loop and choose another ray
    if (result_from_ifIntersect == IntersectionType2D::ON_NODE ||
        result_from_ifIntersect == IntersectionType2D::ON_STARTPOINT) {
#ifdef BENCHMARK
      high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
      d_in_bounding_box.first += duration_cast<duration<double>>(t1 - t0);
      d_out_raybbs.emplace_back(d_out_raybb);
      d_out_line_circles.emplace_back(d_out_line_sphere);
      d_in_bounding_boxs.emplace_back(d_in_bounding_box);
      high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
      d_line_line.first += duration_cast<duration<double>>(t2 - t_s_line_line);
      d_line_lines.emplace_back(d_line_line);
#endif
      return InsideType2D::UNKNOWN;
    }
    if (result_from_ifIntersect == IntersectionType2D::ON_LINE) {
#ifdef BENCHMARK
      high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
      d_in_bounding_box.first += duration_cast<duration<double>>(t1 - t0);
      d_out_raybbs.emplace_back(d_out_raybb);
      d_out_line_circles.emplace_back(d_out_line_sphere);
      d_in_bounding_boxs.emplace_back(d_in_bounding_box);
      high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
      d_line_line.first += duration_cast<duration<double>>(t2 - t_s_line_line);
      d_line_lines.emplace_back(d_line_line);
#endif
      return InsideType2D::ON_SURFACE;
    }
    if (result_from_ifIntersect == IntersectionType2D::INTERSECT) {
      noOfIntersections++;
    }
#ifdef BENCHMARK
    high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    duration_line_line += duration_cast<duration<double>>(t2 - t_s_line_line);
#endif
  }
#ifdef BENCHMARK
  high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  d_line_line.first += duration_line_line;
  d_in_bounding_box.first += duration_cast<duration<double>>(t1 - t0);
  d_out_raybbs.emplace_back(d_out_raybb);
  d_out_line_circles.emplace_back(d_out_line_sphere);
  d_in_bounding_boxs.emplace_back(d_in_bounding_box);
  d_line_lines.emplace_back(d_line_line);
#endif
  if ((noOfIntersections % 2) != 0) {
    return InsideType2D::INSIDE;
  } else {
    return InsideType2D::OUTSIDE;
  }
}

void RayTracer2D::generateRay(int pick, double *result) {
  const auto &boundingBox = bounds_;

  double random_perturb[2] = {
      1 + double(rand()) / RAND_MAX,
      1 + double(rand()) / RAND_MAX
  };
  double x_length = boundingBox.max[0] - boundingBox.min[0];
  double y_length = boundingBox.max[1] - boundingBox.min[1];

  switch (pick) {
    case 0:result[0] = (boundingBox.min[0] - x_length * random_perturb[0]);
      result[1] = (boundingBox.min[1] - y_length * random_perturb[1]);
      break;
    case 1:result[0] = (boundingBox.max[0] + x_length * random_perturb[0]);
      result[1] = (boundingBox.min[1] - y_length * random_perturb[1]);
      break;
    case 2:result[0] = (boundingBox.min[0] - x_length * random_perturb[0]);
      result[1] = (boundingBox.max[1] + y_length * random_perturb[1]);
      break;
    case 3:
    default:result[0] = (boundingBox.max[0] + x_length * random_perturb[0]);
      result[1] = (boundingBox.max[1] + y_length * random_perturb[1]);
      break;
  }
}

void RayTracer2D::generateRay1D(const double *rayend, double *result, int &binno) const {
  const auto &boundingBox = bounds_;

  result[0] = (rayend[0]);
  result[1] = (rayend[1]);
  int ray_dir = (binning_type + 1) % 2;
  double bb_ray_dir_length = boundingBox.max[ray_dir] - boundingBox.min[ray_dir];
  result[ray_dir] = (boundingBox.min[ray_dir] - bb_ray_dir_length);

  const int no_bins_1d = bins_1d.size();
  double bin_width = (bounds_.max[binning_type] - bounds_.min[binning_type]) / no_bins_1d;
  binno = static_cast<int>((rayend[binning_type] - bounds_.min[binning_type]) / bin_width);
}

Box2d RayTracer2D::calcBoundingBox(const RayTracer2D::Line2DData *lines, const size_t numLines) {
  Box2d box;
  for (int d = 0; d < 2; d++) {
    box.min[d] = lines[0].v[0].data[d];
    box.max[d] = box.min[d];
  }

  for (unsigned int i = 0; i < numLines; i++) {
    for (const auto &vert: lines[i].v) {
      for (int dim = 0; dim < 2; dim++) {
        box.min[dim] = std::min(box.min[dim], vert.data[dim]);
        box.max[dim] = std::max(box.max[dim], vert.data[dim]);
      }
    }
  }

  double eps = 1e-6;
  for (int d = 0; d < 2; d++) {
    box.min[d] -= eps;
    box.max[d] += eps;
  }

  return box;
}

#ifdef BENCHMARK
InsideType2D::Type RayTracer2D::ifInside(const double *point, bool forceNoBin) {
#else
  InsideType2D::Type RayTracer2D::ifInside(const double *point, bool forceNoBin) const {
#endif

#ifdef BENCHMARK
  high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
#endif
  if ((binning_type == BinningType2D::Type::X_1D ||
      binning_type == BinningType2D::Type::Y_1D) and (not forceNoBin)) {
    double ray_start[3];
    int binno;
    generateRay1D(point, ray_start, binno);
    InsideType2D::Type res = ifInside(ray_start, point, binno);
    if (res == InsideType2D::UNKNOWN) {
      /// falling back to original ray-tracing
      for (int i = 0; i < 4; i++) {
        res = ifInside(ray_starts_[i].data, point);
        if (res == InsideType2D::UNKNOWN) {
//          std::cout << "Changing ray to no: " << i << std::endl;
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

//      throw std::runtime_error("In Out test failed");
      return InsideType2D::UNKNOWN;
    }
#ifdef BENCHMARK
    high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    d_in_out.first += duration_cast<duration<double>>(t2 - t0);
    d_in_out.second++;
#endif
    return res;
  } else if (forceNoBin) {
    for (int i = 0; i < 4; i++) {
      InsideType2D::Type res = ifInside(ray_starts_[i].data, point);
      if (res == InsideType2D::UNKNOWN) {
//        std::cout << "Changing ray to no: " << i << std::endl;
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
    std::cout << "Points = " << point[0] << " " << point[1] << "\n";
    throw std::runtime_error("In Out test failed");

    return InsideType2D::UNKNOWN;
  }
}

RayTracer2D::RayTracer2D(const Line2DData *lineData, const size_t numLines)
    : lines_(lineData), numLines_(numLines) {

  bounds_ = calcBoundingBox(lines_, numLines_);
  for (int i = 0; i < 4; i++) {
    generateRay(i, ray_starts_[i].data);
  }
  setBinning();
}

RayTracer2D::~RayTracer2D() {
}

//void RayTracer2D::setLines(std::vector<LINE2D::Line> &lines) {
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
//#endif
//  lines_.resize(lines.size());
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//#endif
//  for (unsigned int i = 0; i < lines.size(); i++) {
//    const auto &from = lines[i];
//    auto &to = lines_[i];
//    for (int vert = 0; vert < 2; vert++)
//      to.v[vert] = from.v[vert];
//    to.boundingCircle = Circle::fromLine(to.v);
//  }
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
//  d_bounding_circle = duration_cast<duration<double>>(t2 - t1);
//  high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
//#endif
//  bounds_ = calcBoundingBox(lines_);
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
//  d_bounding_box = duration_cast<duration<double>>(t4 - t3);
//#endif
//  for (int i = 0; i < 4; i++) {
//    generateRay(i, ray_starts_[i].data);
//  }
//  setBinning();
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t5 = std::chrono::high_resolution_clock::now();
//  d_set_line = duration_cast<duration<double>>(t5 - t0);
//#endif
//
//}
//
//void RayTracer2D::setLines(std::vector<RayTracer2D::Line2DData> &&lines) {
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
//#endif
//  lines_ = lines;
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//#endif
//  for (unsigned int i = 0; i < lines.size(); i++) {
//    lines_[i].boundingCircle = Circle::fromLine(lines_[i].v);
//  }
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
//  d_bounding_circle = duration_cast<duration<double>>(t2 - t1);
//  high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
//#endif
//  bounds_ = calcBoundingBox(lines_);
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
//  d_bounding_box = duration_cast<duration<double>>(t4 - t3);
//#endif
//  for (int i = 0; i < 4; i++) {
//    generateRay(i, ray_starts_[i].data);
//  }
//  setBinning();
//#ifdef BENCHMARK
//  high_resolution_clock::time_point t5 = std::chrono::high_resolution_clock::now();
//  d_set_line = duration_cast<duration<double>>(t5 - t0);
//#endif
//}

void RayTracer2D::setBinning() {
  double x_range = bounds_.max[0] - bounds_.min[0];
  double y_range = bounds_.max[1] - bounds_.min[1];
  binning_type = BinningType2D::Type::Y_1D;
  if (x_range <= y_range and binning_type == BinningType2D::Type::NOBIN) {
    binning_type = BinningType2D::Type::Y_1D;
  }
  if (y_range <= x_range and binning_type == BinningType2D::Type::NOBIN) {
    binning_type = BinningType2D::Type::X_1D;
  }
  binning_1d(int(numLines_ / 50));
}

void RayTracer2D::binning_nobin() {
  bins_1d.resize(1);
  bins_1d[0].first.resize(2);
  // lower and upper bound with overlap of 2 * eps
  bins_1d[0].first[0] = bounds_.min[0] - 0.01;
  bins_1d[0].first[1] = bounds_.max[0] + 0.01;
  bins_1d[0].second.reserve(numLines_);
  for (unsigned int lineID = 0; lineID < numLines_; lineID++) {
    bins_1d[0].second.emplace_back(lineID);
  }
}

void RayTracer2D::binning_1d(int no_bins_1d = 1) {
  // classify lines to bins
  int bin_dir_1d = -1;
  if (binning_type == BinningType2D::Type::X_1D ||
      binning_type == BinningType2D::Type::Y_1D) {
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
    bins_1d[i].second.reserve(numLines_ / no_bins_1d * 1.5);
  }

  for (unsigned int lineID = 0; lineID < numLines_; lineID++) {
    const auto &line = lines_[lineID];
    double c0 = line.v[0].data[bin_dir_1d], c1 = line.v[1].data[bin_dir_1d];
    double cmin = std::min(c0, c1);
    double cmax = std::max(c0, c1);

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
      bins_1d[bin_no].second.emplace_back(lineID);
    }
  }
}

#ifdef BENCHMARK
void RayTracer2D::printTimer(const std::string &filename) const {
  std::ofstream logfile;
  if (not(filename.empty())) {
    logfile.open((filename + ".log").c_str(), std::ios::out);
    int sum = 0;
    if (binning_type == BinningType2D::Type::X_1D ||
        binning_type == BinningType2D::Type::Y_1D ||
        binning_type == BinningType2D::Type::NOBIN) {
      logfile << "bins: ";
      for (const auto &bin : bins_1d) {
        sum += bin.second.size();
        logfile << bin.second.size() << " ";
      }
      logfile << "\nTotal lines: " << getNumLines() << ", sum of lines in " << bins_1d.size() << " bins: " << sum << "\n";
    }

    logfile << "Generate ray-tracer: " << d_set_line.count() << "\n";
    logfile << "  - calc bounding spheres: " << d_bounding_circle.count() << "\n";
    logfile << "  - calc bounding box: " << d_bounding_box.count() << "\n";
    logfile << "  - overhead: " << d_set_line.count() - d_bounding_circle.count() - d_bounding_box.count() << "\n";
    logfile << "Testing points: " << d_in_out.second
            << " --- " << d_in_out.first.count()
            << " --- " << d_in_out.second / d_in_out.first.count() << " points/sec.\n";
    logfile << "  - out bbox: " << d_out_bounding_box.second <<
            " --- " << d_out_bounding_box.first.count() <<
            " --- " << d_out_bounding_box.second / d_out_bounding_box.first.count() << " points/sec.\n";
  }
  std::cout << "Generate ray-tracer: " << d_set_line.count() << "\n";
  std::cout << "  - calc bounding spheres: " << d_bounding_circle.count() << "\n";
  std::cout << "  - calc bounding box: " << d_bounding_box.count() << "\n";
  std::cout << "  - overhead: " << d_set_line.count() - d_bounding_circle.count() - d_bounding_box.count() << "\n";
  std::cout << "Testing points: " << d_in_out.second
            << " --- " << d_in_out.first.count()
            << " --- " << d_in_out.second / d_in_out.first.count() << " points/sec.\n";
  std::cout << "  - out bbox: " << d_out_bounding_box.second <<
            " --- " << d_out_bounding_box.first.count() <<
            " --- " << d_out_bounding_box.second / d_out_bounding_box.first.count() << " points/sec.\n";
  assert(d_in_bounding_boxs.size() == d_out_line_circles.size());
  assert(d_line_lines.size() == d_out_line_circles.size());
  assert(d_out_raybbs.size() == d_out_line_circles.size());
  double inbbox_time = 0.0, out_raybb_time = 0.0, out_circle_time = 0.0, out_line_line_time = 0.0;
  unsigned long inbbox_points = d_in_bounding_boxs.size();
  unsigned long inbbox_points_ave_out_raybb = 0;
  unsigned long inbbox_points_ave_out_line_circle = 0;
  unsigned long inbbox_points_ave_line_line = 0;
  for (unsigned int i = 0; i < d_in_bounding_boxs.size(); i++) {
    assert(d_in_bounding_boxs[i].second - d_out_line_circles[i].second - d_line_lines[i].second - d_out_raybbs[i].second == 0);
    inbbox_time += d_in_bounding_boxs[i].first.count();
    inbbox_points_ave_out_raybb += d_out_raybbs[i].second;
    inbbox_points_ave_out_line_circle += d_out_line_circles[i].second;
    inbbox_points_ave_line_line += d_line_lines[i].second;
    out_raybb_time += d_out_raybbs[i].first.count();
    out_circle_time += d_out_line_circles[i].first.count();
    out_line_line_time += d_line_lines[i].first.count();
  }
  if (not(filename.empty())) {
    logfile << "  - in bbox: " << inbbox_points
            << " --- " << inbbox_time
            << " --- " << inbbox_points / inbbox_time << " points/sec.\n";
    logfile << "    - outside raybb: " << inbbox_points_ave_out_raybb
            << " --- " << out_raybb_time
            << " --- " << inbbox_points_ave_out_raybb / out_raybb_time << " rays/sec"
            << " --- " << inbbox_points_ave_out_raybb / (out_raybb_time * inbbox_points) << " points/sec.\n";
    logfile << "    - outside bbcircle: " << inbbox_points_ave_out_line_circle
            << " --- " << out_circle_time
            << " --- " << inbbox_points_ave_out_line_circle / out_circle_time << " rays/sec"
            << " --- " << inbbox_points_ave_out_line_circle / (out_circle_time * inbbox_points) << " points/sec.\n";
    logfile << "    - line-line: " << inbbox_points_ave_line_line
            << " --- " << out_line_line_time
            << " --- " << inbbox_points_ave_line_line / out_line_line_time << " rays/sec"
            << " --- " << inbbox_points_ave_line_line / (out_line_line_time * inbbox_points) << " points/sec.\n";
  }
  logfile.close();
  std::cout << "  - in bbox: " << inbbox_points
            << " --- " << inbbox_time
            << " --- " << inbbox_points / inbbox_time << " points/sec.\n";
  std::cout << "    - outside raybb: " << inbbox_points_ave_out_raybb
            << " --- " << out_raybb_time
            << " --- " << inbbox_points_ave_out_raybb / out_raybb_time << " rays/sec"
            << " --- " << inbbox_points_ave_out_raybb / (out_raybb_time * inbbox_points) << " points/sec.\n";
  std::cout << "    - outside bbcircle: " << inbbox_points_ave_out_line_circle
            << " --- " << out_circle_time
            << " --- " << inbbox_points_ave_out_line_circle / out_circle_time << " rays/sec"
            << " --- " << inbbox_points_ave_out_line_circle / (out_circle_time * inbbox_points) << " points/sec.\n";
  std::cout << "    - line-line: " << inbbox_points_ave_line_line
            << " --- " << out_line_line_time
            << " --- " << inbbox_points_ave_line_line / out_line_line_time << " rays/sec"
            << " --- " << inbbox_points_ave_line_line / (out_line_line_time * inbbox_points) << " points/sec.\n";
}

void RayTracer2D::cleanUpTimer() {
  d_bounding_circle = std::chrono::milliseconds::zero();
  d_bounding_box = std::chrono::milliseconds::zero();
  d_set_line = std::chrono::milliseconds::zero();
  d_in_out.first = std::chrono::milliseconds::zero();
  d_in_out.second = 0;
  d_out_bounding_box.first = std::chrono::milliseconds::zero();
  d_out_bounding_box.second = 0;
  /// for each point tested, we have a breakdown of performance, time/no. of lines
  d_in_bounding_boxs.clear();
  d_out_raybbs.clear();
  d_out_line_circles.clear();
  d_line_lines.clear();
}

#endif
