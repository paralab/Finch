//
// Created by boshun on 7/5/20.
//

#ifndef RAY_TRACING_TESTPOINTS_H
#define RAY_TRACING_TESTPOINTS_H

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "raytracer/stl_reader.h"
#include "raytracer/gmsh_2d_reader.h"
#include "raytracer/ray_tracer.h"
#include "raytracer/ray_tracer_2d.h"

struct test_point {
  double xyz[3];
  InsideType::Type flag;
  test_point() : xyz{}, flag(InsideType::UNKNOWN) {}
  test_point(float xp, float yp, float zp) : xyz{xp, yp, zp}, flag(InsideType::UNKNOWN) {}
};

struct test_point_2d {
  double xy[2];
  InsideType2D::Type flag;
  test_point_2d() : xy{}, flag(InsideType2D::UNKNOWN) {}
  test_point_2d(float xp, float yp) : xy{xp, yp}, flag(InsideType2D::UNKNOWN) {}
};

/**
 * Generate a grid of test points based on the raytracer's bbox,
 * @param raytracer
 * @param Nxyz
 * @param points
 */
void grid_test_points(const RayTracer &raytracer, std::vector<int> Nxyz, std::vector<test_point> &points) {
  points.clear();
  assert(Nxyz.size() == 3);
  std::vector<double> center(3, 0), lb(3, 0), ub(3, 0);
  for (int dim = 0; dim < 3; dim++) {
    center[dim] = (raytracer.bounds().min[dim] + raytracer.bounds().max[dim]) / 2;
    lb[dim] = center[dim] - (center[dim] - raytracer.bounds().min[dim]) * 1.2;
    ub[dim] = center[dim] + (raytracer.bounds().max[dim] - center[dim]) * 1.2;
  }
  points.reserve(Nxyz[0] * Nxyz[1] * Nxyz[2]);
  for (int i = 0; i < Nxyz[0]; i++) {
    for (int j = 0; j < Nxyz[1]; j++) {
      for (int k = 0; k < Nxyz[2]; k++) {
//        test_point temp(i + double(rand()) / RAND_MAX, j + double(rand()) / RAND_MAX, k + double(rand()) / RAND_MAX);
        test_point temp(static_cast<float>(lb[0] + i * (ub[0] - lb[0]) / Nxyz[0]),
                        static_cast<float>(lb[1] + j * (ub[1] - lb[1]) / Nxyz[1]),
                        static_cast<float>(lb[2] + k * (ub[2] - lb[2]) / Nxyz[2]));
        points.emplace_back(temp);
      }
    }
  }
}

/**
 * Generate a grid of test points based on the raytracer's bbox,
 * @param raytracer
 * @param Nxyz
 * @param points
 */
void grid_test_points_2d(const RayTracer2D &raytracer_2d, std::vector<int> Nxy, std::vector<test_point_2d> &points) {
  points.clear();
  assert(Nxy.size() == 2);
  std::vector<double> center(3, 0), lb(3, 0), ub(3, 0);
  for (int dim = 0; dim < 2; dim++) {
    center[dim] = (raytracer_2d.bounds().min[dim] + raytracer_2d.bounds().max[dim]) / 2;
    lb[dim] = center[dim] - (center[dim] - raytracer_2d.bounds().min[dim]) * 1.2;
    ub[dim] = center[dim] + (raytracer_2d.bounds().max[dim] - center[dim]) * 1.2;
  }
  points.reserve(Nxy[0] * Nxy[1]);
  for (int i = 0; i < Nxy[0]; i++) {
    for (int j = 0; j < Nxy[1]; j++) {
//        test_point temp(i + double(rand()) / RAND_MAX, j + double(rand()) / RAND_MAX, k + double(rand()) / RAND_MAX);
      test_point_2d temp(static_cast<float>(lb[0] + i * (ub[0] - lb[0]) / Nxy[0]),
                         static_cast<float>(lb[1] + j * (ub[1] - lb[1]) / Nxy[1]));
      points.emplace_back(temp);
    }
  }
}

void Vec_Sub(double *dest, const double *v1, const double *v2) {
  dest[0] = v1[0] - v2[0];
  dest[1] = v1[1] - v2[1];
  dest[2] = v1[2] - v2[2];
}

void Vec_Sub_2D(double *dest, const double *v1, const double *v2) {
  dest[0] = v1[0] - v2[0];
  dest[1] = v1[1] - v2[1];
}

void Vec_Cross(double *dest, const double *v1, const double *v2) {
  dest[0] = v1[1] * v2[2] - v1[2] * v2[1];
  dest[1] = v1[2] * v2[0] - v1[0] * v2[2];
  dest[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

double Vec_Dot(const double *v1, const double *v2) {
  return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}

double Vec_Dot_2D(const double *v1, const double *v2) {
  return (v1[0] * v2[0] + v1[1] * v2[1]);
}

void boundary_test_points(const RayTracer &raytracer, std::vector<test_point> &points, int no_of_points) {
  points.clear();
  points.reserve(no_of_points + 10);
  unsigned long int step = raytracer.getNumTriangles() / no_of_points;
  if (step == 0) step = 1;
  for (unsigned long i = 0; i < raytracer.getNumTriangles(); i += step) {
    const auto &tri = raytracer.triangles()[i];
    double xc = (tri.v[0].data[0] + tri.v[1].data[0] + tri.v[2].data[0]) / 3;
    double yc = (tri.v[0].data[1] + tri.v[1].data[1] + tri.v[2].data[1]) / 3;
    double zc = (tri.v[0].data[2] + tri.v[1].data[2] + tri.v[2].data[2]) / 3;
    double edge1[3], edge2[3], normal[3];
    Vec_Sub(edge1, tri.v[1].data, tri.v[0].data);
    Vec_Sub(edge2, tri.v[2].data, tri.v[0].data);
    Vec_Cross(normal, edge1, edge2);
    double normal_mag = sqrt(Vec_Dot(normal, normal));
    normal[0] /= normal_mag;
    normal[1] /= normal_mag;
    normal[2] /= normal_mag;

    double shift_x = 0.001 * normal[0]; //double(rand()) / RAND_MAX - 0.5;
    double shift_y = 0.001 * normal[1]; //double(rand()) / RAND_MAX - 0.5;
    double shift_z = 0.001 * normal[2]; //double(rand()) / RAND_MAX - 0.5;
    test_point temp(static_cast<float>(xc + shift_x), static_cast<float>(yc + shift_y), static_cast<float>(zc + shift_z));
    points.emplace_back(temp);
  }
}

void boundary_test_points_2d(const RayTracer2D &raytracer_2d, std::vector<test_point_2d> &points, int no_of_points) {
  points.clear();
  points.reserve(no_of_points + 10);
  unsigned long int step = raytracer_2d.getNumLines() / no_of_points;
  if (step == 0) step = 1;
  for (unsigned long i = 0; i < raytracer_2d.getNumLines(); i += step) {
    const auto &line = raytracer_2d.lines()[i];
    double xc = (line.v[0].data[0] + line.v[1].data[0]) / 2;
    double yc = (line.v[0].data[1] + line.v[1].data[1]) / 2;
    double edge1[2], normal[2];
    Vec_Sub_2D(edge1, line.v[1].data, line.v[0].data);
    normal[0] = edge1[1];
    normal[0] = -edge1[0];
    double normal_mag = sqrt(Vec_Dot_2D(normal, normal));
    normal[0] /= normal_mag;
    normal[1] /= normal_mag;

    double shift_x = 0.001 * normal[0]; //double(rand()) / RAND_MAX - 0.5;
    double shift_y = 0.001 * normal[1]; //double(rand()) / RAND_MAX - 0.5;
    test_point_2d temp(static_cast<float>(xc + shift_x), static_cast<float>(yc + shift_y));
    points.emplace_back(temp);
  }
}
void read_test_points(const std::string &filename, std::vector<test_point> &points) {
  points.clear();
  std::ifstream myFile(filename);
  // Make sure the file is open
  if (!myFile.is_open()) throw std::runtime_error("Could not open file");


  // Helper vars
  std::string line, colname;
  // Read the column names
  if (myFile.good()) {
    // Extract the first line in the file
    std::getline(myFile, line);
  }

  double val = 0.0;

  // Read data, line by line
  while (std::getline(myFile, line)) {
    // Create a stringstream of the current line
    std::stringstream ss(line);

    std::vector<double> coor;
    // Extract each integer
    while (ss >> val) {
      // Add the current integer to the 'colIdx' column's values vector
      coor.push_back(val);
      // If the next token is a comma, ignore it and move on
      if (ss.peek() == ',' or ss.peek() == ' ') {
        ss.ignore();
      }
    }
    test_point temp(coor[0], coor[1], coor[2]);
    points.push_back(temp);
  }

  // Close file
  myFile.close();
}

void output_points(const std::string &name, const std::vector<test_point> &points) {
  std::ofstream pointFile;
  pointFile.open((name + ".3D").c_str(), std::ios::out);
  pointFile << "X Y Z value\n";
  for (const auto &p: points) {
    pointFile.precision(14);
    pointFile << p.xyz[0] << " "
              << p.xyz[1] << " "
              << p.xyz[2] << " "
              << p.flag << std::endl;
  }
  pointFile.close();
}

void output_points_2d(const std::string &name, const std::vector<test_point_2d> &points) {
  std::ofstream pointFile;
  pointFile.open((name + ".3D").c_str(), std::ios::out);
  pointFile << "X Y Z value\n";
  for (const auto &p: points) {
    pointFile.precision(14);
    pointFile << p.xy[0] << " "
              << p.xy[1] << " "
              << 0 << " "
              << p.flag << std::endl;
  }
  pointFile.close();
}

void test_ray_tracing(RayTracer &rayTracer,
                      std::vector<test_point> &points,
                      bool out_file = false,
                      const std::string &pointFilename = "") {
  /// Testing points
  std::cout << "Testing: " << points.size() << " points." << std::endl;
  /// Testing all the points
  high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < points.size(); i++) {
    points[i].flag = rayTracer.ifInside(points[i].xyz);
    if (points[i].flag == InsideType::UNKNOWN) {
      std::cout << "Could not test point " << i << "\n";
      assert(false);
    }
    /// Output progress
    if (i % (points.size() / 10) == 0) {
      std::cout << "Points tested: " << i << "\tPoints remaining: " << points.size() - i << std::endl;
    }
  }
  high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  duration<double> ray_trace_time = duration_cast<duration<double>>(t2 - t1);
  std::cout << "time_span_ray_trace: " << ray_trace_time.count() << "\n";
#if defined(BENCHMARK)
  #if defined(BENCHMARKLOG)
  rayTracer.printTimer(pointFilename);
#else
  rayTracer.printTimer();
#endif
  rayTracer.cleanUpTimer();
#endif
  if (out_file) {
    output_points(pointFilename, points);
  }
}

void test_ray_tracing_2d(RayTracer2D &rayTracer_2d,
                         std::vector<test_point_2d> &points,
                         bool out_file = false,
                         const std::string &pointFilename = "") {
  /// Testing points
  std::cout << "Testing: " << points.size() << " points." << std::endl;
  /// Testing all the points
  high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < points.size(); i++) {
    points[i].flag = rayTracer_2d.ifInside(points[i].xy);
    if (points[i].flag == InsideType2D::UNKNOWN) {
      std::cout << "Could not test point " << i << "\n";
      assert(false);
    }
    /// Output progress
    if (i % (points.size() / 10) == 0) {
      std::cout << "Points tested: " << i << "\tPoints remaining: " << points.size() - i << std::endl;
    }
  }
  high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  duration<double> ray_trace_time = duration_cast<duration<double>>(t2 - t1);
  std::cout << "time_span_ray_trace: " << ray_trace_time.count() << "\n";
#if defined(BENCHMARK)
  #if defined(BENCHMARKLOG)
  rayTracer_2d.printTimer(pointFilename);
#else
  rayTracer_2d.printTimer();
#endif
  rayTracer_2d.cleanUpTimer();
#endif
  if (out_file) {
    output_points_2d(pointFilename, points);
  }
}

#endif //RAY_TRACING_TESTPOINTS_H
