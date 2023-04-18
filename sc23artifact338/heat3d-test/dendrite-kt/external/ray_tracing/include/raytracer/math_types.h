//
// Created by lofquist on 1/18/18.
//

#ifndef RAY_TRACING_MATH_H
#define RAY_TRACING_MATH_H

#include <cmath>

struct Vector3d {
  double data[3];
  Vector3d() : data{} {}
  Vector3d(float xp, float yp, float zp) : data{xp, yp, zp} {}
};

struct Vector2d {
  double data[2];
  Vector2d() : data{} {}
  Vector2d(float xp, float yp) : data{xp, yp} {}
};


struct Box3d {
  double min[3];
  double max[3];
  Box3d() : min{}, max{} {};
  Box3d(float xmin, float ymin, float zmin, float xmax, float ymax, float zmax)
      : min{xmin, ymin, zmin}, max{xmax, ymax, zmax} {}

  /**
   * Returns true if this and other are overlapping.
   * @param other
   * @return
   */
  inline bool intersects(const Box3d& other) const {
      return !((this->min[0] > other.max[0] || this->min[1] > other.max[1] || this->min[2] > other.max[2])
          || (other.min[0] > this->max[0] || other.min[1] > this->max[1] || other.min[2] > this->max[2]));
  }
};

struct Box2d {
  double min[2];
  double max[2];
  Box2d() : min{}, max{} {};
  Box2d(float xmin, float ymin, float xmax, float ymax)
      : min{xmin, ymin}, max{xmax, ymax} {}

  /**
   * Returns true if this and other are overlapping.
   * @param other
   * @return
   */
  inline bool intersects(const Box2d& other) const {
    return !((this->min[0] > other.max[0] || this->min[1] > other.max[1])
        || (other.min[0] > this->max[0] || other.min[1] > this->max[1]));
  }
};

struct Sphere {
  double center[3];
  double radius;
  Sphere() : center{}, radius {} {};
  Sphere(float xc, float yc, float zc, float radius)
      : center{xc, yc, zc}, radius{radius} {}

  static Sphere fromTriangle(const Vector3d* v) {
    double xc = (v[0].data[0] + v[1].data[0] + v[2].data[0]) / 3;
    double yc = (v[0].data[1] + v[1].data[1] + v[2].data[1]) / 3;
    double zc = (v[0].data[2] + v[1].data[2] + v[2].data[2]) / 3;
    double distance_1 = std::sqrt(std::pow(v[0].data[0] - xc, 2) + std::pow(v[0].data[1] - yc, 2) + std::pow(v[0].data[2] - zc, 2));
    double distance_2 = std::sqrt(std::pow(v[1].data[0] - xc, 2) + std::pow(v[1].data[1] - yc, 2) + std::pow(v[1].data[2] - zc, 2));
    double distance_3 = std::sqrt(std::pow(v[2].data[0] - xc, 2) + std::pow(v[2].data[1] - yc, 2) + std::pow(v[2].data[2] - zc, 2));
    double radius = distance_1;
    if (radius < distance_2) {
      radius = distance_2;
    }
    if (radius < distance_3) {
      radius = distance_3;
    }
    double eps = 1e-2;
    return {static_cast<float>(xc), static_cast<float>(yc), static_cast<float>(zc), static_cast<float>(radius * (1 + eps))};
  }
};

struct Circle {
  double center[2];
  double radius;
  Circle() : center{}, radius {} {};
  Circle(float xc, float yc, float radius)
      : center{xc, yc}, radius{radius} {}

  static Circle fromLine(const Vector2d* v) {
    double xc = (v[0].data[0] + v[1].data[0]) / 2;
    double yc = (v[0].data[1] + v[1].data[1]) / 2;
    double distance_1 = std::sqrt(std::pow(v[0].data[0] - xc, 2) + std::pow(v[0].data[1] - yc, 2));
    double distance_2 = std::sqrt(std::pow(v[1].data[0] - xc, 2) + std::pow(v[1].data[1] - yc, 2));
    double radius = distance_1;
    if (radius < distance_2) {
      radius = distance_2;
    }
    double eps = 1e-2;
    return {static_cast<float>(xc), static_cast<float>(yc), static_cast<float>(radius * (1 + eps))};
  }
};


#endif //RAY_TRACING_MATH_H
