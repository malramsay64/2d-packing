/*
 * shapes.h
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>

#include "basis.h"
#include "math.h"

#ifndef SHAPES_H
#define SHAPES_H

/* \class Shape
 *
 * Defines a shape from a set of radially defined points.
 *
 */
class Shape {
public:
  Shape(
      const std::string& name,
      const std::vector<double>& radial_points,
      const int rotational_symmetries,
      const int mirrors);
  Shape(const std::string& name, const std::vector<double>& radial_points);

  std::string name;
  std::vector<double> radial_points;
  int rotational_symmetries;
  int mirrors;
  double min_radius;
  double max_radius;
  double shape_var = 0;

  int resolution() const;
  double angular_step() const;
  double get_point(int index) const;

  void plot(const std::string& filename) const;
  double area() const;

  std::vector<Vect2>
  generate_position_cache(const Vect2& position, double angle_to_shape) const;
  std::vector<Vect2> generate_position_cache_full(const Vect2& position) const;
};

void export_Shape(pybind11::module& m);

#endif /* SHAPES_H */
