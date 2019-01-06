/*
 * shapes.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "shapes.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <pybind11/stl.h>

#include "math.h"

namespace py = pybind11;

Shape::Shape(
    const std::string& name,
    const std::vector<double>& radial_points,
    const int rotational_symmetries,
    const int mirrors)
    : name(name), radial_points(radial_points),
      rotational_symmetries(rotational_symmetries), mirrors(mirrors) {
  auto min_max = std::minmax_element(radial_points.begin(), radial_points.end());
  this->min_radius = *min_max.first;
  this->max_radius = *min_max.second;
}

Shape::Shape(const std::string& name, const std::vector<double>& radial_points)
    : Shape(name, radial_points, 1, 0){};

int Shape::resolution() const {
  return this->radial_points.size();
}

double Shape::angular_step() const {
  return 2 * PI / this->resolution();
}

double Shape::get_point(const int index) const {
  return this->radial_points.at(index);
}

void Shape::plot(const std::string& filename) const {
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);

  for (std::size_t i = 0; i < this->radial_points.size(); ++i) {
    double angle{this->angular_step() * i};
    outfile << std::fixed << std::setprecision(12)
            << this->radial_points[i] * cos(angle) << " "
            << this->radial_points[i] * sin(angle) << std::endl;
  }
  outfile.close();
}

/* This uses the side-angle-side method of calculating the area, that is
 *     Area = side_a * side_b * sin(angle_ab)
 * In this case the angle is the same for every triangle so is precomputed above.
 */
double Shape::area() const {
  double areasum{0.0};
  const double angle{std::sin(this->angular_step())};
  for (std::size_t index = 0; index < this->radial_points.size(); ++index) {
    std::size_t next_index = (index + 1) % this->radial_points.size();
    areasum +=
        0.5 * this->radial_points[index] * this->radial_points[next_index] * angle;
  }
  return areasum;
}

std::vector<Vect2> Shape::generate_position_cache(
    const Vect2& position,
    const double angle_to_shape) const {
  const int resolution{this->resolution()};

  std::vector<Vect2> position_cache;
  // Reserve the expected size of the vector on initialisation
  position_cache.reserve(resolution / 2 + 1);

  const double angular_step = M_2_PI / resolution;
  int angle_int = static_cast<int>(std::round(angle_to_shape / angular_step));

  for (int index = -resolution / 4; index <= (resolution / 4); index++) {
    int compare_index{angle_int + index};
    // Change base to angle between shapes
    double theta{fabs(compare_index * angular_step - angle_to_shape)};
    position_cache.push_back(Vect2(
        this->get_point(compare_index) * cos(theta),
        this->get_point(compare_index) * sin(theta)));
  }
  return position_cache;
}

std::vector<Vect2> Shape::generate_position_cache_full(const Vect2& position) const {
  const int resolution{this->resolution()};

  std::vector<Vect2> position_cache;
  // Reserve the expected size of the vector on initialisation
  position_cache.reserve(resolution / 2 + 1);

  const double angular_step{M_2_PI / resolution};

  for (int index = -resolution / 2; index <= (resolution / 2); index++) {
    // Change base to angle between shapes
    double theta{std::fabs(index * angular_step)};
    position_cache.push_back(Vect2(
        this->get_point(index) * std::cos(theta),
        this->get_point(index) * std::sin(theta)));
  }
  return position_cache;
}

// Export the shape class to python uisng pybind11
void export_Shape(py::module& m) {
  py::class_<Shape> shape(m, "Shape");
  shape
      .def(
          py::init<const std::string&, const std::vector<double>&, int, int>(),
          py::arg("name"),
          py::arg("radial_points"),
          py::arg("rotational_symmetries") = 0,
          py::arg("mirrors") = 0)
      .def("resolution", &Shape::resolution)
      .def("plot", &Shape::plot)
      .def("angular_step", &Shape::angular_step)
      .def("get_point", &Shape::get_point)
      .def("area", &Shape::area)
      .def_readonly("name", &Shape::name)
      .def_readonly("radial_points", &Shape::radial_points)
      .def_readonly("rotational_symmetries", &Shape::rotational_symmetries)
      .def_readonly("mirrors", &Shape::mirrors);
}
