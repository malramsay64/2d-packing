/*
 * shapes.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "shapes.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <pybind11/stl.h>

#include "geometry.h"
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
  // Flipping reverses the direction of the points

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
  // Flipping reverses the direction of the points

  for (int index = -resolution / 2; index <= (resolution / 2); index++) {
    // Change base to angle between shapes
    double theta{std::fabs(index * angular_step)};
    position_cache.push_back(Vect2(
        this->get_point(index) * std::cos(theta),
        this->get_point(index) * std::sin(theta)));
  }
  return position_cache;
}

bool ImageType::operator==(const ImageType& other) const {
  return (
      this->x_coeffs == other.x_coeffs && this->y_coeffs == other.y_coeffs &&
      this->rotation_offset == other.rotation_offset &&
      this->site_mirror == other.site_mirror);
}

Vect2 ImageType::real_to_fractional(const Vect3& real) const {
  /* converts site variables and wyckoff site coefficients into the location
   * of the actual wyckoff image in fractional coordinates */
  Vect2 v(0, 0);
  v.x = this->x_coeffs.x * real.x + this->x_coeffs.y * real.y + this->x_coeffs.z;
  v.y = this->y_coeffs.x * real.x + this->y_coeffs.y * real.y + this->y_coeffs.z;
  return positive_modulo(v, 1.);
}

Vect2 ImageType::real_to_fractional(const Site& site) const {
  /* converts site variables and wyckoff site coefficients into the location
   * of the actual wyckoff image in fractional coordinates */
  Vect2 v(0, 0);
  v.x = this->x_coeffs.x * site.x->get_value() +
        this->x_coeffs.y * site.y->get_value() + this->x_coeffs.z;
  v.y = this->y_coeffs.x * site.x->get_value() +
        this->y_coeffs.y * site.y->get_value() + this->y_coeffs.z;
  return positive_modulo(v, 1.);
}

bool WyckoffType::operator==(const WyckoffType& other) const {
  return (
      this->multiplicity == other.multiplicity && this->letter == other.letter &&
      this->some_variability == other.some_variability &&
      this->site_rotations == other.site_rotations &&
      this->site_mirrors == other.site_mirrors && this->image == other.image);
}

bool ShapeInstance::operator==(const ShapeInstance& other) const {
  return (
      this->shape == other.shape && this->site == other.site &&
      this->image == other.image);
}

Vect2 ShapeInstance::get_fractional_coordinates() const {
  return this->image->real_to_fractional(*this->site);
}

Vect2 ShapeInstance::get_real_coordinates() const {
  return this->site->get_position();
}

double ShapeInstance::get_angle() const {
  return this->site->angle->get_value();
}

double ShapeInstance::get_rotational_offset() const {
  return this->image->rotation_offset;
}

bool ShapeInstance::get_flipped() const {
  return this->image->flipped ^ this->site->flip_site;
}

std::pair<double, double> ShapeInstance::compute_incline(
    const ShapeInstance& other,
    const Vect2& position_other) const {

  const Vect2& position_this{this->get_real_coordinates()};
  const double central_dist{(position_this - position_other).norm()};
  double a_to_b_incline{acos((position_other.x - position_this.x) / central_dist)};

  if (std::isnan(a_to_b_incline)) {
    if (is_close(position_other.x - position_this.x, central_dist, 1e-8)) {
      a_to_b_incline = 0.0;
    } else if (is_close(position_other.x - position_this.x, -central_dist, 1e-8)) {
      a_to_b_incline = M_PI;
    }
  }
  if (position_other.x < position_this.x) {
    a_to_b_incline = M_2_PI - a_to_b_incline;
  }

  // Set reverse incline
  double b_to_a_incline{a_to_b_incline + M_PI};

  // Deal with the flipping of shapes
  if (this->get_flipped()) {
    a_to_b_incline = M_2_PI - a_to_b_incline;
  }
  if (other.get_flipped()) {
    b_to_a_incline = M_2_PI - b_to_a_incline;
  }

  /* now add in the rotation due to the orientation parameters */
  /* This is unaffected by the flip states */
  a_to_b_incline += this->get_angle();
  b_to_a_incline += other.get_angle();

  /* now add in the rotation due to the rotation of this image wrt the other
   * images of the same wyckoff */
  a_to_b_incline += std::pow(-1, this->get_flipped()) * this->get_rotational_offset();
  b_to_a_incline += std::pow(-1, other.get_flipped()) * other.get_rotational_offset();

  a_to_b_incline = positive_modulo(a_to_b_incline, M_2_PI);
  b_to_a_incline = positive_modulo(b_to_a_incline, M_2_PI);
  return std::pair<double, double>{a_to_b_incline, b_to_a_incline};
}

bool ShapeInstance::intersects_with(
    const ShapeInstance& other,
    const Vect2& position_other) const {

  const Vect2& position_this{this->get_real_coordinates()};
  const double central_dist{(position_this - position_other).norm()};
  /* No clash when further apart than the maximum shape radii measures */
  if (central_dist > this->shape->max_radius + other.shape->max_radius) {
    return false;
  }

  double angle_this_to_other, angle_other_to_this;
  std::tie(angle_this_to_other, angle_other_to_this) =
      this->compute_incline(other, position_other);

  std::vector<Vect2> position_a_cache =
      this->shape->generate_position_cache(position_this, angle_this_to_other);
  std::vector<Vect2> position_b_cache =
      other.shape->generate_position_cache(position_other, angle_other_to_this);

  // The final element is the inital previous position providing a closed shape
  // regardless of the number of points checked.
  Vect2& position_a_prev = position_a_cache.back();
  Vect2& position_b_prev = position_b_cache.back();
  for (auto position_a = position_a_cache.begin(); position_a != position_a_cache.end();
       ++position_a) {
    for (auto position_b = position_b_cache.begin();
         position_b != position_b_cache.end();
         ++position_b) {
      if (segments_cross(position_a_prev, *position_a, position_b_prev, *position_b)) {
        return true;
      }
      // Update the previous point
      position_a_prev = *position_a;
      position_b_prev = *position_b;
    }
  }
  return false;
}

std::size_t group_multiplicity(const std::vector<Site>& occupied_sites) {
  std::size_t total_multiplicity = 0;
  for (const Site& site : occupied_sites) {
    total_multiplicity = site.wyckoff->multiplicity;
  }
  return total_multiplicity;
}

std::string create_filename(
    const Shape& shape,
    const WallpaperGroup& group,
    const std::string site_list,
    const std::string& directory) {
  std::stringstream stream_filename;
  stream_filename << directory << "/solution_" << shape.name << "_" << group.label
                  << "_" << site_list << "_" << std::fixed << std::setprecision(3)
                  << shape.shape_var << ".svg";
  return stream_filename.str();
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
      .def("area", &Shape::area);
}
