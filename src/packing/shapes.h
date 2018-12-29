/*
 * shapes.h
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <string>
#include <vector>

#include <pybind11/pybind11.h>

#include "basis.h"
#include "math.h"

#ifndef SHAPES_H
#define SHAPES_H

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

/*! \enum mirror
 *
 *  Detailed description
 */
enum mirror {
  m0 = 0,
  m90 = 90,
  m45 = 45,
  m135 = 135,
  m30 = 30,
  m60 = 60,
  m330 = 330,
  m300 = 300
};

class ImageType {
public:
  ImageType(Vect3 x_coeffs, Vect3 y_coeffs, double rotation_offset);
  Vect3 x_coeffs;
  Vect3 y_coeffs;
  double rotation_offset;
  bool flipped;
  enum mirror site_mirror;

  bool operator==(const ImageType& other) const;

  Vect2 real_to_fractional(const Vect3& real) const;
  Vect2 real_to_fractional(const Site& site) const;
};

class WyckoffType {
public:
  int multiplicity;
  char letter;
  int some_variability;
  int site_rotations;
  int site_mirrors;
  std::vector<ImageType> image;

  bool operator==(const WyckoffType& other) const;
};

class ShapeInstance {
public:
  ShapeInstance(const Shape& shape, Site& site, ImageType& image)
      : shape(&shape), site(&site), image(&image){};
  const std::shared_ptr<const Shape> shape;
  const std::shared_ptr<const Site> site;
  const std::shared_ptr<const ImageType> image;

  bool operator==(const ShapeInstance& other) const;

  Vect2 get_fractional_coordinates() const;
  Vect2 get_real_coordinates() const;
  double get_angle() const;
  double get_rotational_offset() const;
  bool get_flipped() const;
  bool intersects_with(const ShapeInstance& other, const Vect2& coords_other) const;
  std::pair<double, double>
  compute_incline(const ShapeInstance& other, const Vect2& position_other) const;
};

class WallpaperGroup {
public:
  WallpaperGroup(std::string label, std::vector<WyckoffType> wyckoffs);
  const std::string label;
  const bool a_b_equal = false;
  const bool hexagonal = false;
  const bool rectangular = false;
  int num_symmetries = 0;
  std::vector<WyckoffType> wyckoffs;
  int num_wyckoffs;
};

std::size_t group_multiplicity(const std::vector<Site>& occupied_sites);

std::string create_filename(
    const Shape& shape,
    const WallpaperGroup& group,
    const std::string site_list,
    const std::string& directory);

void export_Shape(pybind11::module& m);

#endif /* SHAPES_H */
