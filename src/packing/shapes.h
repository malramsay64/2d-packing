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

/*! \enum mirror
 *
 *  The type of mirror symmetry a SymmetryTransform can have.
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

/** \class SymmetryTransform
 *
 * Define the transformations of particle positions for a WyckoffSite.
 *
 * A collection of operators which define the symmetry transformations.
 *
 */
class SymmetryTransform {
public:
  SymmetryTransform(Vect3 x_coeffs, Vect3 y_coeffs, double rotation_offset);
  Vect3 x_coeffs;
  Vect3 y_coeffs;
  double rotation_offset;
  bool flipped;
  enum mirror site_mirror;

  bool operator==(const SymmetryTransform& other) const;
  friend std::ostream& operator<<(std::ostream&, const SymmetryTransform&);

  Vect2 real_to_fractional(const Vect3& real) const;
};

/** \class WyckoffSite
 *
 * A class which defines a set of Wyckoff parameters.
 *
 * Each set of Wyckoff paramters defines a set of properties, including the positions of
 * particles in Wyckoff Site. A Wyckoff site can have multiple positions, indicated by
 * the muiltiplicity. Each of these positions have some symmetry relationship to each
 * other, either rotaional or mirror.
 *
 * Each Wyckoff instance is labelled by a letter which is an identifier from the
 * Crystallographic Database.
 *
 * A WyckoffSite instance only describes the transformations for taking a position and
 * obtaining the symmetry related positions.
 *
 */
class WyckoffSite {
public:
  char letter;
  int variability;
  int rotations;
  int mirrors;
  std::vector<SymmetryTransform> symmetries;

  std::size_t multiplicity() const;
  bool vary_x() const;
  bool vary_y() const;
  int mirror_type() const;
  std::string str() const;
};

/** \class ShapeInstance
 *
 * A specific instance of a Shape object which has coordinates and orientation.
 *
 * To define the boundary of a shape in cartesian coordinates, a lot of information is
 * required.
 *  - The actual Shape class itself, which radially defines the positions of
 *  particles.
 *  - The OccupiedSite, which defines the WyckoffSite and the coordinates of the site
 *  - The SymmetryTransform, defining which of the symmetry transforms of the
 *  WyckoffSite this particular shape occupies.
 *
 */
class ShapeInstance {
  const std::shared_ptr<const Shape> shape;
  const std::shared_ptr<const OccupiedSite> site;
  const std::shared_ptr<const SymmetryTransform> symmetry_transform;

public:
  ShapeInstance(
      std::shared_ptr<const Shape> shape,
      std::shared_ptr<const OccupiedSite> site,
      std::shared_ptr<const SymmetryTransform> symmetry_transform)
      : shape(shape), site(site), symmetry_transform(symmetry_transform){};
  ShapeInstance(
      const Shape& shape,
      const OccupiedSite& site,
      const SymmetryTransform& symmetry_transform)
      : ShapeInstance(
            std::shared_ptr<const Shape>(&shape),
            std::shared_ptr<const OccupiedSite>(&site),
            std::shared_ptr<const SymmetryTransform>(&symmetry_transform)){};

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

void export_Shape(pybind11::module& m);

#endif /* SHAPES_H */
