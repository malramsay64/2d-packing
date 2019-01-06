/*
 * packing.h
 * Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <memory>
#include <vector>

#include "basis.h"
#include "shapes.h"
#include "wallpaper.h"

#ifndef PACKING_H
#define PACKING_H

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

std::vector<OccupiedSite> initialise_structure(
    const Shape& shape,
    const IsopointalGroup& isopointal_group,
    const WallpaperGroup& group,
    Cell& cell,
    std::vector<Basis>& basis,
    const double step_size);

double calculate_packing_fraction(
    const Shape& shape,
    const Cell& cell,
    const std::vector<OccupiedSite>& occupied_sites);

bool check_for_intersection(
    const ShapeInstance& shape_a,
    const ShapeInstance& shape_b,
    const Cell& cell);

bool check_state_for_intersection(
    const Shape& shape,
    const std::vector<OccupiedSite>& occupied_sites,
    const Cell& cell);

std::size_t calculate_shape_replicas(const std::vector<OccupiedSite>& sites);

#endif /* !PACKING_H */
