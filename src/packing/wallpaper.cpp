/*
 * wallpaper.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "wallpaper.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "basis.h"
#include "geometry.h"
#include "math.h"
#include "random.h"
#include "shapes.h"
#include "util.h"

bool SymmetryTransform::operator==(const SymmetryTransform& other) const {
  return (
      this->x_coeffs == other.x_coeffs && this->y_coeffs == other.y_coeffs &&
      this->rotation_offset == other.rotation_offset &&
      this->site_mirror == other.site_mirror);
}

Vect2 SymmetryTransform::real_to_fractional(const Vect3& real) const {
  /* converts site variables and wyckoff site coefficients into the location
   * of the actual wyckoff image in fractional coordinates */
  Vect2 v(0, 0);
  v.x = this->x_coeffs.x * real.x + this->x_coeffs.y * real.y + this->x_coeffs.z;
  v.y = this->y_coeffs.x * real.x + this->y_coeffs.y * real.y + this->y_coeffs.z;
  return positive_modulo(v, 1.);
}

bool WyckoffSite::operator==(const WyckoffSite& other) const {
  return (
      this->letter == other.letter && this->variability == other.variability &&
      this->rotations == other.rotations && this->mirrors == other.mirrors &&
      this->symmetries == other.symmetries);
}

std::size_t WyckoffSite::multiplicity() const {
  return this->symmetries.size();
}

bool WyckoffSite::vary_x() const {
  // If the first symmetry can vary x, all of them can
  return fabs(this->symmetries[0].x_coeffs.x) > 0.1;
}

bool WyckoffSite::vary_y() const {
  // If the first symmetry can vary x, all of them can
  return fabs(this->symmetries[0].y_coeffs.y) > 0.1;
}

int WyckoffSite::mirror_type() const {
  return this->symmetries[0].site_mirror;
}

/** Check whether two shape instances intersect
 */
bool check_for_intersection(
    const ShapeInstance& shape_a,
    const ShapeInstance& shape_b,
    const Cell& cell) {

  // a is fixed, b is moved to the periodic sites to test for the intersection
  Vect2 fcoords_b{shape_b.get_fractional_coordinates()};
  Vect2 img_fcoords_b{0, 0};

  int shells = 1;
  // For extreme angles, using only the nearest shell fails, so have to look at 2
  // shells. The designator for 'extreme' angle is PI/4 or 45 degrees.
  if (cell.angle->get_value() < M_PI_4) {
    shells = 2;
  } else if (M_2_PI - cell.angle->get_value() < M_PI_4) {
    shells = 2;
  }

  // Loop over the possible periodic positions
  for (int cell_img_x = -shells; cell_img_x <= shells; cell_img_x++) {
    for (int cell_img_y = -shells; cell_img_y <= shells; cell_img_y++) {
      // Intersections with one's self are excluded
      if ((shape_a == shape_b) && (cell_img_x == 0) && (cell_img_y == 0)) {
        continue;
      }
      Vect2 coords_b = cell.fractional_to_real(
          Vect2(fcoords_b.x + cell_img_x, fcoords_b.y + cell_img_y));
      if (shape_a.intersects_with(shape_b, coords_b)) {
        return true;
      }
    }
  }
  return false;
}

bool check_state_for_intersection(
    const Shape& shape,
    const std::vector<OccupiedSite>& occupied_sites,
    const Cell& cell) {
  // Loop over all the occupied sites
  for (auto site_one = occupied_sites.begin(); site_one != occupied_sites.end();
       site_one++) {
    // Loop over all symmetries for the first occupied site
    for (const auto& image_one : site_one->wyckoff->symmetries) {
      const ShapeInstance shape_one{shape, *site_one, image_one};
      // Loop over all occupied sites which haven't already been compared with site_one
      for (auto site_two = std::next(site_one); site_two != occupied_sites.end();
           site_two++) {
        // Loop over all symmetries for the second occupied site
        for (const auto& image_two : site_two->wyckoff->symmetries) {
          const ShapeInstance shape_two{shape, *site_two, image_two};
          /* Finally perform the comparison of shapes here */
          if (check_for_intersection(shape_one, shape_two, cell)) {
            // If the two shapes intersect, return true, breaking out of the loop.
            return true;
          }
        }
      }
    }
  }
  // Should there be no intersections between any shapes, return false
  return false;
}

double calculate_packing_fraction(
    const Shape& shape,
    const Cell& cell,
    const std::vector<OccupiedSite>& occupied_sites) {
  std::size_t shape_replicas(calculate_shape_replicas(occupied_sites));

  double packing_fraction = shape_replicas * shape.area() / cell.area();
  if (std::isnan(packing_fraction)) {
    auto console = spdlog::stdout_color_mt("console");
    console->warn(
        "nan encountered %f %f %f\n",
        cell.x_len->get_value(),
        cell.y_len->get_value(),
        cell.angle->get_value());
    throw "Packing fraction calculation failed";
  }

  return packing_fraction;
}

bool ShapeInstance::operator==(const ShapeInstance& other) const {
  return (
      this->shape == other.shape && this->site == other.site &&
      this->symmetry_transform == other.symmetry_transform);
}

Vect2 ShapeInstance::get_fractional_coordinates() const {
  return this->symmetry_transform->real_to_fractional(this->site->site_variables());
}

Vect2 ShapeInstance::get_real_coordinates() const {
  return this->site->get_position();
}

double ShapeInstance::get_angle() const {
  return this->site->angle->get_value();
}

double ShapeInstance::get_rotational_offset() const {
  return this->symmetry_transform->rotation_offset;
}

bool ShapeInstance::get_flipped() const {
  return this->symmetry_transform->flipped ^ this->site->flip_site;
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
