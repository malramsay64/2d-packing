/*
 * wallpaper.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "wallpaper.h"

#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "basis.h"
#include "geometry.h"
#include "math.h"
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

std::string IsopointalGroup::group_string() const {
  std::stringstream ss;
  for (const auto& site : this->wyckoff_sites) {
    ss << site.letter;
  }
  return ss.str();
}

std::size_t IsopointalGroup::group_multiplicity() const {
  std::size_t multiplicity{0};
  for (const WyckoffSite& site : this->wyckoff_sites) {
    multiplicity += site.multiplicity();
  }
  return multiplicity;
}

std::vector<IsopointalGroup> generate_isopointal_groups(
    const Shape& shape,
    const WallpaperGroup& group,
    std::size_t num_occupied_sites) {
  auto console = spdlog::stdout_color_mt("console");

  std::vector<WyckoffSite> valid_sites;
  for (const WyckoffSite& wyckoff : group.wyckoff_sites) {
    // first test if the shape has the required symmetries for various sites
    // at the moment only tests rotations & mirrors... that's all?
    int rotational_match{shape.rotational_symmetries % wyckoff.rotations};
    int mirror_match{shape.mirrors % wyckoff.mirrors};
    if (rotational_match == 0) {
      if ((wyckoff.mirrors == 0) || ((mirror_match == 0) && (shape.mirrors != 0))) {
        if (wyckoff.variability) {
          console->debug("site %c is valid and variable\n", wyckoff.letter);
        } else {
          printf("site %c is valid\n", wyckoff.letter);
        }
        if (wyckoff.variability) {
          // Variable sites are added with replacement, so add one for each occupied
          // site
          valid_sites.insert(valid_sites.end(), num_occupied_sites, wyckoff);
        } else {
          valid_sites.push_back(wyckoff);
        }
      } else {
        console->debug(
            "site %c does not have the required rotational symmetry\n", wyckoff.letter);
      }
    }
  }

  auto occupied_sites = combinations<WyckoffSite>(valid_sites, num_occupied_sites);
  uniqueify<std::vector<WyckoffSite>>(occupied_sites);

  console->info(
      "enumeration complete: in fact there were %lu ways\n", occupied_sites.size());

  std::vector<IsopointalGroup> isopointal_groups;
  for (auto& combination : occupied_sites) {
    isopointal_groups.push_back(IsopointalGroup(combination));
    console->debug("%s\n", isopointal_groups.back().group_string());
  }
  return isopointal_groups;
}
