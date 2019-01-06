/*
 * isopointal.cpp
 * Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "isopointal.h"

#include <sstream>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "util.h"
#include "wallpaper.h"

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
