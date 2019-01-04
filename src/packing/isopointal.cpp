/*
 * isopointal.cpp
 * Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "isopointal.h"

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "util.h"

std::vector<std::vector<WyckoffType>> generate_isopointal_groups(
    const Shape& shape,
    const WallpaperGroup& group,
    std::size_t num_occupied_sites) {
  auto console = spdlog::stdout_color_mt("console");

  std::vector<WyckoffType> valid_sites;
  for (const WyckoffType& wyckoff : group.wyckoffs) {
    // first test if the shape has the required symmetries for various sites
    // at the moment only tests rotations & mirrors... that's all?
    int rotational_match{shape.rotational_symmetries % wyckoff.site_rotations};
    int mirror_match{shape.mirrors % wyckoff.site_mirrors};
    if (rotational_match == 0) {
      if ((wyckoff.site_mirrors == 0) ||
          ((mirror_match == 0) && (shape.mirrors != 0))) {
        if (wyckoff.some_variability) {
          console->debug("site %c is valid and variable\n", wyckoff.letter);
        } else {
          printf("site %c is valid\n", wyckoff.letter);
        }
        if (wyckoff.some_variability) {
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

  std::vector<std::vector<WyckoffType>> occupied_sites = combinations<WyckoffType>(
      valid_sites.begin(), valid_sites.end(), num_occupied_sites);

  console->info(
      "enumeration complete: in fact there were %lu ways\n", occupied_sites.size());

  for (const auto& combination : occupied_sites) {
    for (const auto& site : combination) {
      console->debug(" %c", site.letter);
    }
    console->debug("\n");
  }
  return occupied_sites;
}
