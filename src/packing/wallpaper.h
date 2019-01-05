/*
 * wallpaper.h
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <vector>

#include "isopointal.h"
#include "shapes.h"

#ifndef WALLPAPER_H
#define WALLPAPER_H

class WallpaperGroup {
public:
  WallpaperGroup(std::string label, std::vector<WyckoffSite> wyckoffs);
  const std::string label;
  const bool a_b_equal = false;
  const bool hexagonal = false;
  const bool rectangular = false;
  int num_symmetries = 0;
  std::vector<WyckoffSite> wyckoff_sites;
  int num_wyckoffs;
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

#endif /* !WALLPAPER_H */
