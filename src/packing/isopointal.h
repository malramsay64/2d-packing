/*
 * isopointal.h
 * Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "shapes.h"

#ifndef ISOPOINTAL_H
#define ISOPOINTAL_H

class WyckoffSite {
  double x;
  double y;
  double theta;
};

class IsopointalGroup {
  std::vector<WyckoffSite> occupied_sites;
};

std::vector<IsopointalGroup> generate_isopointal_groups(
    const Shape& shape,
    const WallpaperGroup& group,
    std::size_t num_occupied_sites);

#endif /* !ISOPOINTAL_H */
