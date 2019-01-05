/*
 * isopointal.h
 * Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */
#include <string>

#include "shapes.h"

#ifndef ISOPOINTAL_H
#define ISOPOINTAL_H

// Forward definition of WallpaperGroup
class WallpaperGroup;

/** \class IsopointalGroup
 *
 * An IsopointalGroup is a collection
 */
class IsopointalGroup {
public:
  std::vector<WyckoffSite> wyckoff_sites;

  IsopointalGroup(std::vector<WyckoffSite>& sites) : wyckoff_sites(sites){};

  std::size_t group_multiplicity() const;
  std::string group_string() const;
};

std::vector<IsopointalGroup> generate_isopointal_groups(
    const Shape& shape,
    const WallpaperGroup& group,
    std::size_t num_occupied_sites);

#endif /* !ISOPOINTAL_H */
