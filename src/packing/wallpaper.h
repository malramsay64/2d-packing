/*
 * wallpaper.h
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <vector>

#include "shapes.h"

#ifndef WALLPAPER_H
#define WALLPAPER_H

std::vector<std::vector<WyckoffType>> generate_isopointal_groups(
    const Shape& shape,
    const WallpaperGroup& group,
    std::size_t num_occupied_sites);

#endif /* !WALLPAPER_H */
