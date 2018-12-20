/*
 * wallpaper.h
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef WALLPAPER_H
#define WALLPAPER_H

#include "random.h"
#include "globals.h"
#include "shapes.h"
#include <cmath>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <string>
#include <vector>

inline bool is_close(float value, float expected, float rel_tol = 1e-8) {
  return fabs(value - expected) < rel_tol * expected;
}

#define MAXNUMOCCSITES 3

#endif /* !WALLPAPER_H */
