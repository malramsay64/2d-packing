/*
 * monte_carlo.h
 * Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <memory>
#include <vector>

#include <pybind11/pybind11.h>

#include "basis.h"
#include "shapes.h"
#include "wallpaper.h"

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

struct MCVars {
  double kT_start = 0.1;
  double kT_finish = 5e-4;
  double max_step_size = 0.01;
  std::size_t num_cycles = 32;
  std::size_t steps = 10000;

  double kT_ratio() const;
};

class PackedState {
  const std::shared_ptr<const WallpaperGroup> wallpaper;
  const std::shared_ptr<const Shape> shape;
  const std::shared_ptr<Cell> cell;
  const std::shared_ptr<std::vector<OccupiedSite>> occupied_sites;
  const std::shared_ptr<std::vector<Basis>> basis;

public:
  PackedState(
      std::shared_ptr<const WallpaperGroup> wallpaper,
      std::shared_ptr<const Shape> shape,
      std::shared_ptr<Cell> cell,
      std::shared_ptr<std::vector<OccupiedSite>> occupied_sites,
      std::shared_ptr<std::vector<Basis>> basis)

      : wallpaper(wallpaper), shape(shape), cell(cell), occupied_sites(occupied_sites),
        basis(basis){};

  PackedState(
      const WallpaperGroup& wallpaper,
      const Shape& shape,
      Cell& cell,
      std::vector<OccupiedSite>& occupied_sites,
      std::vector<Basis> basis)
      : PackedState(
            std::shared_ptr<const WallpaperGroup>(&wallpaper),
            std::shared_ptr<const Shape>(&shape),
            std::shared_ptr<Cell>(&cell),
            std::shared_ptr<std::vector<OccupiedSite>>(&occupied_sites),
            std::shared_ptr<std::vector<Basis>>(&basis)){};

  PackedState(
      const WallpaperGroup& wallpaper,
      const Shape& shape,
      Cell& cell,
      std::vector<OccupiedSite>& occupied_sites)
      : PackedState(wallpaper, shape, cell, occupied_sites, std::vector<Basis>()){};

  std::string str() const;

  friend std::ostream& operator<<(std::ostream& os, const PackedState& packed_state);
};

PackedState initialise_structure(
    const Shape& shape,
    const IsopointalGroup& isopointal,
    const WallpaperGroup& wallpaper,
    const double step_size);

void export_PackedState(pybind11::module& m);

#endif /* !MONTE_CARLO_H */
