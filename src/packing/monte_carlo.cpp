/*
 * monte_carlo.cpp
 * Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "monte_carlo.h"

#include <sstream>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "wallpaper.h"

double MCVars::kT_ratio() const {
  return std::pow(this->kT_finish / this->kT_start, 1.0 / this->steps);
};

std::ostream& operator<<(std::ostream& os, const PackedState& packed_state) {
  os << "Shape: " << packed_state.shape->name << std::endl;
  os << "Cell:" << std::endl;
  os << "  a: " << packed_state.cell->x_len->get_value() << std::endl;
  os << "  b: " << packed_state.cell->y_len->get_value() << std::endl;
  os << "  angle: " << packed_state.cell->angle->get_value() << std::endl;
  os << "Wallpaper Group: " << packed_state.wallpaper->label << std::endl;
  for (const auto& site : *packed_state.occupied_sites) {
    os << site << std::endl;
  }
}

std::string PackedState::str() const {
  std::stringstream ss;
  ss << *this;
  return ss.str();
}

PackedState uniform_best_packing_in_isopointal_group(
    const Shape& shape,
    const WallpaperGroup& wallpaper,
    const IsopointalGroup& isopointal,
    const MCVars& mc_vars) {
  auto console = spdlog::stdout_color_mt("console");
  std::vector<Basis> best_basis;
  std::vector<bool> best_flips;

  std::size_t count_replicas{0};

  double packing{-1};
  double packing_max{0};

  /* Each cycle starts with a new random initialisation */
  std::size_t monte_carlo_steps{0};
  std::size_t rejections{0};
  double kT{mc_vars.kT_start};

  double packing_prev;

  std::vector<Basis> basis;
  Cell cell;
  std::vector<OccupiedSite> sites = initialise_structure(
      shape, isopointal, wallpaper, cell, basis, mc_vars.max_step_size);
  FlipBasis flip_basis{sites};

  packing = calculate_packing_fraction(shape, cell, sites);
  console->info("Initial packing fraction = %f\n", packing);

  while (monte_carlo_steps < mc_vars.steps) {
    kT *= mc_vars.kT_ratio();

    std::size_t vary_index{rand() % basis.size()};
    Basis& basis_current = basis[vary_index];

    /* Occasionally allow flips */
    if (monte_carlo_steps % 100) {
      const double flip_index{flip_basis.get_random_value(kT)};
      flip_basis.set_value(flip_index);
    }

    packing_prev = packing;
    const double new_value{basis_current.get_random_value(kT)};
    basis_current.set_value(new_value);

    if (check_state_for_intersection(shape, sites, cell)) {
      rejections++;
      basis_current.reset_value();
    } else {
      packing = calculate_packing_fraction(shape, cell, sites);
      if (fluke() >
          temperature_distribution(packing_prev, packing, kT, count_replicas)) {
        rejections++;
        basis_current.reset_value();
        flip_basis.reset_value();
        packing = packing_prev;
      }

      if (packing > packing_max) {
        /* best packing seen yet ... save data */
        best_basis = std::vector<Basis>(basis);
        best_flips = std::vector<bool>();
        for (const OccupiedSite& site : sites) {
          best_flips.push_back(site.flip_site);
        }
      }

      if (monte_carlo_steps % 500 == 0) {
        console->debug(
            "step %ld of %d, kT=%g, packing %f, angle %f, "
            "b/a=%f, rejection %f percent\n",
            monte_carlo_steps,
            mc_vars.steps,
            kT,
            packing,
            cell.angle->get_value() * 180.0 / M_PI,
            cell.x_len->get_value() / cell.y_len->get_value(),
            (100.0 * rejections) / monte_carlo_steps);
      }
    }

    packing = calculate_packing_fraction(shape, cell, sites);
    console->info(
        "BEST: cell %f %f angle %6.2f packing %f rejection (%f\%%) ",
        cell.x_len->get_value(),
        cell.y_len->get_value(),
        cell.angle->get_value() * 180.0 / M_PI,
        packing_max,
        (100.0 * rejections) / monte_carlo_steps);
    for (const Basis& b : basis) {
      printf("%f ", b.get_value());
    }
    for (const bool flipped : best_flips) {
      printf("%d ", flipped);
    }
    printf("\n");
  }

  return PackedState(wallpaper, shape, cell, sites);
}
