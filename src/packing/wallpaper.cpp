/*
 * wallpaper.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "wallpaper.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "basis.h"
#include "geometry.h"
#include "math.h"
#include "random.h"
#include "shapes.h"

bool clash_polygon(
    const ShapeInstance& shape_a,
    const ShapeInstance& shape_b,
    const Cell& cell) {

  Vect2 fcoords_b = shape_b.get_fractional_coordinates();

  // a is fixed, copies of b are made to test for the clash
  Vect2 img_fcoords_b(0, 0);
  Vect2 coords_a = cell.fractional_to_real(shape_a.get_fractional_coordinates());

  int shells = 1;
  // For extreme angles, only the nearest shell fails, so have to look at 2 shells.
  // The designator for 'extreme' angle is PI/4 or 45 degrees.
  if (cell.angle->get_value() < M_PI_4) {
    shells = 2;
  } else if (M_2_PI - cell.angle->get_value() < M_PI_4) {
    shells = 2;
  }

  for (int cell_img_x = -shells; cell_img_x <= shells; cell_img_x++) {
    for (int cell_img_y = -shells; cell_img_y <= shells; cell_img_y++) {
      // Intersections with one's self are excluded
      if ((shape_a == shape_b) && (cell_img_x == 0) && (cell_img_y == 0)) {
        continue;
      }
      //
      Vect2 coords_b = cell.fractional_to_real(
          Vect2(fcoords_b.x + cell_img_x, fcoords_b.y + cell_img_y));
      if (shape_a.intersects_with(shape_b, coords_b)) {
        return true;
      }
    }
  }
  return false;
}

double calculate_packing_fraction(
    const Shape& shape,
    const Cell& cell,
    const std::vector<Site>& occupied_sites) {

  auto console = spdlog::stdout_color_mt("console");

  int count_replicas = 0;
  for (const Site& site : occupied_sites) {
    count_replicas += site.wyckoff->multiplicity;
  }

  double packing_fraction = count_replicas * shape.area() / cell.area();
  if (std::isnan(packing_fraction)) {

    console->warn(
        "nan encountered %f %f %f\n",
        cell.x_len->get_value(),
        cell.y_len->get_value(),
        cell.angle->get_value());
    exit(1);
  }

  return packing_fraction;
}

int initialize_structure_in_group(
    Shape& shape,
    WallpaperGroup& group,
    Cell& cell,
    std::vector<Site>& occupied_sites,
    std::vector<Basis>& basis,
    const double step_size) {

  // Logging to console which can be turned off easily
  auto console = spdlog::stdout_color_mt("console");

  int count_replicas = group_multiplicity(occupied_sites);

  // cell sides.
  const double max_cell_size{4 * shape.max_radius * count_replicas};
  if (group.a_b_equal) {
    console->debug("Cell sides equal");
    basis.push_back(CellLengthBasis(max_cell_size, 0.1, max_cell_size, step_size));

    cell.x_len = std::shared_ptr<Basis>(&basis.back());
    cell.y_len = std::shared_ptr<Basis>(&basis.back());
  } else {
    basis.push_back(CellLengthBasis(max_cell_size, 0.1, max_cell_size, step_size));
    cell.x_len = std::shared_ptr<Basis>(&basis.back());

    basis.push_back(CellLengthBasis(max_cell_size, 0.1, max_cell_size, step_size));
    cell.y_len = std::shared_ptr<Basis>(&basis.back());
  }

  // cell angles.
  if (group.hexagonal) {
    console->debug("Hexagonal group");
    cell.angle = std::make_shared<FixedBasis>(M_PI / 3);
  } else if (group.rectangular) {
    console->debug("Rectangular group");
    cell.angle = std::make_shared<FixedBasis>(M_PI_2);
  } else {
    console->debug("Tilted group");
    basis.push_back(CellAngleBasis(
        M_PI_4 + fluke() * M_PI_2,
        M_PI_4,
        3 * M_PI_4,
        step_size,
        cell.x_len,
        cell.y_len));
    cell.angle = std::shared_ptr<Basis>(&basis.back());
  }

  // now position the particles.
  for (Site& site : occupied_sites) {
    count_replicas += site.wyckoff->multiplicity;

    console->debug("Wyckoff site: %c ", site.wyckoff->letter);

    if (fabs(site.wyckoff->image[0].x_coeffs.x) > 0.1) {
      /* x is variable*/
      basis.push_back(Basis(fluke(), 0, 1));
      site.x = std::shared_ptr<Basis>(&basis.back());
      console->debug("Site x variable %f\n", site.x->get_value());
    }
    if (fabs(site.wyckoff->image[0].y_coeffs.y) > 0.1) {
      /* then y is variable*/
      basis.push_back(Basis(fluke(), 0, 1));
      site.y = std::shared_ptr<Basis>(&basis.back());
      console->debug("Site y variable %f\n", site.y->get_value());
    }

    /*choose the orientation of the zeroth image of this particle */
    /* This could also be done within the Wyckoff position SITEROTATION
     * parameter?... */
    if (site.wyckoff->site_mirrors) {
      const int mirrors{site.wyckoff->image[0].site_mirror};
      const double value{M_PI / 180 * mirrors};
      basis.push_back(MirrorBasis(value, 0, 2 * PI, mirrors));
      site.angle = std::shared_ptr<Basis>(&basis.back());
    } else {
      const double value{fluke() * M_2_PI};
      basis.push_back(Basis(value, 0, M_2_PI, step_size));
      site.angle = std::shared_ptr<Basis>(&basis.back());
      console->debug("site offset-angle is variable %f\n", site.angle->get_value());
    }
  }

  console->debug("replicas %d variables %d os ", count_replicas, basis.size());

  return basis.size();
}

bool there_is_collision() {
  return true;
}

void uniform_best_packing_in_isopointal_group(
    Shape& shape,
    WallpaperGroup& group,
    const std::size_t num_cycles,
    const std::size_t max_steps,
    const double max_step_size,
    const double kT_start = 0.1,
    const double kT_finish = 5e-4) {
  auto console = spdlog::stdout_color_mt("console");
  std::vector<Basis> best_basis;
  std::vector<bool> best_flips;

  std::size_t count_replicas{0};

  double packing_fraction{-1};
  double packing_fraction_max{0};

  const double kT_ratio{std::pow(kT_finish / kT_start, 1.0 / max_steps)};

  /* Each cycle starts with a new random initialisation */
  std::size_t monte_carlo_steps{0};
  std::size_t rejections{0};
  double kT{kT_start};

  double packing_fraction_prev;

  std::vector<Basis> basis;
  std::vector<Site> occupied_sites;
  Cell cell;
  FlipBasis flip_basis{occupied_sites};

  initialize_structure_in_group(
      shape, group, cell, occupied_sites, basis, max_step_size);

  packing_fraction = calculate_packing_fraction(shape, cell, occupied_sites);

  console->info("Initial packing fraction = %f\n", packing_fraction);

  while (monte_carlo_steps < max_steps) {
    kT *= kT_ratio;

    std::size_t vary_index{rand() % basis.size()};
    Basis& basis_current = basis[vary_index];

    /* Occasionally allow flips */
    if (monte_carlo_steps % 100) {
      const double flip_index{flip_basis.get_random_value(kT)};
      flip_basis.set_value(flip_index);
    }

    packing_fraction_prev = packing_fraction;
    const double new_value{basis_current.get_random_value(kT)};
    basis_current.set_value(new_value);

    if (there_is_collision()) {
      rejections++;
      basis_current.reset_value();
    } else {
      packing_fraction = calculate_packing_fraction(shape, cell, occupied_sites);
      if (fluke() > temperature_distribution(
                        packing_fraction_prev, packing_fraction, kT, count_replicas)) {
        rejections++;
        basis_current.reset_value();
        flip_basis.reset_value();
        packing_fraction = packing_fraction_prev;
      }

      if (packing_fraction > packing_fraction_max) {
        /* best packing seen yet ... save data */
        best_basis = std::vector<Basis>(basis);
        best_flips = std::vector<bool>();
        for (const Site& site : occupied_sites) {
          best_flips.push_back(site.flip_site);
        }
      }

      if (monte_carlo_steps % 500 == 0) {
        console->debug(
            "step %ld of %d, kT=%g, packing %f, angle %f, "
            "b/a=%f, rejection %f percent\n",
            monte_carlo_steps,
            max_steps,
            kT,
            packing_fraction,
            cell.angle->get_value() * 180.0 / M_PI,
            cell.x_len->get_value() / cell.y_len->get_value(),
            (100.0 * rejections) / monte_carlo_steps);
      }
    }

    packing_fraction = calculate_packing_fraction(shape, cell, occupied_sites);
    console->info(
        "BEST: cell %f %f angle %6.2f packing %f rejection (%f\%%) ",
        cell.x_len->get_value(),
        cell.y_len->get_value(),
        cell.angle->get_value() * 180.0 / M_PI,
        packing_fraction_max,
        (100.0 * rejections) / monte_carlo_steps);
    for (const Basis& b : basis) {
      printf("%f ", b.get_value());
    }
    for (const bool flipped : best_flips) {
      printf("%d ", flipped);
    }
    printf("\n");
  }
}

char compute_chiral_state(const std::vector<Site>& occupied_sites) {
  int chiralsum = 0;
  int totalsum = 0;
  for (const Site& site : occupied_sites) {
    chiralsum += site.get_flip_sign() * site.get_multiplicity();
    totalsum += site.get_multiplicity();
  }
  if (chiralsum != 0) {
    return (chiralsum == totalsum) ? 'c' : 's';
  }
  return 'a';
}
