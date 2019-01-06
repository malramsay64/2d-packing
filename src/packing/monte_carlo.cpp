/*
 * monte_carlo.cpp
 * Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "monte_carlo.h"

#include <sstream>

#include <pybind11/pybind11.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "packing.h"
#include "wallpaper.h"

namespace py = pybind11;

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
  return os;
}

std::string PackedState::str() const {
  std::stringstream ss;
  ss << *this;
  return ss.str();
}

double PackedState::packing_fraction() const {
  return this->num_shapes() * this->shape->area() / this->cell->area();
};

std::size_t PackedState::num_shapes() const {
  std::size_t num_shapes{0};
  for (const auto& site : *this->occupied_sites) {
    num_shapes += site.get_multiplicity();
  }
  return num_shapes;
}

bool PackedState::check_intersection() const {
  // Loop over all the occupied sites
  for (auto site_one = this->occupied_sites->begin();
       site_one != this->occupied_sites->end();
       site_one++) {
    // Loop over all symmetries for the first occupied site
    for (const auto& symmetry_one : site_one->wyckoff->symmetries) {
      const ShapeInstance shape_one{*this->shape, *site_one, symmetry_one};
      // Loop over all occupied sites which haven't already been compared with site_one
      for (auto site_two = std::next(site_one); site_two != this->occupied_sites->end();
           site_two++) {
        // Loop over all symmetries for the second occupied site
        for (const auto& image_two : site_two->wyckoff->symmetries) {
          const ShapeInstance shape_two{*this->shape, *site_two, image_two};
          /* Finally perform the comparison of shapes here */
          if (check_for_intersection(shape_one, shape_two, *this->cell)) {
            // If the two shapes intersect, return true, breaking out of the loop.
            return true;
          }
        }
      }
    }
  }
  // Should there be no intersections between any shapes, return false
  return false;
}

PackedState initialise_structure(
    const Shape& shape,
    const IsopointalGroup& isopointal,
    const WallpaperGroup& wallpaper,
    const double step_size) {

  // Logging to console which can be turned off easily
  auto console = spdlog::stdout_color_mt("console");

  std::vector<Basis> basis;
  Cell cell;

  // cell sides.
  std::size_t count_replicas{isopointal.group_multiplicity()};
  const double max_cell_size{4 * shape.max_radius * count_replicas};
  if (wallpaper.a_b_equal) {
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
  if (wallpaper.hexagonal) {
    console->debug("Hexagonal group");
    cell.angle = std::make_shared<FixedBasis>(M_PI / 3);
  } else if (wallpaper.rectangular) {
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

  std::vector<OccupiedSite> sites;
  // The chosen Wyckoff sites are in the IsopointalGroup class.
  for (const WyckoffSite& wyckoff : isopointal.wyckoff_sites) {
    OccupiedSite site{};

    console->debug("Wyckoff site: %c ", wyckoff.letter);

    // x is not fixed
    if (wyckoff.vary_x()) {
      basis.push_back(Basis(fluke(), 0, 1));
      site.x = std::shared_ptr<Basis>(&basis.back());
      console->debug("WyckoffSite x variable %f\n", site.x->get_value());
    }
    // y is not fixed
    if (wyckoff.vary_y()) {
      /* then y is variable*/
      basis.push_back(Basis(fluke(), 0, 1));
      site.y = std::shared_ptr<Basis>(&basis.back());
      console->debug("WyckoffSite y variable %f\n", site.y->get_value());
    }

    // Setting the angle of the Wyckoff Site.
    // Where there are mirrors the angle has fewer orienations it is allowed to take.
    if (wyckoff.mirrors) {
      const int mirrors{wyckoff.mirror_type()};
      const double value{M_PI / 180 * mirrors};
      basis.push_back(MirrorBasis(value, 0, M_2_PI, mirrors));
      site.angle = std::shared_ptr<Basis>(&basis.back());
    } else {
      const double value{fluke() * M_2_PI};
      basis.push_back(Basis(value, 0, M_2_PI, step_size));
      site.angle = std::shared_ptr<Basis>(&basis.back());
      console->debug("site offset-angle is variable %f\n", site.angle->get_value());
    }
    sites.push_back(site);
  }

  console->debug("replicas %d variables %d os ", count_replicas, basis.size());

  return PackedState(wallpaper, shape, cell, sites, basis);
}

PackedState uniform_best_packing_in_isopointal_group(
    const Shape& shape,
    const WallpaperGroup& wallpaper,
    const IsopointalGroup& isopointal,
    const MCVars& mc_vars) {
  auto console = spdlog::stdout_color_mt("console");
  std::vector<double> best_values;
  std::vector<bool> best_flips;

  std::size_t count_replicas{0};

  double packing{-1};
  double packing_max{0};

  /* Each cycle starts with a new random initialisation */
  std::size_t monte_carlo_steps{0};
  std::size_t rejections{0};
  double kT{mc_vars.kT_start};

  double packing_prev;

  PackedState sim_state =
      initialise_structure(shape, isopointal, wallpaper, mc_vars.max_step_size);
  FlipBasis flip_basis{*sim_state.occupied_sites};

  packing = sim_state.packing_fraction();
  console->info("Initial packing fraction = %f\n", packing);

  while (monte_carlo_steps < mc_vars.steps) {
    kT *= mc_vars.kT_ratio();

    std::size_t vary_index{rand() % sim_state.basis->size()};
    Basis& basis_current = sim_state.basis->at(vary_index);

    /* Occasionally allow flips */
    if (monte_carlo_steps % 100) {
      const double flip_index{flip_basis.get_random_value(kT)};
      flip_basis.set_value(flip_index);
    }

    packing_prev = packing;
    const double new_value{basis_current.get_random_value(kT)};
    basis_current.set_value(new_value);

    if (sim_state.check_intersection()) {
      rejections++;
      basis_current.reset_value();
    } else {
      packing = sim_state.packing_fraction();
      if (fluke() >
          temperature_distribution(packing_prev, packing, kT, count_replicas)) {
        rejections++;
        basis_current.reset_value();
        flip_basis.reset_value();
        packing = packing_prev;
      }

      if (packing > packing_max) {
        best_values = sim_state.save_basis();
        /* best packing seen yet ... save data */
        best_flips = std::vector<bool>();
        for (const OccupiedSite& site : *sim_state.occupied_sites) {
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
            sim_state.cell->angle->get_value() * 180.0 / M_PI,
            sim_state.cell->x_len->get_value() / sim_state.cell->y_len->get_value(),
            (100.0 * rejections) / monte_carlo_steps);
      }
    }

    sim_state.load_basis(best_values);
    packing = sim_state.packing_fraction();
    console->info(
        "BEST: cell %f %f angle %6.2f packing %f rejection (%f\%%) ",
        sim_state.cell->x_len->get_value(),
        sim_state.cell->y_len->get_value(),
        sim_state.cell->angle->get_value() * 180.0 / M_PI,
        packing_max,
        (100.0 * rejections) / monte_carlo_steps);
  }

  return sim_state;
}

void export_PackedState(py::module& m) {
  py::class_<PackedState>(m, "PackedState").def("__str__", &PackedState::str);
}
