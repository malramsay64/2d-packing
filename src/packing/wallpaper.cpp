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
#include "util.h"

/** Check whether two shape instances intersect
 */
bool check_for_intersection(
    const ShapeInstance& shape_a,
    const ShapeInstance& shape_b,
    const Cell& cell) {

  // a is fixed, b is moved to the periodic sites to test for the intersection
  Vect2 fcoords_b{shape_b.get_fractional_coordinates()};
  Vect2 img_fcoords_b{0, 0};

  int shells = 1;
  // For extreme angles, using only the nearest shell fails, so have to look at 2
  // shells. The designator for 'extreme' angle is PI/4 or 45 degrees.
  if (cell.angle->get_value() < M_PI_4) {
    shells = 2;
  } else if (M_2_PI - cell.angle->get_value() < M_PI_4) {
    shells = 2;
  }

  // Loop over the possible periodic positions
  for (int cell_img_x = -shells; cell_img_x <= shells; cell_img_x++) {
    for (int cell_img_y = -shells; cell_img_y <= shells; cell_img_y++) {
      // Intersections with one's self are excluded
      if ((shape_a == shape_b) && (cell_img_x == 0) && (cell_img_y == 0)) {
        continue;
      }
      Vect2 coords_b = cell.fractional_to_real(
          Vect2(fcoords_b.x + cell_img_x, fcoords_b.y + cell_img_y));
      if (shape_a.intersects_with(shape_b, coords_b)) {
        return true;
      }
    }
  }
  return false;
}

bool check_state_for_intersection(
    const Shape& shape,
    const std::vector<OccupiedSite>& occupied_sites,
    const Cell& cell) {
  // Loop over all the occupied sites
  for (auto site_one = occupied_sites.begin(); site_one != occupied_sites.end();
       site_one++) {
    // Loop over all symmetries for the first occupied site
    for (const auto& image_one : site_one->wyckoff->symmetries) {
      const ShapeInstance shape_one{shape, *site_one, image_one};
      // Loop over all occupied sites which haven't already been compared with site_one
      for (auto site_two = std::next(site_one); site_two != occupied_sites.end();
           site_two++) {
        // Loop over all symmetries for the second occupied site
        for (const auto& image_two : site_two->wyckoff->symmetries) {
          const ShapeInstance shape_two{shape, *site_two, image_two};
          /* Finally perform the comparison of shapes here */
          if (check_for_intersection(shape_one, shape_two, cell)) {
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

double calculate_packing_fraction(
    const Shape& shape,
    const Cell& cell,
    const std::vector<OccupiedSite>& occupied_sites) {
  std::size_t shape_replicas(calculate_shape_replicas(occupied_sites));

  double packing_fraction = shape_replicas * shape.area() / cell.area();
  if (std::isnan(packing_fraction)) {
    auto console = spdlog::stdout_color_mt("console");
    console->warn(
        "nan encountered %f %f %f\n",
        cell.x_len->get_value(),
        cell.y_len->get_value(),
        cell.angle->get_value());
    throw "Packing fraction calculation failed";
  }

  return packing_fraction;
}

std::vector<OccupiedSite> initialise_structure(
    const Shape& shape,
    const IsopointalGroup& isopointal_group,
    const WallpaperGroup& group,
    Cell& cell,
    std::vector<Basis>& basis,
    const double step_size) {

  // Logging to console which can be turned off easily
  auto console = spdlog::stdout_color_mt("console");

  // cell sides.
  std::size_t count_replicas{isopointal_group.group_multiplicity()};
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

  std::vector<OccupiedSite> sites;
  // now position the particles.
  for (const WyckoffSite& wyckoff : isopointal_group.wyckoff_sites) {
    OccupiedSite site{};

    console->debug("Wyckoff site: %c ", wyckoff.letter);

    if (wyckoff.vary_x()) {
      /* x is variable*/
      basis.push_back(Basis(fluke(), 0, 1));
      site.x = std::shared_ptr<Basis>(&basis.back());
      console->debug("WyckoffSite x variable %f\n", site.x->get_value());
    }
    if (wyckoff.vary_y()) {
      /* then y is variable*/
      basis.push_back(Basis(fluke(), 0, 1));
      site.y = std::shared_ptr<Basis>(&basis.back());
      console->debug("WyckoffSite y variable %f\n", site.y->get_value());
    }

    /*choose the orientation of the zeroth image of this particle */
    /* This could also be done within the Wyckoff position SITEROTATION
     * parameter?... */
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

  return sites;
}
