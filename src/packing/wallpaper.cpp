#include "wallpaper.h"

// To find orientation of ordered triplet (a, b, c).
// The function returns following values
// 0 --> a, b and c are colinear
// 1 --> Clockwise
// -1 --> Counterclockwise
int triplet_orientation(const Vect2& a, const Vect2& b, const Vect2& c) {
  // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
  // for details of below formula.
  Vect2 tmp = (b - a) * (c - b);
  return sign(tmp.x - tmp.y);
}

// Given three colinear points a, b, c, the function checks if
// point b lies on line segment 'ac'
bool on_segment(const Vect2& a, const Vect2& b, const Vect2& c) {
  if (b.x <= std::max(a.x, c.x) && b.x >= std::min(a.x, c.x) &&
      b.y <= std::max(a.y, c.y) && b.y >= std::min(a.y, c.y))
    return true;

  return false;
}

bool segments_cross(
    const Vect2& A1,
    const Vect2& A2,
    const Vect2& B1,
    const Vect2& B2) {

  // Find the four orientations needed for general and
  // special cases
  int o1 = triplet_orientation(A1, B1, A2);
  int o2 = triplet_orientation(A1, B1, B2);
  int o3 = triplet_orientation(A2, B2, A1);
  int o4 = triplet_orientation(A2, B2, B1);

  // General case
  if (o1 != o2 && o3 != o4)
    return true;

  // Special Cases
  // A1, B1 and A2 are colinear and A2 lies on segment p1q1
  if (o1 == 0 && on_segment(A1, A2, B1))
    return true;

  // p1, q1 and q2 are colinear and q2 lies on segment p1q1
  if (o2 == 0 && on_segment(A1, B2, B1))
    return true;

  // p2, q2 and p1 are colinear and p1 lies on segment p2q2
  if (o3 == 0 && on_segment(A2, A1, B2))
    return true;

  // p2, q2 and q1 are colinear and q1 lies on segment p2q2
  if (o4 == 0 && on_segment(A2, B1, B2))
    return true;

  return false; // Doesn't fall in any of the above cases
}

std::vector<Vect2> generate_position_cache(
    const Shape& shape,
    const Vect2& position,
    double angle_to_shape) {
  int resolution = shape.resolution();

  std::vector<Vect2> position_cache;
  // Reserve the expected size of the vector on initialisation
  position_cache.reserve(resolution / 2 + 1);

  double angular_step = 2.0 * PI / resolution;
  int angle_int = static_cast<int>(std::round(angle_to_shape / angular_step));
  // Flipping reverses the direction of the points

  for (int index = -resolution / 4; index <= (resolution / 4); index++) {
    int compare_index = angle_int + index;
    // Change base to angle between shapes
    double theta = fabs(compare_index * angular_step - angle_to_shape);
    position_cache.push_back(Vect2(
        shape.get_point(compare_index) * cos(theta),
        shape.get_point(compare_index) * sin(theta)));
  }
  return position_cache;
}

bool indirect_clash(
    const Shape& shape_a,
    const Shape& shape_b,
    const Vect2& position_a,
    const Vect2& position_b,
    double angle_a_to_b,
    double angle_b_to_a) {

  std::vector<Vect2> position_a_cache =
      generate_position_cache(shape_a, position_a, angle_a_to_b);
  std::vector<Vect2> position_b_cache =
      generate_position_cache(shape_b, position_b, angle_b_to_a);

  // The final element is the inital previous position providing a closed shape
  // regardless of the number of points checked.
  Vect2& position_a_prev = position_a_cache.back();
  Vect2& position_b_prev = position_b_cache.back();
  for (auto position_a = position_a_cache.begin(); position_a != position_a_cache.end();
       ++position_a) {
    for (auto position_b = position_b_cache.begin();
         position_b != position_b_cache.end();
         ++position_b) {
      if (segments_cross(position_a_prev, *position_a, position_b_prev, *position_b)) {
        return true;
      }
      // Update the previous point
      position_a_prev = *position_a;
      position_b_prev = *position_b;
    }
  }
  return false;
}

bool pair_clash(
    const ShapeInstance& shape_a,
    const ShapeInstance& shape_b,
    const Vect2& coords_a,
    const Vect2& coords_b) {

  double mathtol = 1e-8;

  double central_dist = (coords_a - coords_b).norm();

  /* No clash when further apart than the maximum shape radii measures */
  if (central_dist > shape_a.shape->max_radius + shape_b.shape->max_radius) {
    return false;
  }

  /* we need a polygon overlap calculation */
  /* first calculate the direct overlap */
  double a_to_b_incline = acos((coords_b.x - coords_a.x) / central_dist);

  if (std::isnan(a_to_b_incline)) {
    if (is_close(std::pow(coords_b.x - coords_a.x, 2), central_dist, mathtol)) {
      a_to_b_incline = 0.0;
    } else if (is_close(coords_b.x - coords_a.x, -central_dist, mathtol)) {
      a_to_b_incline = M_PI;
    }
  }
  if (coords_b.x < coords_a.x) {
    a_to_b_incline = 2.0 * M_PI - a_to_b_incline;
  }

  // Set reverse incline
  double b_to_a_incline = a_to_b_incline + M_PI;

  // Deal with the flipping of shapes
  if (shape_a.get_flipped()) {
    a_to_b_incline = 2 * M_PI - a_to_b_incline;
  }
  if (shape_b.get_flipped()) {
    b_to_a_incline = 2 * M_PI - b_to_a_incline;
  }

  /* now add in the rotation due to the orientation parameters */
  /* This is unaffected by the flip states */
  a_to_b_incline += shape_a.get_angle();
  b_to_a_incline += shape_b.get_angle();

  /* now add in the rotation due to the rotation of this image wrt the other
   * images of the same wyckoff */
  a_to_b_incline +=
      std::pow(-1, shape_a.get_flipped()) * shape_a.get_rotational_offset();
  b_to_a_incline +=
      std::pow(-1, shape_b.get_flipped()) * shape_b.get_rotational_offset();

  a_to_b_incline = positive_modulo(a_to_b_incline, 2.0 * PI);
  b_to_a_incline = positive_modulo(b_to_a_incline, 2.0 * PI);

  return indirect_clash(
      *shape_a.shape,
      *shape_b.shape,
      coords_a,
      coords_b,
      a_to_b_incline,
      b_to_a_incline);
}

bool clash_polygon(
    const ShapeInstance& shape_a,
    const ShapeInstance& shape_b,
    const Cell& cell) {

  Vect2 fcoords_a = shape_a.get_fractional_coordinates();
  Vect2 fcoords_b = shape_b.get_fractional_coordinates();

  /* a is fixed, copies of b are made to test for the clash */
  Vect2 img_fcoords_b(0, 0);

  /*
   * The method of only looking at the 3x3 nearest cells fails for those with
   * extreme angles, so is expanded to 5x5.
   */
  for (int cell_img_x = -2; cell_img_x <= 2; cell_img_x++) {

    img_fcoords_b.x = fcoords_b.x + cell_img_x;

    for (int cell_img_y = -2; cell_img_y <= 2; cell_img_y++) {

      // Intersections with one's self are excluded
      if ((shape_a == shape_b) && (cell_img_x == 0) && (cell_img_y == 0)) {
        continue;
      }

      img_fcoords_b.y = fcoords_b.y + cell_img_y;

      Vect2 coords_a = cell.fractional_to_real(fcoords_a);
      Vect2 coords_b = cell.fractional_to_real(img_fcoords_b);

      if (pair_clash(shape_a, shape_b, coords_a, coords_b)) {
        return 1;
      }
    }
  }
  return 0;
}

double calculate_packing_fraction(
    Shape& shape,
    Cell& cell,
    std::vector<Site>& occupied_sites) {

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
    double step_size) {

  // Logging to console which can be turned off easily
  auto console = spdlog::stdout_color_mt("console");

  int count_replicas = group_multiplicity(occupied_sites);

  // cell sides.
  double max_cell_size = 4 * shape.max_radius * count_replicas;
  if (group.a_b_equal) {
    console->debug("Cell sides equal");
    basis.push_back(CellLengthBasis(max_cell_size, 0.1, max_cell_size, step_size));

    cell.x_len = &basis.back();
    cell.y_len = &basis.back();
  } else {
    basis.push_back(CellLengthBasis(max_cell_size, 0.1, max_cell_size, step_size));
    cell.x_len = &basis.back();

    basis.push_back(CellLengthBasis(max_cell_size, 0.1, max_cell_size, step_size));
    cell.y_len = &basis.back();
  }

  // cell angles.
  if (group.hexagonal) {
    console->debug("Hexagonal group");
    cell.angle = new FixedBasis(PI / 3);
  } else if (group.rectangular) {
    console->debug("Rectangular group");
    cell.angle = new FixedBasis(PI / 2);
  } else {
    console->debug("Tilted group");
    basis.push_back(CellAngleBasis(
        PI / 4 * fluke() * PI / 2,
        PI / 4,
        3 * PI / 4,
        step_size,
        cell.x_len,
        cell.y_len));
    cell.angle = &basis.back();
  }

  // now position the particles.
  for (Site& site : occupied_sites) {
    count_replicas += site.wyckoff->multiplicity;

    console->debug("Wyckoff site: %c ", site.wyckoff->letter);

    if (fabs(site.wyckoff->image[0].x_coeffs.x) > 0.1) {
      /* x is variable*/
      basis.push_back(Basis(fluke(), 0, 1));
      site.x = &basis.back();
      console->debug("Site x variable %f\n", site.x->get_value());
    }
    if (fabs(site.wyckoff->image[0].y_coeffs.y) > 0.1) {
      /* then y is variable*/
      basis.push_back(Basis(fluke(), 0, 1));
      site.y = &basis.back();
      console->debug("Site y variable %f\n", site.y->get_value());
    }

    /*choose the orientation of the zeroth image of this particle */
    /* This could also be done within the Wyckoff position SITEROTATION
     * parameter?... */
    if (site.wyckoff->site_mirrors) {
      int mirrors = site.wyckoff->image[0].site_mirror;
      double value = M_PI / 180 * mirrors;
      basis.push_back(MirrorBasis(value, 0, 2 * PI, mirrors));
      site.angle = &basis.back();
    } else {
      double value = fluke() * 2.0 * M_PI;
      basis.push_back(Basis(value, 0, 2 * PI, step_size));
      site.angle = &basis.back();
      console->debug("site offset-angle is variable %f\n", site.angle->get_value());
    }
  }

  console->debug("replicas %d variables %d os ", count_replicas, basis.size());

  return basis.size();
}

bool there_is_collision() { return true; }

void uniform_best_packing_in_isopointal_group(
    Shape& shape,
    WallpaperGroup& group,
    size_t num_cycles,
    size_t max_steps,
    double max_step_size,
    double kT_start = 0.1,
    double kT_finish = 5e-4) {
  auto console = spdlog::stdout_color_mt("console");
  std::vector<Basis> best_basis;
  std::vector<bool> best_flips;

  size_t count_replicas = 0;

  double packing_fraction = -1;
  double packing_fraction_max = 0.0;

  double kT_ratio = std::pow(kT_finish / kT_start, 1.0 / max_steps);

  /* Each cycle starts with a new random initialisation */
  size_t monte_carlo_steps = 0;
  size_t rejections = 0;
  double kT = kT_start;

  double packing_fraction_prev;

  std::vector<Basis> basis;
  std::vector<Site> occupied_sites;
  Cell cell;
  FlipBasis flip_basis{&occupied_sites};

  initialize_structure_in_group(
      shape, group, cell, occupied_sites, basis, max_step_size);

  packing_fraction = calculate_packing_fraction(shape, cell, occupied_sites);

  console->info("Initial packing fraction = %f\n", packing_fraction);

  while (monte_carlo_steps < max_steps) {
    kT *= kT_ratio;

    size_t vary_index = rand() % basis.size();
    Basis& basis_current = basis[vary_index];

    /* Occasionally allow flips */
    if (monte_carlo_steps % 100) {
      double flip_index = flip_basis.get_random_value(kT);
      flip_basis.set_value(flip_index);
    }

    packing_fraction_prev = packing_fraction;
    double new_value = basis_current.get_random_value(kT);
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
  for (auto site : occupied_sites) {
    chiralsum += site.get_flip_sign() * site.get_multiplicity();
    totalsum += site.get_multiplicity();
  }
  if (chiralsum != 0) {
    return (chiralsum == totalsum) ? 'c' : 's';
  }
  return 'a';
}
