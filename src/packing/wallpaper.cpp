#include "wallpaper.h"

int enumerate(
    int OS[][MAXNUMOCCSITES],
    int maxcombinations,
    int n_comb,
    int os,
    int bag_i,
    int combination_bag[],
    int bag_size,
    int numoccsites) {
  int comb_new;
  int comb_count = 0;

  if (os == numoccsites)
    return 1;
  else {
    for (int i = bag_i; i < bag_size - (numoccsites - os - 1); i++) {
      if ((i != bag_i) && (combination_bag[i] == combination_bag[i - 1]))
        continue;
      comb_new = enumerate(
          OS,
          maxcombinations,
          n_comb + comb_count,
          os + 1,
          i + 1,
          combination_bag,
          bag_size,
          numoccsites);
      for (int j = 0; j < comb_new; j++) {
        OS[n_comb + comb_count + j][os] = combination_bag[i];
        // printf("saving [%d][%d]=%d (new %d, count %d)\n",
        // n_comb+j+count, os, i, new, count);
      }
      comb_count += comb_new;
    }
  }
  return comb_count;
}

Vect2 fractcoords(
    const Vect3& x_coeffs,
    const Vect3& y_coeffs,
    const Vect3 sitevariables) {
  /* converts site variables and wyckoff site coefficients into the location
   * of the actual wyckoff image in fractional coordinates */
  Vect2 v(0, 0);
  v.x = positive_modulo(
      x_coeffs.x * sitevariables.x + x_coeffs.y * sitevariables.y + x_coeffs.z, 1.);
  v.y = positive_modulo(
      y_coeffs.x * sitevariables.x + y_coeffs.y * sitevariables.y + y_coeffs.z, 1.);
  return v;
}

bool segments_cross(
    double Ax,
    double Ay,
    double Bx,
    double By,
    double Cx,
    double Cy,
    double Dx,
    double Dy) {
  //  Determines the intersection point of the line segment
  //  defined by points A and B with the line segment defined
  //  by points C and D.
  //
  //  Returns YES if the intersection point was found, and
  //  stores that point in X,Y. Returns NO if there is no
  //  determinable intersection point, in which case X,Y will
  //  be unmodified.

  double distAB, theCos, theSin, newX, ABpos;
  // double X, Y;

  //  Fail if either line segment is zero-length.
  // if (Ax==Bx && Ay==By || Cx==Dx && Cy==Dy) return NO;

  //  Count as a clash if the segments share an end-point.
  if (((Ax == Cx) && (Ay == Cy)) || ((Bx == Cx) && (By == Cy)) ||
      ((Ax == Dx) && (Ay == Dy)) || ((Bx == Dx) && (By == Dy))) {
    return 1;
  }

  //  (1) Translate the system so that point A is on the
  //  origin.
  Bx -= Ax;
  By -= Ay;
  Cx -= Ax;
  Cy -= Ay;
  Dx -= Ax;
  Dy -= Ay;

  //  Discover the length of segment A-B.
  distAB = sqrt(Bx * Bx + By * By);

  //  (2) Rotate the system so that point B is on the positive
  //  X axis.
  theCos = Bx / distAB;
  theSin = By / distAB;
  newX = Cx * theCos + Cy * theSin;
  Cy = Cy * theCos - Cx * theSin;
  Cx = newX;
  newX = Dx * theCos + Dy * theSin;
  Dy = Dy * theCos - Dx * theSin;
  Dx = newX;

  //  Fail if segment C-D doesn't cross line A-B.
  if (((Cy < 0.) && (Dy < 0.)) || ((Cy >= 0.) && (Dy >= 0.)))
    return 0;

  //  (3) Discover the position of the intersection point
  //  along line A-B.
  ABpos = Dx + (Cx - Dx) * Dy / (Dy - Cy);

  //  Fail if segment C-D crosses line A-B outside of segment
  //  A-B.
  if (ABpos < 0. || ABpos > distAB)
    return 0;

  //  (4) Apply the discovered position to line A-B in the
  //  original coordinate system.
  // X=Ax+ABpos*theCos;
  // Y=Ay+ABpos*theSin;

  //  Success (i.e. clash).
  return 1;
}

void dump_clash(
    const Shape& shape_a,
    const Shape& shape_b,
    const Vect2& coords_a,
    const Vect2& coords_b,
    double angle_a,
    double angle_b,
    bool flip_a,
    bool flip_b,
    double central_distsq) {
  if (shape_a.resolution() != shape_b.resolution()) {
    throw "Resolutions of shapes are not equal";
  }
  const int resolution = shape_a.resolution();
  double angular_step = 2.0 * PI / resolution;

  double theta_ia0 = 0.0;
  double theta_ia1;
  std::vector<double> theta_tab(resolution / 2);
  std::vector<double> xb_tab(resolution / 2);
  std::vector<double> yb_tab(resolution / 2);
  double xa0, xa1, ya0, ya1;
  double centraldist;
  int angle_a_int, angle_b_int;

  angle_a_int = (int)(angle_a / angular_step + 0.5);
  angle_b_int = (int)(angle_b / angular_step + 0.5);

  std::ofstream fp_test;
  fp_test.open("testclash.xy", std::ios::out);

  // Have a positive then a negative number using power pow(-1, 0) = 1,
  // pow(-1, 1) = 1
  for (int scan_sign_power = 0; scan_sign_power < 2; scan_sign_power++) {
    int scan_sign = pow(-1, scan_sign_power);
    for (int ib = -1; ib < (resolution / 4) + 1; ib++) {
      theta_tab[ib + 1] = fabs(
          (angle_b_int - (1 - 2 * flip_b) * scan_sign * ib) * 2.0 * M_PI / resolution -
          angle_b);
      // sine_theta_ib_tab[ib] = sin(theta_tab[ib]);
      xb_tab[ib + 1] =
          shape_b.radial_points
              [(angle_b_int - (1 - 2 * flip_b) * scan_sign * ib + 5 * resolution) %
               resolution] *
          cos(theta_tab[ib + 1]);
      yb_tab[ib + 1] =
          shape_b.radial_points
              [(angle_b_int - (1 - 2 * flip_b) * scan_sign * ib + 5 * resolution) %
               resolution] *
          sin(theta_tab[ib + 1]);
      if (ib == -1)
        yb_tab[ib + 1] *= -1.0; /* to correct an issue with the fabs of
                                   the very first point */

      fp_test << xb_tab[ib + 1] << "\t" << yb_tab[ib + 1] << std::endl;
    }

    fp_test << 0.0 << "\t" << 0.0 << std::endl;

    theta_ia1 = fabs((angle_a_int + (1 - 2 * flip_a) * scan_sign * -1) * 2.0 * angle_a);
    xa1 = sqrt(central_distsq) -
          shape_a.radial_points
                  [(angle_a_int + (1 - 2 * flip_a) * scan_sign * -1 + 5 * resolution) %
                   resolution] *
              cos(theta_ia0);
    ya1 = shape_a.radial_points
              [(angle_a_int + (1 - 2 * flip_a) * scan_sign * -1 + 5 * resolution) %
               resolution] *
          sin(theta_ia0);
    ya1 *= -1.0; /* to correct an issue with the fabs of the very first
                    point */
    for (int ia = -1; ia < (resolution / 4) + 1; ia++) {

      theta_ia0 = theta_ia1;
      theta_ia1 = fabs(
          (angle_a_int + (1 - 2 * flip_a) * scan_sign * (ia + 1)) * 2.0 * M_PI /
              resolution -
          angle_a);
      ya0 = ya1;
      xa0 = xa1;
      xa1 = centraldist -
            shape_a.radial_points
                    [(angle_a_int + (1 - 2 * flip_a) * scan_sign * (ia + 1) +
                      5 * resolution) %
                     resolution] *
                cos(theta_ia1);
      ya1 = shape_a.radial_points
                [(angle_a_int + (1 - 2 * flip_a) * scan_sign * (ia + 1) +
                  5 * resolution) %
                 resolution] *
            sin(theta_ia1);

      fp_test << xa0 << "\t" << ya0 << std::endl;
    }
  }
  fp_test.close();
}

bool indirect_clash(
    Shape& shape_a,
    Shape& shape_b,
    const Vect2& coords_a,
    const Vect2& coords_b,
    double angle_a,
    double angle_b,
    bool flip_a,
    bool flip_b,
    double central_distsq,
    int monte_carlo_steps) {

  if (shape_a.resolution() != shape_b.resolution()) {
    throw "Resolutions of shapes are not equal";
  }
  const int resolution = shape_a.resolution();
  double angular_step = 2.0 * PI / resolution;
  int iflip_a = flip_a ? -1 : 1;
  int iflip_b = flip_b ? -1 : 1;

  double theta_ia0 = 0.0, theta_ia1;
  std::vector<double> theta_tab(resolution / 2);
  std::vector<double> xb_tab(resolution / 2);
  std::vector<double> yb_tab(resolution / 2);
  double xa0, xa1, ya0, ya1;
  double centraldist;
  int compare_index_a, compare_index_b;

  int angle_a_int, angle_b_int;

  /* this is the most expensive routine of all - try to optimize (especially
   * inside the double loops) */

  /* TSH THIS DOES NOT UTILIZE THE FLIP_A and FLIP_B parameters
     .... it sort of does now, not sure if it's dealing with the round offs
     correctly */

  centraldist = sqrt(central_distsq);

  angle_a_int = (int)(angle_a / angular_step + 0.5);
  angle_b_int = (int)(angle_b / angular_step + 0.5);

  for (int scan_sign_power = 0; scan_sign_power < 2; scan_sign_power++) {
    int scan_sign = pow(-1, scan_sign_power);

    if (monte_carlo_steps == 0) {
      // If we have just initialised (before doing any moves), and are
      // checking for clashes, we want to check all the way around the
      // edge of each shape, to make sure that by chance the two shapes
      // weren't positioned on top of each other.
      for (int ib = -1; ib < resolution / 2; ib++) {
        compare_index_b =
            positive_modulo(angle_b_int + iflip_b * scan_sign * ib, resolution);
        theta_tab[ib + 1] = theta_tab[ib + 1] =
            fabs(compare_index_b * angular_step - angle_b);
        // sine_theta_ib_tab[ib] = sin(theta_tab[ib]);
        xb_tab[ib + 1] = shape_b.get_point(compare_index_b) * cos(theta_tab[ib + 1]);
        yb_tab[ib + 1] = shape_b.get_point(compare_index_b) * sin(theta_tab[ib + 1]);
        if (ib == -1)
          yb_tab[ib + 1] *= -1.0; /* to correct an issue with the fabs
                                     of the very first point */
      }

      compare_index_a = angle_a_int + (1 - 2 * flip_a) * scan_sign * -1;
      theta_ia1 = fabs(compare_index_a * angular_step - angle_a);
      xa1 = sqrt(central_distsq) - shape_a.get_point(compare_index_a) * cos(theta_ia0);
      ya1 = shape_a.get_point(compare_index_a) * sin(theta_ia0);
      ya1 *= -1.0; /* to correct an issue with the fabs of the very first
                      point */
      for (int ia = -1; ia < resolution / 2; ia++) {
        theta_ia0 = theta_ia1;
        compare_index_a =
            positive_modulo(angle_a_int + iflip_a * scan_sign * (ia + 1), resolution);
        theta_ia1 = fabs(compare_index_a * angular_step - angle_a);
        xa0 = xa1;
        ya0 = ya1;
        xa1 = centraldist - shape_a.get_point(compare_index_a) * cos(theta_ia1);
        ya1 = shape_a.get_point(compare_index_a) * sin(theta_ia1);

        for (int ib = -1; ib < resolution / 2; ib++) {
          if (segments_cross(
                  xa0,
                  ya0,
                  xa1,
                  ya1,
                  xb_tab[ib + 1],
                  yb_tab[ib + 1],
                  xb_tab[ib + 2],
                  yb_tab[ib + 2])) {
            // printf("segments cross 1: ib %d\n",ib);
            return 1;
          }
        }
      }
    } else {
      for (int ib = -1; ib < (resolution / 4) + 1; ib++) {
        compare_index_b = angle_b_int + iflip_b * scan_sign * ib;
        theta_tab[ib + 1] = fabs(compare_index_b * angular_step * angle_b);
        // sine_theta_ib_tab[ib] = sin(theta_tab[ib]);
        xb_tab[ib + 1] = shape_b.get_point(compare_index_b) * cos(theta_tab[ib + 1]);
        yb_tab[ib + 1] = shape_b.get_point(compare_index_b) * sin(theta_tab[ib + 1]);
        if (ib == -1)
          yb_tab[ib + 1] *= -1.0; /* to correct an issue with the fabs
                                     of the very first point */
      }

      compare_index_a =
          positive_modulo(angle_a_int + iflip_a * scan_sign * -1, resolution);
      theta_ia1 = fabs(compare_index_a * angular_step - angle_a);
      xa1 = sqrt(central_distsq) - shape_a.get_point(compare_index_a) * cos(theta_ia0);
      ya1 = shape_a.get_point(compare_index_a) * sin(theta_ia0);
      ya1 *= -1.0; /* to correct an issue with the fabs of the very first
                      point */
      for (int ia = -1; ia < (resolution / 4) + 1; ia++) {
        compare_index_a =
            positive_modulo(angle_a_int + iflip_a * scan_sign * (ia + 1), resolution);
        theta_ia0 = theta_ia1;
        theta_ia1 = fabs(compare_index_a * angular_step - angle_a);
        xa0 = xa1;
        ya0 = ya1;
        xa1 = centraldist - shape_a.get_point(compare_index_a) * cos(theta_ia1);
        ya1 = shape_a.get_point(compare_index_a) * sin(theta_ia1);

        for (int ib = -1; ib < (resolution / 4) + 1; ib++) {
          if (segments_cross(
                  xa0,
                  ya0,
                  xa1,
                  ya1,
                  xb_tab[ib + 1],
                  yb_tab[ib + 1],
                  xb_tab[ib + 2],
                  yb_tab[ib + 2])) {
            // printf("segments cross 2: ib %d %f %f %f %f %f %f %f
            // %f\n",ib, xa0, ya0, xa1, ya1, xb_tab[ib+1], yb_tab[ib+1],
            // xb_tab[ib+2], yb_tab[ib+2]);

            // dump_clash(shape_a, shape_b, coords_a, coords_b, angle_a,
            // angle_b, flip_a, flip_b, central_distsq);
            return 1;
          }
        }
      }
    }
  }

  return 0;
}

bool pair_clash(
    Shape& shape_a,
    Shape& shape_b,
    const Vect2& coords_a,
    const Vect2& coords_b,
    double orient_a,
    double orient_b,
    int rotation_a,
    int rotation_b,
    bool flip_a,
    bool flip_b,
    int monte_carlo_steps) {

  if (shape_a.resolution() != shape_b.resolution()) {
    throw "Resolutions of shapes are not equal";
  }
  const int resolution = shape_a.resolution();
  double angular_step = 2.0 * PI / resolution;
  double maxshapecontact_distsq;
  double a_to_b_incline, b_to_a_incline;

  double mathtol = 1e-8;

  double central_distsq = (coords_a - coords_b).norm_sq();

  maxshapecontact_distsq = shape_a.max_radius + shape_b.max_radius;
  maxshapecontact_distsq *= maxshapecontact_distsq;

  /* No clash when urther apart than the maximum shape radii measures */
  if (central_distsq > maxshapecontact_distsq) {
    return false;
  }

  /* we need a polygon overlap calculation */
  /* first calculate the direct overlap */
  a_to_b_incline = acos((coords_b.x - coords_a.x) / sqrt(central_distsq));

  if (isnan(a_to_b_incline)) {
    // printf("nan atob incline 1: %f %f %f =>>> acos(%f)\n", coords_b[0],
    // coords_a[0],
    // central_distsq,(coords_b[0]-coords_a[0])/sqrt(central_distsq));
    if (is_close(coords_b.x - coords_a.x, sqrt(central_distsq), mathtol)) {
      a_to_b_incline = 0.0;
    } else if (is_close(coords_b.x - coords_a.x, -sqrt(central_distsq), mathtol)) {
      a_to_b_incline = M_PI;
    }
  }
  if (coords_b.x < coords_a.x) {
    a_to_b_incline = 2.0 * M_PI - a_to_b_incline;
  }
  b_to_a_incline = a_to_b_incline + M_PI;
  if (flip_a) {
    a_to_b_incline = 2 * M_PI - a_to_b_incline;
  }
  if (flip_b) {
    b_to_a_incline = 2 * M_PI - b_to_a_incline;
  }

  /* now add in the rotation due to the orientation parameters */
  /* This is unaffected by the flip states */
  a_to_b_incline += orient_a;
  b_to_a_incline += orient_b;

  /* now add in the rotation due to the rotation of this image wrt the other
   * images of the same wyckoff */
  a_to_b_incline += (1 - 2 * flip_a) * rotation_a * angular_step;
  b_to_a_incline += (1 - 2 * flip_b) * rotation_b * angular_step;

  a_to_b_incline = positive_modulo(a_to_b_incline, 2.0 * PI);
  b_to_a_incline = positive_modulo(b_to_a_incline, 2.0 * PI);

  return indirect_clash(
      shape_a,
      shape_b,
      coords_a,
      coords_b,
      a_to_b_incline,
      b_to_a_incline,
      flip_a,
      flip_b,
      central_distsq,
      monte_carlo_steps);
}

bool clash_polygon(
    Shape& shape_a,
    Shape& shape_b,
    WallpaperGroup& group,
    Cell& cell,
    std::vector<Site>& occupied_sites,
    Int2& rep_a,
    Int2& rep_b,
    int monte_carlo_steps) {

  ImageType& image_a = group.get_image(occupied_sites[rep_a.x]).image[rep_a.y];
  ImageType& image_b = group.get_image(occupied_sites[rep_b.x]).image[rep_b.y];

  Vect2 fcoords_a = fractcoords(
      image_a.x_coeffs, image_a.y_coeffs, occupied_sites[rep_a.x].site_variables());
  Vect2 fcoords_b = fractcoords(
      image_b.x_coeffs, image_b.y_coeffs, occupied_sites[rep_b.x].site_variables());

  double orient_a = occupied_sites[rep_a.x].angle;
  double orient_b = occupied_sites[rep_b.x].angle;

  /* This means extra rotations due to the symmetries that generated this
   * image from the site basis, done as an integer offset */
  int rotation_a = image_a.rotation_offset;
  int rotation_b = image_b.rotation_offset;

  bool flip_a = image_a.flipped ^ occupied_sites[rep_a.x].flip_site;
  bool flip_b = image_b.flipped ^ occupied_sites[rep_b.x].flip_site;

  Vect2 img_fcoords_b(0, 0);
  /* This method og only looking at the 9 nearest cells fails for those with
   * extreme angles, so instead expand this to 5x5.
   */
  for (int cell_img_x = -2; cell_img_x <= 2; cell_img_x++) {

    /* a is fixed, copies of b are made to test for the clash */
    img_fcoords_b.x = fcoords_b.x + cell_img_x;

    for (int cell_img_y = -2; cell_img_y <= 2; cell_img_y++) {

      if ((rep_a == rep_b) && (cell_img_x == 0) && (cell_img_y == 0)) {
        continue;
      }

      img_fcoords_b.y = fcoords_b.y + cell_img_y;

      Vect2 coords_a = cell.fractional_to_real(fcoords_a);
      Vect2 coords_b = cell.fractional_to_real(fcoords_b);

      if (pair_clash(
              shape_a,
              shape_b,
              coords_a,
              coords_b,
              orient_a,
              orient_b,
              rotation_a,
              rotation_b,
              flip_a,
              flip_b,
              monte_carlo_steps)) {
        return 1;
      }
    }

    return 0;
  }
  return 0;
}

double calculate_packing_fraction(
    Shape& shape,
    WallpaperGroup& group,
    Cell& cell,
    std::vector<Site>& occupied_sites,
    size_t monte_carlo_steps,
    size_t max_steps,
    size_t max_iters = 20) {
  std::vector<Int2> replicaindex;
  int m, r, s;
  bool currclash = 1;
  int iter = 0;

  auto console = spdlog::stdout_color_mt("console");

  int countreplicas = 0;
  for (const Site& site : occupied_sites) {
    const WyckoffType& wyckoff_site = group.get_wyckoff(site);
    for (int m = 0; m < wyckoff_site.multiplicity; m++) {
      replicaindex[countreplicas + m].x = site.wyckoff_index;
      replicaindex[countreplicas + m].y = m;
    }
    countreplicas += wyckoff_site.multiplicity;
  }

  while ((currclash) && (iter < max_iters)) {
    iter++;
    currclash = 0;

    for (int r = 0; r < countreplicas; r++) {
      // this line is new on 2012-12-03. It means we only check one replica of
      // a shape against the other replicas (copies in the same cell, either
      // from the multiplicity of the same occupied site or at other occupied
      // sites and their multiplicities) when that replica is the 'first'
      // (m=0) of its occupied site. If there is no clash from one replica to
      // the others, then there won't be a clash between replicas on the same
      // occupied site.
      if (replicaindex[r].y == 0) {
        // if (1) {
        for (s = 0; s < countreplicas; s++) {
          // s begins at double_count, which starts at 0. this is to avoid
          // doubling up when checking pairs. the first time we take a
          // shape, double_count is 0 so we check it for clashes against
          // every replica. we then add 1 to double_count. we continue
          // (without checking clashes) until we find a replica on a new
          // occupied site, at which point m=0. we then check this replica
          // against all replicas EXCEPT the original replica we checked
          // when double_count was 0, as this pairing has already been
          // checked for clashes.
          if (monte_carlo_steps == 0) {
            currclash = clash_polygon(
                shape,
                shape,
                group,
                cell,
                replicaindex[r],
                replicaindex[s],
                monte_carlo_steps);
            if (currclash) {
              console->warn("found a clash when initialising");
              throw ClashRejection;
            }

          } else {
            currclash = clash_polygon(
                shape,
                shape,
                group,
                cell,
                replicaindex[r],
                replicaindex[s],
                monte_carlo_steps);
            if (currclash) {
              /* if we have found any clash, and we are not in the
               * initial phases (monte_carlo_steps>0), then there's no need to
               * search through for any more. exit the two for loops
               * and the while loop immediately. also, flag there has
               * been a clash so we can
               * reject the move that caused it */
              throw ClashRejection;
            }
          }
        }
      }
    }
  }

  double packing_fraction = countreplicas * shape.area() / cell.area();
  if (isnan(packing_fraction)) {
    printf("nan encountered %f %f %f\n", cell.x_len, cell.y_len, cell.angle);
    exit(1);
  }

  return packing_fraction;
}

int initialize_structure_in_group(
    Shape& shape,
    WallpaperGroup& group,
    Cell& cell,
    std::vector<Site>& occupied_sites,
    std::vector<Basis>& basis) {

  int count_replicas = 0;
  Basis b;

  // Logging to console which can be turned off easily
  auto console = spdlog::stdout_color_mt("console");

  count_replicas += group.group_multiplicity();

  // cell angles.
  b = Basis(&cell.angle);
  b.cell_angle = true;
  if (group.hexagonal) {
    console->debug("Hexagonal group");
    *b.value = PI / 3;
  } else if (group.rectangular) {
    console->debug("Rectangular group");
    *b.value = PI / 2;
  } else {
    console->debug("Tilted group");
    b.fixed = false;
    b.min_val = PI / 4;
    b.max_val = 3 * PI / 4;
    *b.value = b.min_val + fluke() * b.value_range();
  }

  basis.push_back(b);

  // cell sides.
  double max_cell_size = 4 * shape.max_radius * count_replicas;
  if (group.a_b_equal) {
    console->debug("Cell sides equal");
    basis.push_back(Basis(max_cell_size, 0.1, max_cell_size));
    b = Basis(&
    basis[1].cell_side = true;
    cell.x_len = basis[1].value;
    cell.y_len = basis[1].value;
  } else {
    basis.push_back(Basis(max_cell_size, 0.1, max_cell_size));
    basis.push_back(Basis(max_cell_size, 0.1, max_cell_size));
    basis[1].cell_side = true;
    basis[2].cell_side = true;
    cell.x_len = basis[1].value;
    cell.y_len = basis[2].value;
  }

  // now position the particles.
  for (Site& site : occupied_sites) {
    WyckoffType& wyckoff_site = group.get_wyckoff(site);
    count_replicas += wyckoff_site.multiplicity;

    console->debug("Wyckoff site: %c ", wyckoff_site.letter);

    if (fabs(wyckoff_site.image[0].x_coeffs.x) > 0.1) {
      /* x is variable*/
      basis.push_back(Basis(fluke(), 0, 1));
      site.x = basis.back().value;
      console->debug("Site x variable %f\n", site.x);
    }
    if (fabs(wyckoff_site.image[0].y_coeffs.y) > 0.1) {
      /* then y is variable*/
      basis.push_back(Basis(fluke(), 0, 1));
      site.y = basis.back().value;
      console->debug("Site y variable %f\n", site.y);
    }

    /*choose the orientation of the zeroth image of this particle */
    /* This could also be done within the Wyckoff position SITEROTATION
     * parameter?... */
    if (wyckoff_site.site_mirrors) {
      double value = M_PI / 180 * wyckoff_site.image[0].site_mirror;
      basis.push_back(Basis(value, true));
      site.angle = basis.back().value;
    } else {
      double value = fluke() * 2.0 * M_PI;
      basis.push_back(Basis(value, 0, 2 * PI));
      site.angle = basis.back().value;
      console->debug("site offset-angle is variable %f\n", site.angle);
    }
  }

  console->debug("replicas %d variables %d os ", count_replicas, basis.size());

  return basis.size();
}

void uniform_best_packing_in_isopointal_group(
    Shape& shape,
    WallpaperGroup& group,
    size_t num_cycles,
    size_t max_steps,
    double max_step_size,
    double kT_start = 0.1,
    double kT_finish = 5e-4) {
  int vary_index; /* which of our basis values are we currently trying to vary
                   */
  double old_value;

  double max_phi_in_cycle;
  // double old_angle;
  double old_cell0, old_cell1;
  int i, os, m;

  bool flips[MAXNUMOCCSITES];

  auto console = spdlog::stdout_color_mt("console");
  std::vector<double> best_basis;

  Vect2 fcoords; /*fractional coordinates (scaled between 0 and 1) */
  Vect2 coords;

  double image_scale, shape_scale;
  int cell_img_x, cell_img_y; /* indexes of x and y periodic cells used to check
                                 for overlap */
  double image_width, image_height;
  double origin_translate_x, origin_translate_y;

  int flip_site = -1;

  int clash_rejection = 0;

  size_t countreplicas = group.group_multiplicity();

  double packing_fraction = -1;
  double packing_fraction_max = 0.0;

  double kT_ratio = pow(kT_finish / kT_start, 1.0 / max_steps);

  for (size_t cycle = 0; cycle < num_cycles; cycle++) {
    /* Each cycle starts with a new random initialisation */
    size_t monte_carlo_steps = 0;
    size_t rejections = 0;
    double kT = kT_start;
    int flip_index;

    Basis basis_prev;
    Site site_prev;
    Cell cell_prev;
    double packing_fraction_prev;

    std::vector<Basis> basis;
    std::vector<Site> occupied_sites;
    Cell cell{};
    initialize_structure_in_group(shape, group, cell, occupied_sites, basis);

    packing_fraction = calculate_packing_fraction(
        shape, group, cell, occupied_sites, monte_carlo_steps, max_steps, 20);

    console->info("Initial packing fraction = %f\n", packing_fraction);

    double packing_fraction_max_cycle = 0.0;

    while (monte_carlo_steps < max_steps) {
      kT *= kT_ratio;
      clash_rejection = 0;

      /* Search until a variable basis is found */
      do {
        vary_index = rand() % basis.size();
      } while (!basis[vary_index].fixed);
      Basis& basis_current = basis[vary_index];

      /* Occasionally allow flips */
      if (monte_carlo_steps % 100) {
        flip_index = rand() % occupied_sites.size();
        flips[flip_index] ^= 1;
      } else {
        flip_index = -1;
      }

      site_prev = Site(occupied_sites[flip_index]);
      cell_prev = Cell(cell);
      basis_prev = Basis(basis_current);
      packing_fraction_prev = packing_fraction;

      if (basis_current.fixed) {
        /* This is a regular structural parameter (non cell-side) which can be
         * varied within its given range */
        if (!basis_current.cell_side) {
          basis_current.value +=
              max_step_size * basis_current.value_range() * (fluke() - 0.5);
        } else {
          /* This is a cell side.  Since these have a direct effect on the
           * packing fraction, make their step size in proportion to kT */
          basis_current.value *= (1.0 + std::min(3.0 * kT, 0.1) * (fluke() - 0.5));
        }
        basis_current.validate();

        if (vary_index == 0) {
          /* this is the cell angle. stretch cell sides so that the total area
           * is preserved (to prevent runaway angle collapse) */
          Cell.x_len *= sqrt(sin(basis_prev.value) / sin(basis_current.value));
          Cell.y_len *= sqrt(sin(basis_prev.value) / sin(basis_current.value));
        }

      } else {
        /* basis_fixed[vary_index]==2 which means it's a particle orientation
         * which must vary in a quantized manner to ensure site mirrors still
         * lie on symmetry mirror planes */

        if ((shape->mirrors % 2 == 0) && (Fluke() < 0.5)) {
          /* turn it 90 degrees to switch the x and y mirror planes */
          if (basis[vary_index] < M_PI * 3.0 / 4.0)
            basis[vary_index] += M_PI / (double)shape->mirrors;
          else
            basis[vary_index] -= M_PI / (double)shape->mirrors;
        } else {
          /* turn it 180 degrees so that all mirror planes are preserved (just
           * reversed) */
          if (basis[vary_index] < M_PI)
            basis[vary_index] += M_PI;
          else
            basis[vary_index] -= M_PI;
        }
      }

      tests++; /* increase tests here so that when calculate packing fraction on
                  first test, tests=1 */

      phi = calculate_packing_fraction(
          shape,
          group,
          cellsides,
          cellangles,
          occupiedsites,
          sitevariables,
          flips,
          tests,
          &clash_rejection);

      if ((clash_rejection == 1) ||
          (Fluke() > exp(((1.0 / old_phi - 1.0 / phi) / kT) +
                         countreplicas * log(old_phi / phi)))) {
        // if (phi<old_phi) {
        /* some probability of rejecting move if packing fraction decreased, or
         * definitely reject it if it caused a clash */
        rejections++;
        basis[vary_index] = old_value;
        *cellsides[0] = old_cell0;
        *cellsides[1] = old_cell1;
        if (flip_site != -1)
          flips[flip_site] ^= 1;
        phi = old_phi;
        /*if (clash_rejection!=1) {
          printf("######not a clash rejection, return to old values\n");
        }
        if (clash_rejection==1) {
          printf("############################a clash rejection, return to old
        values\n");
          }*/
      }

      if (phi > max_phi) {
        /* best packing seen yet ... save data */
        for (i = 0; i < basis_size; i++) {
          best_basis[i] = basis[i];
        }
        for (i = 0; i < numoccsites; i++) {
          best_flips[i] = flips[i];
        }
        max_phi = phi;
      }

      if (phi > max_phi_in_cycle) {
        /* best packing seen yet in current cycle ... save data */
        for (i = 0; i < basis_size; i++) {
          best_basis_in_cycle[i][cycle] = basis[i];
        }
        for (i = 0; i < numoccsites; i++) {
          best_flips_in_cycle[i][cycle] = flips[i];
        }
        max_phi_in_cycle = phi;
      }

      if (!EMPTY_OUTPUT) {
        if (tests % 500 == 0)
          printf(
              "cycle %d of %d, step %ld of %d, kT=%g, packing %f, angle %f, "
              "b/a=%f, rejection %f percent\n",
              cycle + 1,
              CYCLES,
              monte_carlo_steps,
              MAXSTEPS,
              kT,
              phi,
              *cellangles[0] * 180.0 / M_PI,
              (*cellsides[1]) / (*cellsides[0]),
              (100.0 * rejections) / monte_carlo_steps);
      }
    }

    monte_carlo_steps = MAXSTEPS + 1;
    /* add 1 to monte_carlo_steps so that monte_carlo_steps=MAXSTEPS+1, to signify to
     * packing_fraction function not to check for clashes when reloading best structures
     */

    /*now reload the best structure ever seen in the current cycle*/
    for (i = 0; i < basis_size; i++) {
      basis[i] = best_basis_in_cycle[i][cycle];
    }
    for (i = 0; i < numoccsites; i++) {
      flips[i] = best_flips_in_cycle[i][cycle];
    }

    phi = calculate_packing_fraction(
        group,
        cellsides,
        cellangles,
        occupiedsites,
        sitevariables,
        flips,
        monte_carlo_steps,
        &clash_rejection);
    if (!EMPTY_OUTPUT) {
      printf(
          "BEST: cell %f %f angle %6.2f packing %f rejection (%f\%%) ",
          *cellsides[0],
          *cellsides[1],
          *cellangles[0] * 180.0 / M_PI,
          max_phi_in_cycle,
          (100.0 * rejections) / monte_carlo_steps);
      for (i = 0; i < basis_size; i++) {
        printf("%f ", basis[i]);
      }
      for (i = 0; i < numoccsites; i++) {
        printf("%d ", flips[i]);
      }
      printf("\n");
    }
#define btoa(x) ((x) ? "true" : "false")
    char chiral_state;
    int chiralsum = 0;
    int totalsum = 0;
    for (os = 0; os < numoccsites; os++) {
      chiralsum +=
          (best_flips[os] * 2 - 1) * group.wyckoffs[occupiedsites[os]].multiplicity;
      totalsum += abs(
          (best_flips[os] * 2 - 1) * group.wyckoffs[occupiedsites[os]].multiplicity);
    }
    if (chiralsum != 0) {
      chiral_state = 'c';
      if (chiralsum != totalsum) {
        chiral_state = 's';
      }
    } else {
      chiral_state = 'a';
    }

    fprintf(
        fpout,
        "%8s\t%s\t%d\t%.6f\t%.7f\t%.3f\t%s\t%d\t%.5f\t%d\t%.3f\t%d\t%.4f\t%"
        ".4f\t%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%c",
        group.label,
        oslist,
        cycle + 1,
        max_phi,
        area(shape),
        (100.0 * rejections) / monte_carlo_steps,
        shape->name,
        MAXSTEPS,
        max_step_size,
        SHAPE_RESOLUTION,
        shape_var,
        numoccsites,
        shape->minr,
        shape->maxr,
        shape->minr / shape->maxr,
        max_phi_in_cycle,
        *cellsides[0],
        *cellsides[1],
        *cellangles[0] * 180.0 / M_PI,
        chiral_state);

    for (os = 0; os < numoccsites; os++) {
      fprintf(
          fpout,
          "\t%s\t%d",
          btoa(best_flips[os]),
          group.wyckoffs[occupiedsites[os]].multiplicity);
    }
    fprintf(fpout, "\n");
    fprintf(
        fp_isopointal_out,
        "%.3f\t%.6f\t%.7f\t%.3f\t%s\t%8s\t%s\t%d\t%d\t%.5f\t%d\t%d\t%.4f\t%"
        ".4f\t%.4f\t%.6f\t%.6f\t%.6f\t%.6f\n",
        shape_var,
        max_phi,
        area(shape),
        (100.0 * rejections) / monte_carlo_steps,
        shape->name,
        group.label,
        oslist,
        cycle + 1,
        MAXSTEPS,
        max_step_size,
        SHAPE_RESOLUTION,
        numoccsites,
        shape->minr,
        shape->maxr,
        shape->minr / shape->maxr,
        max_phi_in_cycle,
        *cellsides[0],
        *cellsides[1],
        *cellangles[0] * 180.0 / M_PI);

    fflush(NULL);

    /*now reload the best structure ever seen, for the image*/
    for (i = 0; i < basis_size; i++) {
      basis[i] = best_basis[i];
    }
    for (i = 0; i < numoccsites; i++) {
      flips[i] = best_flips[i];
    }
    phi = calculate_packing_fraction(
        shape,
        group,
        cellsides,
        cellangles,
        occupiedsites,
        sitevariables,
        flips,
        monte_carlo_steps,
        &clash_rejection);

    image_scale = (MAXGRAPHICSIZE - 2 * GRAPHICFRAME) /
                  (3.0 * MAX(*cellsides[0] + fabs(*cellsides[1] * cos(*cellangles[0])),
                             *cellsides[1] * sin(*cellangles[0])));
    shape_scale = 1.0; /* (image_scale-1.0)/image_scale;  this was originally an attempt
                          to correct for the 1px stroke width, but since that applies to
                          all directions additively, it is better to remove it in the
                          definition of the shape itself rather than try to scale it */
    origin_translate_x =
        GRAPHICFRAME + MAX(0.0, -image_scale * *cellsides[1] * cos(*cellangles[0]));
    origin_translate_y = GRAPHICFRAME;
    image_width = 2 * GRAPHICFRAME +
                  (MAXGRAPHICSIZE - 2 * GRAPHICFRAME) *
                      (*cellsides[0] + fabs(*cellsides[1] * cos(*cellangles[0]))) /
                      MAX(*cellsides[0] + fabs(*cellsides[1] * cos(*cellangles[0])),
                          *cellsides[1] * sin(*cellangles[0]));
    image_height = 2 * GRAPHICFRAME +
                   (MAXGRAPHICSIZE - 2 * GRAPHICFRAME) *
                       (*cellsides[1] * sin(*cellangles[0])) /
                       MAX(*cellsides[0] + fabs(*cellsides[1] * cos(*cellangles[0])),
                           *cellsides[1] * sin(*cellangles[0]));

    fp = fopen(filename, "w");
    fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
    fprintf(
        fp,
        "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
        "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
    fprintf(
        fp,
        "<svg width=\"%f\" height=\"%f\" version=\"1.1\" "
        "xmlns=\"http://www.w3.org/2000/svg\" "
        "xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n",
        image_width,
        image_height);
    fprintf(fp, "<defs>\n");
    fprintf(fp, "<g id=\"cell\" opacity=\"1.0\" >\n");
    fprintf(
        fp,
        "<polygon points=\"0,0 %f,0 %f,%f %f,%f\" "
        "style=\"stroke:#000000;stroke-width:2\"/>\n",
        image_scale * *cellsides[0],
        image_scale * *cellsides[0] + image_scale * *cellsides[1] * cos(*cellangles[0]),
        image_scale * *cellsides[1] * sin(*cellangles[0]),
        image_scale * *cellsides[1] * cos(*cellangles[0]),
        image_scale * *cellsides[1] * sin(*cellangles[0]));
    fprintf(fp, "</g>\n");
    fprintf(fp, "<g id=\"shape\" opacity=\"0.7\" >\n");
    if (POLYGON_REP == 1) {
      fprintf(fp, "<polygon points=\"");
      for (i = 0; i < SHAPE_RESOLUTION; i++) {
        fprintf(
            fp,
            "%f,%f ",
            (image_scale * shape->r[i] - 0.5) * cos(i * 2.0 * M_PI / SHAPE_RESOLUTION),
            (image_scale * shape->r[i] - 0.5) * sin(i * 2.0 * M_PI / SHAPE_RESOLUTION));
      }
      fprintf(fp, "\" style=\"stroke:#000000;stroke-width:1\"/>\n");
    } else if (POLYGON_REP == 2) {
      fprintf(fp, "<polygon points=\"");
      for (i = 0; i < SHAPE_RESOLUTION; i++) {
        fprintf(
            fp,
            "%f,%f ",
            image_scale * shape->r[i] * cos(i * 2.0 * M_PI / SHAPE_RESOLUTION),
            image_scale * shape->r[i] * sin(i * 2.0 * M_PI / SHAPE_RESOLUTION));
      }
      fprintf(fp, "\" style=\"stroke:#000000;stroke-width:0\"/>\n");
    } else {
      /* draw radial lines */
      for (i = 0; i < SHAPE_RESOLUTION; i++) {
        fprintf(
            fp,
            "<line x1=\"0\" y1=\"0\" x2=\"%f\" y2=\"%f\" stroke=\"purple\" "
            "stroke-width=\"1\" />\n",
            image_scale * shape->r[i] * cos(i * 2.0 * M_PI / SHAPE_RESOLUTION),
            image_scale * shape->r[i] * sin(i * 2.0 * M_PI / SHAPE_RESOLUTION));
      }
    }
    if (AXES_SHOWN) {
      fprintf(
          fp,
          "<line x1=\"0\" y1=\"0\" x2=\"20\" y2=\"0\"  "
          "stroke=\"black\" stroke-width=\"1\"  />\n");
      fprintf(
          fp,
          "<line x1=\"0\" y1=\"0\" x2=\"0\" y2=\"20\"  stroke=\"red\" "
          "stroke-width=\"1\"  />\n");
    }
    if (CENTRE_SPOT) {
      fprintf(
          fp,
          "<circle cx=\"0.0\" cy=\"0.0\" r=\"3.0\" fill=\"black\"  "
          "stroke=\"black\" stroke-width=\"0\" />\n");
    }
    fprintf(fp, "</g>\n");
    fprintf(fp, "</defs>\n");
    for (cell_img_x = 0; cell_img_x < 3; cell_img_x++) {
      for (cell_img_y = 0; cell_img_y < 3; cell_img_y++) {
        fcoords[0] = 1.0 * cell_img_x;
        fcoords[1] = 1.0 * cell_img_y;
        realcoords(fcoords, coords, cellsides, cellangles);
        if (cell_img_x == 1 && cell_img_y == 1)
          fprintf(
              fp,
              "<use xlink:href=\"#cell\" transform=\"translate(%f,%f) "
              "scale(1)\" style=\"fill:grey\"/>\n",
              origin_translate_x + image_scale * coords[0],
              origin_translate_y + image_scale * coords[1]);
        else
          fprintf(
              fp,
              "<use xlink:href=\"#cell\" transform=\"translate(%f,%f) "
              "scale(1)\" style=\"fill:white\"/>\n",
              origin_translate_x + image_scale * coords[0],
              origin_translate_y + image_scale * coords[1]);
      }
    }
    for (os = 0; os < numoccsites; os++) {
      for (m = 0; m < group.wyckoffs[occupiedsites[os]].multiplicity; m++) {
        fractcoords(
            group.wyckoffs[occupiedsites[os]].image[m].coord_coeffs,
            sitevariables[os],
            fcoords);
        for (cell_img_x = 0; cell_img_x < 3; cell_img_x++) {
          for (cell_img_y = 0; cell_img_y < 3; cell_img_y++) {
            fcoords[0] += 1.0 * cell_img_x;
            fcoords[1] += 1.0 * cell_img_y;
            realcoords(fcoords, coords, cellsides, cellangles);
            fcoords[0] -= 1.0 * cell_img_x;
            fcoords[1] -= 1.0 * cell_img_y;
            fprintf(
                fp,
                "<use xlink:href=\"#shape\" transform=\"translate(%f,%f) "
                "rotate(%f, 0, 0) scale(%f %f)\" style=\"fill:",
                origin_translate_x + image_scale * coords[0],
                origin_translate_y + image_scale * coords[1],
                -((group.wyckoffs[occupiedsites[os]].image[m].flipped ^ flips[os]) *
                      (-2) +
                  1) * (*sitevariables[os][2] * 180.0 / M_PI) -
                    group.wyckoffs[occupiedsites[os]].image[m].rotation_offset * 360.0 /
                        SHAPE_RESOLUTION,
                shape_scale,
                ((group.wyckoffs[occupiedsites[os]].image[m].flipped ^ flips[os]) *
                     (-2) +
                 1) *
                    shape_scale);
            if (cell_img_x == 1 && cell_img_y == 1) {
              fprintf(
                  fp,
                  "rgb(%d,%d,%d)",
                  (int)(0.1 * 256 / (os + 1)),
                  (int)(0.3 * 256 / (os + 1)),
                  (int)(0.6 * 256 / (os + 1)));
            } else {
              fprintf(
                  fp,
                  "rgb(%d,%d,%d)",
                  (int)(0.1 * 256 / (os + 1)),
                  (int)(0.6 * 256 / (os + 1)),
                  (int)(0.3 * 256 / (os + 1)));
            }
            fprintf(fp, "\"/>\n");
          }
        }
      }
    }
    fprintf(fp, "</svg>\n");
    fprintf(fp, "<!--Copyright: Toby Hudson\n");
    for (i = 0; i < basis_size; i++) {
      fprintf(fp, "%f ", basis[i]);
    }
    for (i = 0; i < numoccsites; i++) {
      fprintf(fp, "%d ", flips[i]);
    }
    fprintf(
        fp,
        "\ncell %f %f angle %6.2f packing %f (max %f) rejection (%f\%%) ",
        *cellsides[0],
        *cellsides[1],
        *cellangles[0] * 180.0 / M_PI,
        phi,
        max_phi,
        (100.0 * rejections) / monte_carlo_steps);
    fprintf(fp, "-->\n");
    fclose(fp);
  }
}
}
}
}
