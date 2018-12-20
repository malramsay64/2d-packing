#include "shapes.h"
#include <cassert>

Shape::Shape(
    const std::string& name,
    const std::vector<double>& radial_points,
    const int rotational_symmetries,
    const int mirrors)
    : name(name), radial_points(radial_points),
      rotational_symmetries(rotational_symmetries), mirrors(mirrors) {
  auto min_max = std::minmax_element(radial_points.begin(), radial_points.end());
  this->min_radius = *min_max.first;
  this->max_radius = *min_max.second;
};

Shape::Shape(const std::string& name, const std::vector<double>& radial_points)
    : Shape(name, radial_points, 1, 0){};

void Shape::plot(const std::string& filename) const {
  double angle;

  std::ofstream outfile;
  outfile.open(filename, std::ios::out);

  for (size_t i = 0; i < this->radial_points.size(); ++i) {
    angle = this->angular_step() * i;
    outfile << std::fixed << std::setprecision(12)
            << this->radial_points[i] * cos(angle) << " "
            << this->radial_points[i] * sin(angle) << std::endl;
  }
  outfile.close();
};

double Shape::area() const {
  /* a calculation of the area of the polygon */
  double areasum = 0.0;
  for (size_t i = 0; i < this->radial_points.size(); ++i) {
    areasum += 0.5 * this->radial_points[i] *
               this->radial_points[(i + 1) % this->resolution()] * this->angular_step();
  }
  return areasum;
}

Vect2 Cell::fractional_to_real(const Vect2& fractional) const {
  Vect2 v(0, 0);
  v.x = fractional.x * this->x_len->get_value() +
        fractional.y * this->y_len->get_value() * cos(this->angle->get_value());
  v.y = fractional.y * this->y_len->get_value() * sin(this->angle->get_value());
  return v;
}

Vect2 ImageType::real_to_fractional(const Vect3& real) const {
  /* converts site variables and wyckoff site coefficients into the location
   * of the actual wyckoff image in fractional coordinates */
  Vect2 v(0, 0);
  v.x = this->x_coeffs.x * real.x + this->x_coeffs.y * real.y + this->x_coeffs.z;
  v.y = this->y_coeffs.x * real.x + this->y_coeffs.y * real.y + this->y_coeffs.z;
  return positive_modulo(v, 1.);
}

Vect2 ImageType::real_to_fractional(const Site& site) const {
  /* converts site variables and wyckoff site coefficients into the location
   * of the actual wyckoff image in fractional coordinates */
  Vect2 v(0, 0);
  v.x = this->x_coeffs.x * site.x->get_value() +
        this->x_coeffs.y * site.y->get_value() + this->x_coeffs.z;
  v.y = this->y_coeffs.x * site.x->get_value() +
        this->y_coeffs.y * site.y->get_value() + this->y_coeffs.z;
  return positive_modulo(v, 1.);
}

size_t group_multiplicity(const std::vector<Site>& occupied_sites) {
  size_t total_multiplicity = 0;
  for (const Site& site : occupied_sites) {
    total_multiplicity = site.wyckoff->multiplicity;
  }
  return total_multiplicity;
}

std::string create_filename(
    const Shape& shape,
    const WallpaperGroup& group,
    const std::string site_list,
    const std::string& directory) {
  std::stringstream stream_filename;
  stream_filename << directory << "/solution_" << shape.name << "_" << group.label
                  << "_" << site_list << "_" << std::fixed << std::setprecision(3)
                  << shape.shape_var << ".svg";
  return stream_filename.str();
}

void Basis::set_value(double new_value) {
  if (new_value < this->min_val) {
    new_value = this->min_val;
  } else if (new_value > this->max_val) {
    new_value = this->max_val;
  }
  this->value_previous = this->value;
  this->value = new_value;
}

double Basis::get_random_value(double kT) const {
  return this->value + this->step_size * this->value_range() * (fluke() - 0.5);
}

double CellLengthBasis::get_random_value(double kT) const {
  return this->get_value() + (1.0 + std::min(3.0 * kT, 0.1) * (fluke() - 0.5));
}

double CellAngleBasis::get_random_value(double kT) const {
  return this->value + this->step_size * this->value_range() * (fluke() - 0.5);
}

void CellAngleBasis::set_value(double new_value) {
  if (new_value < this->min_val) {
    new_value = this->min_val;
  } else if (new_value > this->max_val) {
    new_value = this->max_val;
  }
  this->value_previous = this->value;
  this->value = new_value;
  this->update_cell_lengths();
};

void CellAngleBasis::reset_value() {
  this->value = this->value_previous;
  this->reset_cell_lengths();
};

void CellAngleBasis::update_cell_lengths() {
  this->cell_x_len->set_value(
      this->cell_x_len->get_value() *
      sqrt(sin(this->value_previous) / sin(this->value)));
  this->cell_y_len->set_value(
      this->cell_y_len->get_value() *
      sqrt(sin(this->value_previous) / sin(this->value)));
};

void CellAngleBasis::reset_cell_lengths() {
  this->cell_x_len->reset_value();
  this->cell_y_len->reset_value();
};

double MirrorBasis::get_random_value(double kT) const {
  if ((this->mirrors % 2 == 0) && (fluke() < 0.5)) {
    /* turn it 90 degrees to switch the x and y mirror planes */
    if (this->value < M_PI * 3.0 / 4.0) {
      return this->value + PI / this->mirrors;
    }
    return this->value - PI / this->mirrors;
  }
  /* turn it 180 degrees so that all mirror planes are preserved */
  return positive_modulo(this->value + PI, 2 * PI);
}

double FlipBasis::get_random_value(double kT) const {
  return rand() % occupied_sites->size();
}

void FlipBasis::set_value(double new_value) {
  this->value_previous = static_cast<int>(std::round(new_value));
  occupied_sites->at(this->value_previous).flip_site ^= 1;
}

void FlipBasis::reset_value() {
  if (this->value_previous == -1) {
    return;
  }
  occupied_sites->at(this->value_previous).flip_site ^= 1;
  // Only allow the reset_value to occur once
  this->value_previous = -1;
}

bool ShapeInstance::pair_clash(ShapeInstance& other) const {
  // We can only compare shapes with the same resolution
  if (this->shape->resolution() != other.shape->resolution()) {
    throw "Resolutions of shapes are not equal";
  }

  double max_contact_dist_squared =
      std::pow(this->shape->max_radius + other.shape->max_radius, 2);

  Vect2 coords_this = this->get_real_coordinates();
  Vect2 coords_other = other.get_real_coordinates();

  double central_dist_sqaured = (coords_this - coords_other).norm_sq();
  // Where the distance is greater than both maximum radii, there is definitely no
  // overlap, so we can avoid the expensive checks.
  if (max_contact_dist_squared > central_dist_sqaured) {
    return false;
  }
  return true;
}
