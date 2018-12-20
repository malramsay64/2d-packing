#include "basis.h"
#include "shapes.h"
#include <cassert>
#include <pybind11/stl.h>

namespace py = pybind11;

double Cell::area() const {
  return this->x_len->get_value() * this->y_len->get_value() *
         fabs(sin(this->angle->get_value()));
};

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

int Shape::resolution() const { return this->radial_points.size(); };
double Shape::angular_step() const { return 2 * PI / this->resolution(); };
double Shape::get_point(int index) const { return this->radial_points.at(index); }

void Shape::plot(const std::string& filename) const {
  double angle;

  std::ofstream outfile;
  outfile.open(filename, std::ios::out);

  for (std::size_t i = 0; i < this->radial_points.size(); ++i) {
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
  for (std::size_t i = 0; i < this->radial_points.size(); ++i) {
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

bool ShapeInstance::operator==(const ShapeInstance& other) const {
  return (
      this->shape == other.shape && this->site == other.site &&
      this->image == other.image);
};

Vect2 ShapeInstance::get_fractional_coordinates() const {
  return this->image->real_to_fractional(*this->site);
}
Vect2 ShapeInstance::get_real_coordinates() const {
  return this->site->get_position();
};
double ShapeInstance::get_angle() const { return this->site->angle->get_value(); };
double ShapeInstance::get_rotational_offset() const {
  return this->image->rotation_offset;
};
bool ShapeInstance::get_flipped() const {
  return this->image->flipped ^ this->site->flip_site;
};

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

std::size_t group_multiplicity(const std::vector<Site>& occupied_sites) {
  std::size_t total_multiplicity = 0;
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

void export_Shape(py::module& m) {
  py::class_<Shape> shape(m, "Shape");
  shape.def(py::init<const std::string&, const std::vector<double>&, int, int>())
      .def("resolution", &Shape::resolution)
      .def("plot", &Shape::plot)
      .def("angular_step", &Shape::angular_step)
      .def("get_point", &Shape::get_point);
}
