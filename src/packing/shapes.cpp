#include "shapes.h"

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

Vect2 Cell::fractional_to_real(const Vect2& fractional) {
  Vect2 v(0, 0);
  v.x = fractional.x * this->x_len + fractional.y * this->y_len * cos(this->angle);
  v.y = fractional.y * this->y_len * sin(this->angle);
  return v;
}

int WallpaperGroup::group_multiplicity() const {
  int total_multiplicity = 0;
  for (int site_index : occupied_sites) {
    total_multiplicity += this->get_occupied_site(site_index).multiplicity;
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
