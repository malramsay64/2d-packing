#include "math.h"
#include "random.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#ifndef SHAPES_H
#define SHAPES_H

// Forward declaration of Site class
class Site;

class Basis {
protected:
  double value_previous;
  double value;

public:
  double min_val;
  double max_val;
  double step_size = 1;

  Basis(double value, double min_val, double max_val, double step_size)
      : value(value), min_val(min_val), max_val(max_val), step_size(step_size){};
  Basis(double value, double min_val, double max_val)
      : Basis(value, min_val, max_val, 1){};

  double value_range() const { return this->max_val - this->min_val; };
  double get_value() const { return this->value; };
  void set_value(double new_value);
  void reset_value() { this->value = this->value_previous; };
  double get_random_value(double kT) const;
};

class CellLengthBasis : public Basis {
public:
  double step_size;

  CellLengthBasis(double value, double min_val, double max_val, double step_size)
      : Basis(value, min_val, max_val), step_size(step_size){};

  double get_random_value(double kT) const;
};

class CellAngleBasis : public Basis {
private:
  void update_cell_lengths();
  void reset_cell_lengths();

public:
  double step_size;
  Basis* cell_x_len;
  Basis* cell_y_len;

  CellAngleBasis(
      double value,
      double min_val,
      double max_val,
      double step_size,
      Basis* cell_x_len,
      Basis* cell_y_len)
      : Basis(value, min_val, max_val), step_size(step_size), cell_x_len(cell_x_len),
        cell_y_len(cell_y_len){};

  void set_value(double new_value);
  void reset_value();
  double get_random_value(double kT) const;
};

class FixedBasis : public Basis {
public:
  FixedBasis(double value) : Basis(value, value, value){};

  void set_value(double new_value){};
  void reset_value(){};
};

class MirrorBasis : public Basis {
public:
  int mirrors;

  MirrorBasis(double value, double min_val, double max_val, int mirrors)
      : Basis(value, min_val, max_val), mirrors(mirrors){};

  double get_random_value(double kT) const;
};

class FlipBasis : public Basis {
private:
  std::vector<Site>* occupied_sites;
  int value_previous;

public:
  FlipBasis(std::vector<Site>* occupied_sites)
      : Basis(0, 0, occupied_sites->size()), occupied_sites(occupied_sites){};

  double get_random_value(double kT) const;
  void set_value(double new_value);
  void reset_value();
};

struct Cell {
  Basis* x_len;
  Basis* y_len;
  Basis* angle;

  Vect2 fractional_to_real(const Vect2&) const;
  double area() const {
    return this->x_len->get_value() * this->y_len->get_value() *
           fabs(sin(this->angle->get_value()));
  };
};

class Shape {
public:
  Shape(
      const std::string& name,
      const std::vector<double>& radial_points,
      const int rotational_symmetries,
      const int mirrors);
  Shape(const std::string& name, const std::vector<double>& radial_points);
  std::string name;
  std::vector<double> radial_points;
  int rotational_symmetries;
  int mirrors;
  double min_radius;
  double max_radius;
  double shape_var = 0;

  inline int resolution() const { return this->radial_points.size(); };
  inline double angular_step() const { return 2 * PI / this->resolution(); };
  inline double get_point(const int index) const { return this->radial_points[index]; }

  void plot(const std::string& filename) const;
  double area() const;
};

/*! \enum mirror
 *
 *  Detailed description
 */
enum mirror {
  m0 = 0,
  m90 = 90,
  m45 = 45,
  m135 = 135,
  m30 = 30,
  m60 = 60,
  m330 = 330,
  m300 = 300
};

class ImageType {
public:
  ImageType(Vect3 x_coeffs, Vect3 y_coeffs, double rotation_offset);
  Vect3 x_coeffs;
  Vect3 y_coeffs;
  /* first index is the output coordinate,
   second index is input coordinate (+1 is for the constant that is added)

   for example, to work out the new values:
   x_new = coord_coeffs[0][0]*x_old + coord_coeffs[0][1]*y_old +
   coord_coeffs[0][2]; y_new = coord_coeffs[1][0]*x_old +
   coord_coeffs[1][1]*y_old + coord_coeffs[1][2];
  */
  double rotation_offset; /*angles move clockwise*/
  bool flipped; /* =1, when flipped using the x-axis as mirror, rotation_offset
                   is then employed if mirror is at y-axis for example */
  enum mirror site_mirror;

  Vect2 real_to_fractional(const Vect3& real) const;
  Vect2 real_to_fractional(const Site& site) const;
};

class WyckoffType {
public:
  int multiplicity;
  char letter;
  int some_variability;
  int site_rotations;
  int site_mirrors;
  std::vector<ImageType> image;
};

class Site {
public:
  WyckoffType* wyckoff;
  Basis* x;
  Basis* y;
  Basis* angle;
  bool flip_site = false;

  Vect3 site_variables() const {
    return Vect3(this->x->get_value(), this->y->get_value(), this->angle->get_value());
  };
  Vect2 get_position() const {
    return Vect2(this->x->get_value(), this->y->get_value());
  };
  int get_flip_sign() const { return this->flip_site ? 1 : -1; };
  int get_multiplicity() const { return this->wyckoff->multiplicity; };
};

class ShapeInstance {
public:
  ShapeInstance(const Shape& shape, Site& site, ImageType& image)
      : shape(&shape), site(&site), image(&image){};
  const Shape* shape;
  Site* site;
  ImageType* image;

  bool operator==(const ShapeInstance& other) const {
    return (
        this->shape == other.shape && this->site == other.site &&
        this->image == other.image);
  };

  Vect2 get_fractional_coordinates() const {
    return this->image->real_to_fractional(*this->site);
  }
  Vect2 get_real_coordinates() const { return this->site->get_position(); };
  double get_angle() const { return this->site->angle->get_value(); };
  double get_rotational_offset() const { return this->image->rotation_offset; };
  bool get_flipped() const { return this->image->flipped ^ this->site->flip_site; };
  bool pair_clash(ShapeInstance& other) const;
};

class WallpaperGroup {
public:
  WallpaperGroup(std::string label, std::vector<WyckoffType> wyckoffs);
  std::string label;
  bool a_b_equal = false;
  bool hexagonal = false;
  bool rectangular = false;
  int num_symmetries = 0;
  std::vector<WyckoffType> wyckoffs;
  int num_wyckoffs;

  size_t group_multiplicity() const;
};

std::string create_filename(
    const Shape& shape,
    const WallpaperGroup& group,
    const std::string site_list,
    const std::string& directory);

#endif /* SHAPES_H */
