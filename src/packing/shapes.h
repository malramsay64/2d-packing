#include "fluke.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

#ifndef SHAPES_H
#define SHAPES_H

const double PI = std::atan(1.0) * 4;

double positive_modulo(double i, double n) { return std::fmod(std::fmod(i, n) + n, n); }

int positive_modulo(int i, int n) { return ((i % n) + n) % n; }

struct Vect2 {
  double x;
  double y;
  Vect2(double x, double y) : x(x), y(y){};
  Vect2() : x(0), y(0){};

  Vect2 operator+(const Vect2& other) const {
    return Vect2(this->x + other.x, this->y + other.y);
  }

  Vect2 operator-(const Vect2& other) const {
    return Vect2(this->x - other.x, this->y - other.y);
  }

  Vect2 operator*(const Vect2& other) const {
    return Vect2(this->x * other.x, this->y * other.y);
  }

  inline double norm_sq() { return this->x * this->x + this->y * this->y; }
  inline double norm() { return sqrt(this->norm_sq()); }
};

struct Vect3 {
  double x;
  double y;
  double z;

  Vect3(double x, double y, double z) : x(x), y(y), z(z){};
};

struct Int2 {
  int x;
  int y;

  bool operator==(const Int2& other) const {
    return this->x == other.x && this->y == other.y;
  }
};

struct Basis {
  double* value;
  double min_val;
  double max_val;
  bool fixed = false;
  bool mirror = false;
  bool cell_side = false;
  bool cell_angle = false;

  Basis(
      double* value,
      double min_val,
      double max_val,
      bool fixed,
      bool mirror,
      bool cell_side)
      : value(value), min_val(min_val), max_val(max_val), fixed(fixed), mirror(mirror),
        cell_side(cell_side){};
  Basis(double* value, bool mirror) : Basis(value, *value, *value, true, true, false){};
  Basis(double* value) : Basis(value, *value, *value, true, false, false){};
  Basis() : Basis(0){};

  double value_range() const { return this->max_val - this->min_val; };
  void validate() {
    if (*this->value < min_val) {
      *this->value = min_val;
    } else if (*this->value > max_val) {
      *this->value = max_val;
    }
  }
};

struct Cell {
  double x_len;
  double y_len;
  double angle;

  Cell() : x_len(0), y_len(0), angle(0){};

  Vect2 fractional_to_real(const Vect2&);
  double area() const { return this->x_len * this->y_len * fabs(sin(this->angle)); };
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

class ShapeInstance {
public:
  ShapeInstance(Shape* shape);
  Shape* shape;
  double x, y;
  double theta;
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
  ImageType(Vect3 x_coeffs, Vect3 y_coeffs, int rotation_offset);
  Vect3 x_coeffs;
  Vect3 y_coeffs;
  /* first index is the output coordinate,
   second index is input coordinate (+1 is for the constant that is added)

   for example, to work out the new values:
   x_new = coord_coeffs[0][0]*x_old + coord_coeffs[0][1]*y_old +
   coord_coeffs[0][2]; y_new = coord_coeffs[1][0]*x_old +
   coord_coeffs[1][1]*y_old + coord_coeffs[1][2];
  */
  int rotation_offset; /*angles move clockwise*/
  bool flipped;        /* =1, when flipped using the x-axis as mirror, rotation_offset
                          is then employed if mirror is at y-axis for example */
  mirror site_mirror;
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

struct Site {
  size_t wyckoff_index;
  double x = 0;
  double y = 0;
  double angle = 0;
  bool flip_site = false;

  Vect3 site_variables() const { return Vect3(this->x, this->y, this->angle); };
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

  const WyckoffType& get_wyckoff(const Site& site) const {
    return this->wyckoffs[site.wyckoff_index];
  }

  WyckoffType& get_wyckoff(const Site& site) {
    return this->wyckoffs[site.wyckoff_index];
  }

  size_t group_multiplicity() const;
};

std::string create_filename(
    const Shape& shape,
    const WallpaperGroup& group,
    const std::string site_list,
    const std::string& directory);

#endif /* SHAPES_H */
