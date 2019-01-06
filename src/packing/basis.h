/*
 * basis.h
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <memory>
#include <vector>

#include <pybind11/pybind11.h>

#include "math.h"
#include "random.h"

#ifndef BASIS_H
#define BASIS_H

// Forward declaration of WyckoffSite
class WyckoffSite;

class Basis {
protected:
  double value_previous;
  double value;

public:
  const double min_val;
  const double max_val;
  const double step_size = 0.01;

  Basis(
      const double value,
      const double min_val,
      const double max_val,
      const double step_size)
      : value(value), min_val(min_val), max_val(max_val), step_size(step_size){};
  Basis(const double value, const double min_val, const double max_val)
      : Basis(value, min_val, max_val, 0.01){};

  double value_range() const;
  double get_value() const;
  void set_value(double new_value);
  void reset_value();
  double get_random_value(const double kT) const;
};

class OccupiedSite {
public:
  std::shared_ptr<WyckoffSite> wyckoff;
  std::shared_ptr<Basis> x;
  std::shared_ptr<Basis> y;
  std::shared_ptr<Basis> angle;

  Vect3 site_variables() const;
  Vect2 get_position() const;
  int get_multiplicity() const;

  std::string str() const;
  friend std::ostream& operator<<(std::ostream&, const OccupiedSite&);
};

struct Cell {
  std::shared_ptr<Basis> x_len;
  std::shared_ptr<Basis> y_len;
  std::shared_ptr<Basis> angle;

  Vect2 fractional_to_real(const Vect2&) const;
  double area() const;
};

class CellLengthBasis : public Basis {
public:
  const double step_size;

  CellLengthBasis(
      const double value,
      const double min_val,
      const double max_val,
      const double step_size)
      : Basis(value, min_val, max_val), step_size(step_size){};

  double get_random_value(const double kT) const;
};

class CellAngleBasis : public Basis {
private:
  void update_cell_lengths();
  void reset_cell_lengths();

public:
  const double step_size;
  std::shared_ptr<Basis> cell_x_len;
  std::shared_ptr<Basis> cell_y_len;

  CellAngleBasis(
      const double value,
      const double min_val,
      const double max_val,
      const double step_size,
      std::shared_ptr<Basis> cell_x_len,
      std::shared_ptr<Basis> cell_y_len)
      : Basis(value, min_val, max_val), step_size(step_size), cell_x_len(cell_x_len),
        cell_y_len(cell_y_len){};

  void set_value(double new_value);
  void reset_value();
  double get_random_value(const double kT) const;
};

class FixedBasis : public Basis {
public:
  FixedBasis(const double value) : Basis(value, value, value){};

  void set_value(double new_value){};
  void reset_value(){};
};

class MirrorBasis : public Basis {
public:
  const int mirrors;

  MirrorBasis(
      const double value,
      const double min_val,
      const double max_val,
      const int mirrors)
      : Basis(value, min_val, max_val), mirrors(mirrors){};

  double get_random_value(const double kT) const;
};

void export_Basis(pybind11::module& m);

#endif /* !BASIS_H */
