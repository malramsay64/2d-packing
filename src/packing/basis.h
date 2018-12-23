/*
 * basis.h
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <vector>

#include <pybind11/pybind11.h>

#include "math.h"
#include "random.h"
#include "shapes.h"

#ifndef BASIS_H
#define BASIS_H

class Basis {
protected:
  double value_previous;
  double value;

public:
  const double min_val;
  const double max_val;
  const double step_size = 1;

  Basis(
      const double value,
      const double min_val,
      const double max_val,
      const double step_size)
      : value(value), min_val(min_val), max_val(max_val), step_size(step_size){};
  Basis(const double value, const double min_val, const double max_val)
      : Basis(value, min_val, max_val, 1){};

  double value_range() const;
  double get_value() const;
  void set_value(double new_value);
  void reset_value();
  double get_random_value(const double kT) const;
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
  Basis* const cell_x_len;
  Basis* const cell_y_len;

  CellAngleBasis(
      const double value,
      const double min_val,
      const double max_val,
      const double step_size,
      Basis* const cell_x_len,
      Basis* const cell_y_len)
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

class FlipBasis : public Basis {
private:
  std::vector<Site>* const occupied_sites;
  int value_previous;

public:
  FlipBasis(std::vector<Site>* const occupied_sites)
      : Basis(0, 0, occupied_sites->size()), occupied_sites(occupied_sites){};

  double get_random_value(const double kT) const;
  void set_value(double new_value);
  void reset_value();
};

void export_Basis(pybind11::module& m);

#endif /* !BASIS_H */
