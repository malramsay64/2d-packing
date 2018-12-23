/*
 * basis.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "basis.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace py = pybind11;

double Basis::value_range() const { return this->max_val - this->min_val; };

void Basis::set_value(double new_value) {
  if (new_value < this->min_val) {
    new_value = this->min_val;
  } else if (new_value > this->max_val) {
    new_value = this->max_val;
  }
  this->value_previous = this->value;
  this->value = new_value;
}
double Basis::get_value() const { return this->value; };

void Basis::reset_value() { this->value = this->value_previous; };

double Basis::get_random_value(const double kT) const {
  return this->value + this->step_size * this->value_range() * (fluke() - 0.5);
}

double CellLengthBasis::get_random_value(const double kT) const {
  return this->get_value() + (1.0 + std::min(3.0 * kT, 0.1) * (fluke() - 0.5));
}

double CellAngleBasis::get_random_value(const double kT) const {
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

double MirrorBasis::get_random_value(const double kT) const {
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

double FlipBasis::get_random_value(const double kT) const {
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

void export_Basis(py::module& m) {}
