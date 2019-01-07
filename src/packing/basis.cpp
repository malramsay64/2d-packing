/*
 * basis.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "basis.h"

#include <algorithm>
#include <cmath>
#include <ostream>
#include <sstream>
#include <string>

#include "shapes.h"
#include "wallpaper.h"

namespace py = pybind11;

double Basis::value_range() const {
  return this->max_val - this->min_val;
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
double Basis::get_value() const {
  return this->value;
}

void Basis::reset_value() {
  this->value = this->value_previous;
}

double Basis::get_random_value(const double kT) const {
  return this->value + this->step_size * this->value_range() * (fluke() - 0.5);
}

Vect3 OccupiedSite::site_variables() const {
  return Vect3(this->x->get_value(), this->y->get_value(), this->angle->get_value());
}

Vect2 OccupiedSite::get_position() const {
  return Vect2(this->x->get_value(), this->y->get_value());
}

int OccupiedSite::get_multiplicity() const {
  return this->wyckoff->multiplicity();
}

std::string OccupiedSite::str() const {
  std::stringstream ss;
  ss << *this;
  return ss.str();
}

std::ostream& operator<<(std::ostream& os, const OccupiedSite& site) {
  os << "Site: " << site.wyckoff->letter << std::endl;
  for (const auto& symmetry : site.wyckoff->symmetries) {
    os << " - " << symmetry.real_to_fractional(site.get_position())
       << "angle: " << site.angle->get_value() + symmetry.rotation_offset << std::endl;
  }
  return os;
}

double Cell::area() const {
  return this->x_len->get_value() * this->y_len->get_value() *
         std::fabs(std::sin(this->angle->get_value()));
}

Vect2 Cell::fractional_to_real(const Vect2& fractional) const {
  Vect2 v(0, 0);
  v.x = fractional.x * this->x_len->get_value() +
        fractional.y * this->y_len->get_value() * cos(this->angle->get_value());
  v.y = fractional.y * this->y_len->get_value() * sin(this->angle->get_value());
  return v;
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
}

void CellAngleBasis::reset_value() {
  this->value = this->value_previous;
  this->reset_cell_lengths();
}

void CellAngleBasis::update_cell_lengths() {
  this->cell_x_len->set_value(
      this->cell_x_len->get_value() *
      sqrt(sin(this->value_previous) / sin(this->value)));
  this->cell_y_len->set_value(
      this->cell_y_len->get_value() *
      sqrt(sin(this->value_previous) / sin(this->value)));
}

void CellAngleBasis::reset_cell_lengths() {
  this->cell_x_len->reset_value();
  this->cell_y_len->reset_value();
}

double MirrorBasis::get_random_value(const double kT) const {
  if ((this->mirrors % 2 == 0) && (fluke() < 0.5)) {
    /* turn it 90 degrees to switch the x and y mirror planes */
    if (this->value < 3 * M_PI_4) {
      return this->value + PI / this->mirrors;
    }
    return this->value - PI / this->mirrors;
  }
  /* turn it 180 degrees so that all mirror planes are preserved */
  return positive_modulo(this->value + M_PI, M_2_PI);
}

void export_Basis(py::module& m) {
  py::class_<Basis, std::shared_ptr<Basis>> basis(m, "Basis");
  basis
      .def(
          py::init<const double, const double, const double, const double>(),
          py::arg("value"),
          py::arg("min_val"),
          py::arg("max_val"),
          py::arg("step_size") = 0.01)
      .def_property("value", &Basis::get_value, &Basis::set_value)
      .def("value_range", &Basis::value_range)
      .def("reset_value", &Basis::reset_value)
      .def("get_random_value", &Basis::get_random_value);

  py::class_<CellLengthBasis, Basis, std::shared_ptr<CellLengthBasis>>
      cell_length_basis(m, "CellLengthBasis");
  cell_length_basis
      .def(
          py::init<const double, const double, const double, const double>(),
          py::arg("value"),
          py::arg("min_val"),
          py::arg("max_val"),
          py::arg("step_size") = 0.01)
      .def_property("value", &CellLengthBasis::get_value, &CellLengthBasis::set_value)
      .def("get_random_value", &CellLengthBasis::get_random_value);

  py::class_<CellAngleBasis, Basis, std::shared_ptr<CellAngleBasis>> cell_angle_basis(
      m, "CellAngleBasis");
  cell_angle_basis
      .def(
          py::init<
              const double,
              const double,
              const double,
              const double,
              std::shared_ptr<Basis>,
              std::shared_ptr<Basis>>(),
          py::arg("value"),
          py::arg("min_val"),
          py::arg("max_val"),
          py::arg("step_size"),
          py::arg("cell_x_len"),
          py::arg("cell_y_len"))
      .def_property("value", &CellAngleBasis::get_value, &CellAngleBasis::set_value)
      .def_readonly("cell_x_len", &CellAngleBasis::cell_x_len)
      .def_readonly("cell_y_len", &CellAngleBasis::cell_y_len)
      .def("reset_value", &CellAngleBasis::reset_value);

  py::class_<FixedBasis, Basis, std::shared_ptr<FixedBasis>> fixed_basis(
      m, "FixedBasis");
  fixed_basis.def(py::init<const double>(), py::arg("value"));

  py::class_<MirrorBasis, Basis, std::shared_ptr<MirrorBasis>> mirror_basis(
      m, "MirrorBasis");
  mirror_basis.def(
      py::init<const double, const double, const double, const int>(),
      py::arg("value"),
      py::arg("min_val"),
      py::arg("max_val"),
      py::arg("mirrors"));
}
