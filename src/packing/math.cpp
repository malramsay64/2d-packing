/*
 * math.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "math.h"

#include <cmath>
#include <cstddef>

#include <pybind11/operators.h>

namespace py = pybind11;

double positive_modulo(const double i, const double n) {
  return std::fmod(std::fmod(i, n) + n, n);
}

int positive_modulo(const int i, const int n) { return ((i % n) + n) % n; }

int sign(double val) { return (double{0} < val) - (val < double{0}); }

Vect2 Vect2::operator+(const Vect2& other) const {
  return Vect2(this->x + other.x, this->y + other.y);
}
Vect2 Vect2::operator-(const Vect2& other) const {
  return Vect2(this->x - other.x, this->y - other.y);
}
Vect2 Vect2::operator*(const Vect2& other) const {
  return Vect2(this->x * other.x, this->y * other.y);
}
Vect2 Vect2::operator*(const float other) const {
  return Vect2(this->x * other, this->y * other);
}
Vect2 Vect2::operator==(const Vect2& other) const {
  return Vect2(this->x == other.x, this->y == other.y);
}
Vect2 Vect2::operator%(const Vect2& other) const {
  return Vect2(std::fmod(this->x, other.x), std::fmod(this->y, other.y));
}
Vect2 Vect2::operator%(const double other) const {
  return Vect2(std::fmod(this->x, other), std::fmod(this->y, other));
}

double Vect2::norm_sq() const { return this->x * this->x + this->y * this->y; }
double Vect2::norm() const { return std::sqrt(this->norm_sq()); }

Vect2& positive_modulo(Vect2& v, const double modulo) {
  v.x = positive_modulo(v.x, modulo);
  v.y = positive_modulo(v.y, modulo);
  return v;
}

double temperature_distribution(
    const double old_val,
    const double new_val,
    const double kT,
    const std::size_t replicas) {
  return std::exp(
      ((1 / old_val - 1 / new_val) / kT) + replicas * std::log(old_val / new_val));
}

bool is_close(const float value, const float expected, const float rel_tol) {
  return fabs(value - expected) < rel_tol * expected;
}

void export_Vect2(py::module& m) {
  py::class_<Vect2> vect2(m, "Vect2");
  vect2.def(py::init<double, double>(), py::arg("x") = 0, py::arg("y") = 0)
      .def_readwrite("x", &Vect2::x)
      .def_readwrite("y", &Vect2::y)
      .def("norm", &Vect2::norm)
      .def("norm_sq", &Vect2::norm_sq)
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self * py::self)
      .def(py::self * float());
}
