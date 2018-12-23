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

double positive_modulo(double i, double n) { return std::fmod(std::fmod(i, n) + n, n); }

int positive_modulo(int i, int n) { return ((i % n) + n) % n; }

template <typename T> int sign(const T val) { return (T(0) < val) - (val < T(0)); }

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
double Vect2::norm() const { return sqrt(this->norm_sq()); }

Vect2& positive_modulo(Vect2& v, double modulo) {
  v.x = positive_modulo(v.x, modulo);
  v.y = positive_modulo(v.y, modulo);
  return v;
}

double temperature_distribution(
    double old_val,
    double new_val,
    double kT,
    std::size_t replicas) {
  return std::exp(
      ((1.0 / old_val - 1.0 / new_val) / kT) + replicas * std::log(old_val / new_val));
}

bool is_close(float value, float expected, float rel_tol) {
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
      .def(py::self * float())
      .def(py::self == py::self);
}
