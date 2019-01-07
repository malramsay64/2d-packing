/*
 * math.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "math.h"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <ostream>
#include <sstream>

#include <pybind11/operators.h>

namespace py = pybind11;

double positive_modulo(const double i, const double n) {
  return std::fmod(std::fmod(i, n) + n, n);
}

int positive_modulo(const int i, const int n) {
  return ((i % n) + n) % n;
}

int sign(double val) {
  return (double{0} < val) - (val < double{0});
}

bool Vect2::operator==(const Vect2& other) const {
  return this->x == other.x && this->y == other.y;
}

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

Vect2 Vect2::operator%(const Vect2& other) const {
  return Vect2(std::fmod(this->x, other.x), std::fmod(this->y, other.y));
}

Vect2 Vect2::operator%(const double other) const {
  return Vect2(std::fmod(this->x, other), std::fmod(this->y, other));
}

bool Vect3::operator==(const Vect3& other) const {
  return this->x == other.x && this->y == other.y && this->z == other.z;
}

std::ostream& operator<<(std::ostream& os, const Vect2& vect2) {
  os << "x: " << vect2.x << " y: " << vect2.y;
  return os;
}

std::string Vect2::str() const {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(5);
  ss << "Vect2(x=" << this->x << ", y=" << this->y << ")";
  return ss.str();
}

double Vect2::norm_sq() const {
  return this->x * this->x + this->y * this->y;
}

double Vect2::norm() const {
  return std::sqrt(this->norm_sq());
}

Vect2& positive_modulo(Vect2& v, const double modulo) {
  v.x = positive_modulo(v.x, modulo);
  v.y = positive_modulo(v.y, modulo);
  return v;
}

std::string Vect3::str() const {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(5);
  ss << "Vect3(x=" << this->x << ", y=" << this->y << ", z=" << this->z << ")";
  return ss.str();
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
      .def("__repr__", [](const Vect2& v) { return "<" + v.str() + ">"; })
      .def(py::self == py::self)
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self * py::self)
      .def(py::self * float());
}

void export_Vect3(py::module& m) {
  py::class_<Vect3> vect3(m, "Vect3");
  vect3
      .def(
          py::init<double, double, double>(),
          py::arg("x") = 0,
          py::arg("y") = 0,
          py::arg("z") = 0)
      .def_readwrite("x", &Vect3::x)
      .def_readwrite("y", &Vect3::y)
      .def_readwrite("z", &Vect3::z)
      .def("__repr__", [](const Vect3& v) { return "<" + v.str() + ">"; })
      .def(py::self == py::self);
}
