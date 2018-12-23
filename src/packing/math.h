/*
 * math.h
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <pybind11/pybind11.h>

#ifndef MATH_H
#define MATH_H

const double PI = M_PI;

double positive_modulo(double i, double n);

int positive_modulo(int i, int n);

bool is_close(float value, float expected, float rel_tol = 1e-8);

template <typename T> int sign(const T val);

struct Vect2 {
  double x;
  double y;
  Vect2(double x, double y) : x(x), y(y){};
  Vect2() : x(0), y(0){};

  Vect2 operator+(const Vect2& other) const;
  Vect2 operator-(const Vect2& other) const;
  Vect2 operator*(const Vect2& other) const;
  Vect2 operator*(const float other) const;
  Vect2 operator==(const Vect2& other) const;
  Vect2 operator%(const Vect2& other) const;
  Vect2 operator%(const double other) const;

  double norm_sq() const;
  double norm() const;
};

Vect2& positive_modulo(Vect2& v, double modulo);

struct Vect3 {
  double x;
  double y;
  double z;

  Vect3(double x, double y, double z) : x(x), y(y), z(z){};
};

double temperature_distribution(
    double old_val,
    double new_val,
    double kT,
    std::size_t replicas);

void export_Vect2(pybind11::module& m);

#endif /* !MATH_H */
