/*
 * math.h
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef MATH_H
#define MATH_H

#include <cmath>

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
  Vect2 operator*(const float other) const {
    return Vect2(this->x * other, this->y * other);
  }
  Vect2 operator==(const Vect2& other) const {
    return Vect2(this->x == other.x, this->y == other.y);
  }

  inline double norm_sq() { return this->x * this->x + this->y * this->y; }
  inline double norm() { return sqrt(this->norm_sq()); }
};

Vect2& positive_modulo(Vect2& v, double modulo) {
  v.x = positive_modulo(v.x, modulo);
  v.y = positive_modulo(v.y, modulo);
  return v;
}

struct Vect3 {
  double x;
  double y;
  double z;

  Vect3(double x, double y, double z) : x(x), y(y), z(z){};
};

#endif /* !MATH_H */
