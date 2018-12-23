/*
 * geometry.h
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "math.h"

#ifndef GEOMETRY_H
#define GEOMETRY_H

// To find orientation of ordered triplet (a, b, c).
// The function returns following values
// 0 --> a, b and c are colinear
// 1 --> Clockwise
// -1 --> Counterclockwise
int triplet_orientation(const Vect2& a, const Vect2& b, const Vect2& c);

// Given three colinear points a, b, c, the function checks if
// point b lies on line segment 'ac'
bool on_segment(const Vect2& a, const Vect2& b, const Vect2& c);

bool segments_cross(const Vect2& A1, const Vect2& A2, const Vect2& B1, const Vect2& B2);

#endif /* !GEOMETRY_H */