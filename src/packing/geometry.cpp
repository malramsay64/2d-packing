/*
 * geometry.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "geometry.h"

#include <algorithm>
#include <cmath>

namespace py = pybind11;

// To find orientation of ordered triplet (a, b, c).
// The function returns following values
// 0 --> a, b and c are colinear
// 1 --> Clockwise
// -1 --> Counterclockwise
int triplet_orientation(const Vect2& a, const Vect2& b, const Vect2& c) {
  // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
  // for details of below formula.
  Vect2 tmp = (b - a) * (c - b);
  return sign(tmp.x - tmp.y);
}

// Given three colinear points a, b, c, the function checks if
// point b lies on line segment 'ac'
bool on_segment(const Vect2& a, const Vect2& b, const Vect2& c) {
  if (b.x <= std::max(a.x, c.x) && b.x >= std::min(a.x, c.x) &&
      b.y <= std::max(a.y, c.y) && b.y >= std::min(a.y, c.y))
    return true;

  return false;
}

bool segments_cross(
    const Vect2& A1,
    const Vect2& A2,
    const Vect2& B1,
    const Vect2& B2) {

  // Find the four orientations needed for general and
  // special cases
  int o1 = triplet_orientation(A1, B1, A2);
  int o2 = triplet_orientation(A1, B1, B2);
  int o3 = triplet_orientation(A2, B2, A1);
  int o4 = triplet_orientation(A2, B2, B1);

  // General case
  if (o1 != o2 && o3 != o4)
    return true;

  // Special Cases
  // A1, B1 and A2 are colinear and A2 lies on segment p1q1
  if (o1 == 0 && on_segment(A1, A2, B1))
    return true;

  // p1, q1 and q2 are colinear and q2 lies on segment p1q1
  if (o2 == 0 && on_segment(A1, B2, B1))
    return true;

  // p2, q2 and p1 are colinear and p1 lies on segment p2q2
  if (o3 == 0 && on_segment(A2, A1, B2))
    return true;

  // p2, q2 and q1 are colinear and q1 lies on segment p2q2
  if (o4 == 0 && on_segment(A2, B1, B2))
    return true;

  return false; // Doesn't fall in any of the above cases
}

void export_geometry(py::module& m) {
  m.def(
      "triplet_orientation",
      &triplet_orientation,
      py::arg("a"),
      py::arg("b"),
      py::arg("c"));
  m.def("on_segment", &on_segment, py::arg("a"), py::arg("b"), py::arg("c"));
  m.def(
      "segments_cross",
      &segments_cross,
      py::arg("A1"),
      py::arg("A2"),
      py::arg("B1"),
      py::arg("B2"));
}
