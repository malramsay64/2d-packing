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

/** Find the orientation of an ordered triplet of points; a, b, and c.
 *
 * This is adapted from a post which has additional details on the algorithm.
 * https://www.geeksforgeeks.org/orientation-3-ordered-points/
 *
 * \param a First point in ordered triplet
 * \param b Second point in ordered triplet
 * \param c Third point in ordered triplet
 *
 * \return Integer describing the orietation of the triplet.
 *          - 0 --> a, b, and c are colinear (on the same line)
 *          - 1 --> A clockwise orientation
 *          - -1 --> Anti-Clockwise orientation
 */
int triplet_orientation(const Vect2& a, const Vect2& b, const Vect2& c) {
  return sign((b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y));
}

/** Checks whether the point b lies on the line segment ac.
 *
 * This function assumes the points are on the same line, so should be used after
 * triplet_orientation has returned the value 0 indicating co-linearity.
 *
 * \param a The first point of the line segment
 * \param b The point we are checking is on the line segment
 * \paarm c The second point of the line segment
 */
bool on_segment(const Vect2& a, const Vect2& b, const Vect2& c) {
  if (b.x <= std::max(a.x, c.x) && b.x >= std::min(a.x, c.x) &&
      b.y <= std::max(a.y, c.y) && b.y >= std::min(a.y, c.y)) {
    return true;
  }

  return false;
}

/** Evaluate whether the line segment A1B1 crosses the line segment A2B2
 *
 * \param A1 The first point of the first line segment
 * \param B1 The second point of the first line segment
 * \param A2 The first point of the second line segment
 * \param B2 The second point of the second line segment
 *
 * \returns bool indicating whether the line segments cross at some point.
 */
bool segments_cross(
    const Vect2& A1,
    const Vect2& B1,
    const Vect2& A2,
    const Vect2& B2) {

  // Find the four orientations needed for general and
  // special cases
  const int o1{triplet_orientation(A1, B1, A2)};
  const int o2{triplet_orientation(A1, B1, B2)};
  const int o3{triplet_orientation(A2, B2, A1)};
  const int o4{triplet_orientation(A2, B2, B1)};

  // General case
  if (o1 != o2 && o3 != o4)
    return true;

  // Special Cases
  // A1, B1 and A2 are colinear and A2 lies on segment A1B1
  if (o1 == 0 && on_segment(A1, A2, B1))
    return true;

  // A1, B1 and B2 are colinear and B2 lies on segment A1B1
  if (o2 == 0 && on_segment(A1, B2, B1))
    return true;

  // A2, B2 and A1 are colinear and A1 lies on segment A2B2
  if (o3 == 0 && on_segment(A2, A1, B2))
    return true;

  // A2, B2 and A1 are colinear and B1 lies on segment A2B2
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
      py::arg("B1"),
      py::arg("A2"),
      py::arg("B2"));
}
