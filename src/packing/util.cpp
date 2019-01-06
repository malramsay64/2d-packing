/*
 * util.cpp
 * Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "util.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "shapes.h"

namespace py = pybind11;

std::size_t calculate_shape_replicas(const std::vector<OccupiedSite>& occupied_sites) {
  std::size_t replicas{0};
  for (const OccupiedSite& site : occupied_sites) {
    replicas += site.get_multiplicity();
  }
  return replicas;
}

char compute_chiral_state(const std::vector<OccupiedSite>& occupied_sites) {
  int chiralsum{0};
  std::size_t totalsum{0};
  for (const OccupiedSite& site : occupied_sites) {
    chiralsum += site.get_flip_sign() * site.get_multiplicity();
    totalsum += site.get_multiplicity();
  }
  if (chiralsum != 0) {
    return (chiralsum == totalsum) ? 'c' : 's';
  }
  return 'a';
}

void export_combinations(py::module& m) {
  m.def("combinations", &combinations<int>, py::arg("values"), py::arg("take"));
  m.def("combinations", &combinations<double>);
  m.def("uniqueify", [](std::vector<int> v) {
    uniqueify<int>(v);
    return v;
  });
  m.def("uniqueify", [](std::vector<double> v) {
    uniqueify<double>(v);
    return v;
  });
  m.def("uniqueify", [](std::vector<std::vector<int>> v) {
    uniqueify<std::vector<int>>(v);
    return v;
  });
  m.def("uniqueify", [](std::vector<std::vector<double>> v) {
    uniqueify<std::vector<double>>(v);
    return v;
  });
}
