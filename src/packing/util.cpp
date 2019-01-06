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
