/*
 * module.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <pybind11/pybind11.h>

double fluked() { return 0.; }

PYBIND11_MODULE(_packing, m) {
  m.attr("__name__") = "pypacking._packing";
  m.doc() = R"pbdoc(
      Pybind11 test
      -------------
      )pbdoc";
  m.def("fluke", &fluked, R"pbdoc(
      Generate a random number
      )pbdoc");

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}
