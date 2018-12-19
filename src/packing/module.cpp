/*
 * module.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "random.h"
#include <pybind11/pybind11.h>

PYBIND11_MODULE(packing, m) {
  m.doc() = R"pbdoc(
      Pybind11 test
      -------------
      )pbdoc";
  m.def("fluke", &fluke, R"pbdoc(
      Generate a random number
      )pbdoc");

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}
