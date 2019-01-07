/*
 * module.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <pybind11/pybind11.h>

#include "basis.h"
#include "geometry.h"
#include "math.h"
#include "random.h"
#include "shapes.h"
#include "util.h"
#include "wallpaper.h"

PYBIND11_MODULE(_packing, m) {
  m.attr("__name__") = "pypacking._packing";
  m.doc() = R"pbdoc(
      Pybind11 test
      -------------
      )pbdoc";
  export_fluke(m);
  export_Shape(m);
  export_Vect2(m);
  export_Vect3(m);
  export_geometry(m);
  export_Basis(m);
  export_combinations(m);
  export_Mirror(m);
  export_SymmetryTransform(m);
  export_WyckoffSite(m);
  export_WallpaperGroup(m);

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}
