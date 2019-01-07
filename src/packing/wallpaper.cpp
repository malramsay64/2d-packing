/*
 * wallpaper.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "wallpaper.h"

#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "basis.h"
#include "geometry.h"
#include "math.h"
#include "shapes.h"
#include "util.h"

namespace py = pybind11;

bool SymmetryTransform::operator==(const SymmetryTransform& other) const {
  return (
      this->x_coeffs == other.x_coeffs && this->y_coeffs == other.y_coeffs &&
      this->rotation_offset == other.rotation_offset &&
      this->site_mirror == other.site_mirror);
}

std::ostream& operator<<(std::ostream& os, const SymmetryTransform& symmetry) {
  os << " - " << symmetry.real_to_fractional(Vect2{0, 0})
     << "angle: " << symmetry.rotation_offset << std::endl;
  return os;
}

Vect2 SymmetryTransform::real_to_fractional(const Vect2& real) const {
  /* converts site variables and wyckoff site coefficients into the location
   * of the actual wyckoff image in fractional coordinates */
  Vect2 v(
      this->x_coeffs.x * real.x + this->x_coeffs.y * real.y + this->x_coeffs.z,
      this->y_coeffs.x * real.x + this->y_coeffs.y * real.y + this->y_coeffs.z);
  return positive_modulo(v, 1.);
}

std::string SymmetryTransform::str() const {
  std::stringstream ss;
  ss << *this;
  return ss.str();
}

bool WyckoffSite::operator==(const WyckoffSite& other) const {
  return (
      this->letter == other.letter && this->variability == other.variability &&
      this->rotations == other.rotations && this->mirrors == other.mirrors &&
      this->symmetries == other.symmetries);
}

std::size_t WyckoffSite::multiplicity() const {
  return this->symmetries.size();
}

bool WyckoffSite::vary_x() const {
  // If the first symmetry can vary x, all of them can
  return fabs(this->symmetries[0].x_coeffs.x) > 0.1;
}

bool WyckoffSite::vary_y() const {
  // If the first symmetry can vary x, all of them can
  return fabs(this->symmetries[0].y_coeffs.y) > 0.1;
}

int WyckoffSite::mirror_type() const {
  return static_cast<int>(this->symmetries.front().site_mirror);
}

std::string WyckoffSite::str() const {
  std::stringstream ss;
  ss << "WyckoffSite(letter= " << this->letter << ", variability=" << this->variability
     << ", rotations=" << this->rotations << ", mirrors=" << this->mirrors << ")"
     << std::endl;
  for (const auto& symmetry : this->symmetries) {
    ss << symmetry;
  }
  return ss.str();
}

std::string IsopointalGroup::group_string() const {
  std::stringstream ss;
  for (const auto& site : this->wyckoff_sites) {
    ss << site.letter;
  }
  return ss.str();
}

std::size_t IsopointalGroup::group_multiplicity() const {
  std::size_t multiplicity{0};
  for (const WyckoffSite& site : this->wyckoff_sites) {
    multiplicity += site.multiplicity();
  }
  return multiplicity;
}

std::vector<IsopointalGroup> generate_isopointal_groups(
    const Shape& shape,
    const WallpaperGroup& group,
    std::size_t num_occupied_sites) {
  auto console = spdlog::stdout_color_mt("console");

  std::vector<WyckoffSite> valid_sites;
  for (const WyckoffSite& wyckoff : group.wyckoff_sites) {
    // first test if the shape has the required symmetries for various sites
    // at the moment only tests rotations & mirrors... that's all?
    int rotational_match = shape.rotational_symmetries % wyckoff.rotations;
    int mirror_match = shape.mirrors % wyckoff.mirrors;
    if (rotational_match == 0) {
      if ((wyckoff.mirrors == 0) || ((mirror_match == 0) && (shape.mirrors != 0))) {
        if (wyckoff.variability) {
          console->debug("site %c is valid and variable\n", wyckoff.letter);
        } else {
          printf("site %c is valid\n", wyckoff.letter);
        }
        if (wyckoff.variability) {
          // Variable sites are added with replacement, so add one for each occupied
          // site
          valid_sites.insert(valid_sites.end(), num_occupied_sites, wyckoff);
        } else {
          valid_sites.push_back(wyckoff);
        }
      } else {
        console->debug(
            "site %c does not have the required rotational symmetry\n", wyckoff.letter);
      }
    }
  }

  auto occupied_sites = combinations<WyckoffSite>(valid_sites, num_occupied_sites);
  uniqueify<std::vector<WyckoffSite>>(occupied_sites);

  console->info(
      "enumeration complete: in fact there were %lu ways\n", occupied_sites.size());

  std::vector<IsopointalGroup> isopointal_groups;
  for (auto& combination : occupied_sites) {
    isopointal_groups.push_back(IsopointalGroup(combination));
    console->debug("%s\n", isopointal_groups.back().group_string());
  }
  return isopointal_groups;
}

void export_Mirror(py::module& m) {
  py::enum_<Mirror>(m, "Mirror", py::arithmetic())
      .value("m0", Mirror::m0)
      .value("m30", Mirror::m30)
      .value("m45", Mirror::m45)
      .value("m60", Mirror::m60)
      .value("m90", Mirror::m90)
      .value("m135", Mirror::m135)
      .value("m300", Mirror::m300)
      .value("m330", Mirror::m330)
      .export_values();
}

void export_SymmetryTransform(py::module& m) {
  py::class_<SymmetryTransform> symmetry_transform(m, "SymmetryTransform");
  symmetry_transform
      .def(
          py::init<const Vect3&, const Vect3&, const double, Mirror>(),
          py::arg("x_coeffs"),
          py::arg("y_coeffs"),
          py::arg("rotation_offset") = 0,
          py::arg("mirror") = Mirror::m0)
      .def("real_to_fractional", &SymmetryTransform::real_to_fractional);
}

void export_WyckoffSite(py::module& m) {
  py::class_<WyckoffSite> wyckoff_site(m, "WyckoffSite");
  wyckoff_site
      .def(
          py::init<
              const char,
              const std::vector<SymmetryTransform>&,
              const std::size_t,
              const std::size_t,
              const std::size_t>(),
          py::arg("letter"),
          py::arg("symmetries"),
          py::arg("variability") = 0,
          py::arg("rotations") = 1,
          py::arg("mirrors") = 0)
      .def("multiplicity", &WyckoffSite::multiplicity)
      .def("vary_x", &WyckoffSite::vary_x)
      .def("vary_y", &WyckoffSite::vary_y)
      .def("mirror_type", &WyckoffSite::mirror_type)
      .def("__str__", &WyckoffSite::str)
      .def_readonly("letter", &WyckoffSite::letter)
      .def_readonly("symmetries", &WyckoffSite::symmetries);
}

void export_WallpaperGroup(py::module& m) {
  py::class_<WallpaperGroup> wallpapergroup(m, "WallpaperGroup");
  wallpapergroup
      .def(
          py::init<
              const std::string,
              const std::vector<WyckoffSite>,
              const std::size_t,
              const bool,
              const bool,
              const bool>(),
          py::arg("label"),
          py::arg("wyckoff_sites"),
          py::arg("num_symmetries") = 1,
          py::arg("a_b_equal") = false,
          py::arg("rectangular") = false,
          py::arg("hexagonal") = false)
      .def_readonly("label", &WallpaperGroup::label)
      .def_readonly("wyckoff_sites", &WallpaperGroup::wyckoff_sites)
      .def_readonly("num_symmetries", &WallpaperGroup::num_symmetries)
      .def("num_wyckoffs", &WallpaperGroup::num_wyckoffs);
}
