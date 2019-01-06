/*
 * wallpaper.h
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <vector>

#include <pybind11/pybind11.h>

#include "shapes.h"

#ifndef WALLPAPER_H
#define WALLPAPER_H

/*! \enum mirror
 *
 *  The type of mirror symmetry a SymmetryTransform can have.
 *
 *  The values of the mirror enum resolve to the values of the rotation, so are used for
 *  the calculations involving the rotations.
 */
enum class Mirror {
  m0 = 0,
  m30 = 30,
  m45 = 45,
  m60 = 60,
  m90 = 90,
  m135 = 135,
  m300 = 300,
  m330 = 330,
};

/** \class SymmetryTransform
 *
 * Define the transformations of particle positions for a WyckoffSite.
 *
 * A collection of operators which define the symmetry transformations.
 *
 */
class SymmetryTransform {
public:
  const Vect3 x_coeffs;
  const Vect3 y_coeffs;
  const double rotation_offset;
  const Mirror site_mirror;

  SymmetryTransform(
      const Vect3& x_coeffs,
      const Vect3& y_coeffs,
      const double rotation_offset,
      const Mirror mirror)
      : x_coeffs(x_coeffs), y_coeffs(y_coeffs), rotation_offset(rotation_offset),
        site_mirror(mirror){};

  SymmetryTransform(const SymmetryTransform& other)
      : x_coeffs(other.x_coeffs), y_coeffs(other.y_coeffs),
        rotation_offset(other.rotation_offset), site_mirror(other.site_mirror){};

  SymmetryTransform operator=(const SymmetryTransform& other) {
    return SymmetryTransform(other);
  };
  bool operator==(const SymmetryTransform& other) const;
  friend std::ostream& operator<<(std::ostream&, const SymmetryTransform&);

  Vect2 real_to_fractional(const Vect3& real) const;
  std::string str() const;
};

/** \class WyckoffSite
 *
 * A class which defines a set of Wyckoff parameters.
 *
 * Each set of Wyckoff paramters defines a set of properties, including the positions of
 * particles in Wyckoff Site. A Wyckoff site can have multiple positions, indicated by
 * the muiltiplicity. Each of these positions have some symmetry relationship to each
 * other, either rotaional or mirror.
 *
 * Each Wyckoff instance is labelled by a letter which is an identifier from the
 * Crystallographic Database.
 *
 * A WyckoffSite instance only describes the transformations for taking a position and
 * obtaining the symmetry related positions.
 *
 */
class WyckoffSite {
public:
  const char letter;
  const std::vector<SymmetryTransform> symmetries;
  const std::size_t variability;
  const std::size_t rotations;
  const std::size_t mirrors;

  WyckoffSite(
      const char letter,
      const std::vector<SymmetryTransform>& symmetries,
      const std::size_t variability,
      const std::size_t rotations,
      const std::size_t mirrors)
      : letter(letter), symmetries(symmetries), variability(variability),
        rotations(rotations), mirrors(mirrors){};

  WyckoffSite(const WyckoffSite& site)
      : letter(site.letter),
        symmetries(std::vector<SymmetryTransform>(site.symmetries)),
        variability(site.variability), rotations(site.rotations),
        mirrors(site.mirrors){};

  std::size_t multiplicity() const;
  bool vary_x() const;
  bool vary_y() const;
  int mirror_type() const;
  std::string str() const;

  WyckoffSite operator=(const WyckoffSite& other) {
    return WyckoffSite(other);
  }
  bool operator==(const WyckoffSite& other) const;
};

/** \class IsopointalGroup
 *
 * An IsopointalGroup specifies the occupation of Wyckoff sites for a given number of
 * occupied sites. The entire collection of isopointal groups, is all the unique
 * combinations of occupied sites.
 */
class IsopointalGroup {
public:
  std::vector<WyckoffSite> wyckoff_sites;

  IsopointalGroup(const std::vector<WyckoffSite>& sites) : wyckoff_sites(sites){};

  std::size_t group_multiplicity() const;
  std::string group_string() const;
};

/** \class WallpaperGroup
 *
 * The WallpaperGroup class defines one of the Crystallographic wallpaper groups.
 *
 * The wallpaper group is the highest level description of the symmetry of a crystal
 * structure.
 *
 */
class WallpaperGroup {
public:
  const std::string label;
  const std::vector<WyckoffSite> wyckoff_sites;
  const std::size_t num_symmetries = 0;
  const bool a_b_equal = false;
  const bool rectangular = false;
  const bool hexagonal = false;

  WallpaperGroup(
      const std::string label,
      const std::vector<WyckoffSite>& wyckoffs,
      const std::size_t num_symmetries,
      const bool a_b_equal,
      const bool rectangular,
      const bool hexagonal)
      : label(label), wyckoff_sites(wyckoffs), num_symmetries(num_symmetries),
        a_b_equal(a_b_equal), rectangular(rectangular), hexagonal(hexagonal){};

  WallpaperGroup(const std::string label, const std::vector<WyckoffSite>& wyckoffs)
      : WallpaperGroup(label, wyckoffs, 1, false, false, false){};

  std::size_t num_wyckoffs() {
    return this->wyckoff_sites.size();
  };
};

std::vector<IsopointalGroup> generate_isopointal_groups(
    const Shape& shape,
    const WallpaperGroup& group,
    std::size_t num_occupied_sites);

void export_Mirror(pybind11::module& m);
void export_SymmetryTransform(pybind11::module& m);
void export_WyckoffSite(pybind11::module& m);
void export_WallpaperGroup(pybind11::module& m);

#endif /* !WALLPAPER_H */
