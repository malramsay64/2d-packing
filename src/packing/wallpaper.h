/*
 * wallpaper.h
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <vector>

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
enum mirror {
  m0 = 0,
  m90 = 90,
  m45 = 45,
  m135 = 135,
  m30 = 30,
  m60 = 60,
  m330 = 330,
  m300 = 300
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
  SymmetryTransform(Vect3 x_coeffs, Vect3 y_coeffs, double rotation_offset);
  Vect3 x_coeffs;
  Vect3 y_coeffs;
  double rotation_offset;
  enum mirror site_mirror;

  bool operator==(const SymmetryTransform& other) const;
  friend std::ostream& operator<<(std::ostream&, const SymmetryTransform&);

  Vect2 real_to_fractional(const Vect3& real) const;
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
  char letter;
  int variability;
  int rotations;
  int mirrors;
  std::vector<SymmetryTransform> symmetries;

  std::size_t multiplicity() const;
  bool vary_x() const;
  bool vary_y() const;
  int mirror_type() const;
  std::string str() const;

  bool operator==(const WyckoffSite& other) const;
};

/** \class IsopointalGroup
 *
 * An IsopointalGroup is a collection
 */
class IsopointalGroup {
public:
  std::vector<WyckoffSite> wyckoff_sites;

  IsopointalGroup(std::vector<WyckoffSite>& sites) : wyckoff_sites(sites){};

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
  WallpaperGroup(std::string label, std::vector<WyckoffSite> wyckoffs);
  const std::string label;
  const bool a_b_equal = false;
  const bool hexagonal = false;
  const bool rectangular = false;
  int num_symmetries = 0;
  std::vector<WyckoffSite> wyckoff_sites;
  int num_wyckoffs;
};

std::vector<IsopointalGroup> generate_isopointal_groups(
    const Shape& shape,
    const WallpaperGroup& group,
    std::size_t num_occupied_sites);

#endif /* !WALLPAPER_H */
