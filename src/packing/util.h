/*
 * util.h
 * Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <algorithm>
#include <set>
#include <vector>

#include <pybind11/pybind11.h>

#include "basis.h"

#ifndef UTIL_H
#define UTIL_H

template <typename T> void uniqueify(std::vector<T>& v) {
  auto end = v.end();
  for (auto it = v.begin(); it != end; ++it) {
    end = std::remove(it + 1, end, *it);
  }
  v.erase(end, v.end());
}

std::size_t calculate_shape_replicas(const std::vector<OccupiedSite>& sites);

template <typename T>
std::vector<std::vector<T>> combinations_iter(
    typename std::vector<T>::iterator iter_begin,
    typename std::vector<T>::iterator iter_end,
    int num_picked) {
  std::vector<std::vector<T>> combination_list;

  // Stopping condition: Only picking a single site
  if (num_picked == 1) {
    // Convert the 1D input array into a 2D array
    for (auto& index = iter_begin; index != iter_end; iter_begin++) {
      combination_list.push_back(std::vector<T>{*index});
    }
    // Return all the sites as a 2D array
    return combination_list;
  }

  for (auto& index = iter_begin; index != iter_end; index++) {
    // Temporary variable for each loop
    std::vector<std::vector<T>> subsets;
    // Wyckoff Type removed
    subsets = combinations_iter<T>(index + 1, iter_end, num_picked - 1);
    // Adding the current value to the front of the vector
    for (auto& sub : subsets) {
      sub.insert(sub.begin(), *index);
    }
    // Add the subset for the current value to the end of all values
    combination_list.insert(combination_list.end(), subsets.begin(), subsets.end());
  }
  return combination_list;
}

template <typename T>
std::vector<std::vector<T>>
combinations(std::vector<T>& values, std::size_t num_picked) {
  return combinations_iter<T>(values.begin(), values.end(), num_picked);
}

void export_combinations(pybind11::module& m);

#endif /* !UTIL_H */
