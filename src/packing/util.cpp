/*
 * util.cpp
 * Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "util.h"

#include <algorithm>

#include "shapes.h"

std::size_t calculate_count_replicas(const std::vector<Site>& occupied_sites) {
  std::size_t replicas{0};
  for (const Site& site : occupied_sites) {
    replicas += site.get_multiplicity();
  }
  return replicas;
}

char compute_chiral_state(const std::vector<Site>& occupied_sites) {
  int chiralsum{0};
  std::size_t totalsum{0};
  for (const Site& site : occupied_sites) {
    chiralsum += site.get_flip_sign() * site.get_multiplicity();
    totalsum += site.get_multiplicity();
  }
  if (chiralsum != 0) {
    return (chiralsum == totalsum) ? 'c' : 's';
  }
  return 'a';
}

template <typename T>
std::vector<std::vector<T>> combinations(
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
    subsets = combinations<T>(index + 1, iter_end, num_picked - 1);
    // Adding the current value to the front of each vector
    for (auto& sub : subsets) {
      sub.insert(sub.begin(), *index);
    }
    // Add the subset for the current value to the end of all values
    combination_list.insert(combination_list.end(), subsets.begin(), subsets.end());
  }
  std::unique(combination_list.begin(), combination_list.end());
  return combination_list;
}
