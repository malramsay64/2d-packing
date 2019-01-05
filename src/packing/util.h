/*
 * util.h
 * Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <vector>

#include "basis.h"

#ifndef UTIL_H
#define UTIL_H

std::size_t calculate_shape_replicas(const std::vector<OccupiedSite>& sites);

template <typename T>
std::vector<std::vector<T>> combinations(
    typename std::vector<T>::iterator iter_begin,
    typename std::vector<T>::iterator iter_end,
    std::size_t to_pick);

#endif /* !UTIL_H */
