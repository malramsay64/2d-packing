/*
 * random.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "random.h"

#include <random>

std::mt19937_64 generator;
std::uniform_real_distribution<double> distribution(0.0, 1.0);

double fluke() { return distribution(generator); };

void export_fluke(pybind11::module& m) { m.def("fluke", &fluke); };
