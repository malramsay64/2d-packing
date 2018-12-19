/*
 * module.cpp
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "random.h"
#include <pybind11/pybind11>

PYBIND11_MODULE(_packing, m) { export_fluke(m); }
