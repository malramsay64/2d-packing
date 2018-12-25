/*
 * random.h
 * Copyright (C) 2018 Malcolm Ramsay <malramsay64@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */
#include <pybind11/pybind11.h>

#ifndef FLUKE_H
#define FLUKE_H

double fluke();
void export_fluke(pybind11::module& m);

#endif /* FLUKE_H */
