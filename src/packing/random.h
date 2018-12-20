#include <pybind11/pybind11.h>
#include <random>

#ifndef FLUKE_H
#define FLUKE_H

double fluke();
void export_fluke(pybind11::module& m);

#endif /* FLUKE_H */
