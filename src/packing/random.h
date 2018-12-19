#include <pybind11/pybind11.h>
#include <stdio.h>

#ifndef FLUKE_H
#define FLUKE_H

double fluke(void);
void FlukeInit(void);

void export_fluke(pybind11::module& m);

#endif /* FLUKE_H */
