cmake_minimum_required (VERSION 2.8)

# projectname is the same as the main-executable
project(packing)

find_package(pybind11 REQUIRED)

set(SOURCES src/packing/module.cpp
            src/packing/random.cpp
            src/packing/shapes.cpp
            src/packing/wallpaper.cpp
            )

pybind11_add_module(_packing, ${SOURCES})
