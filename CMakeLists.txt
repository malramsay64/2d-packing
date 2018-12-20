cmake_minimum_required (VERSION 3.12)

# projectname is the same as the main-executable
project(packing)

find_package(Python3 COMPONENTS Interpreter REQUIRED)

message("Found python3 executable: ${Python3_EXECUTABLE}")

find_package(pybind11 REQUIRED)

set(SOURCES src/packing/module.cpp
            src/packing/math.cpp
            src/packing/random.cpp
            src/packing/shapes.cpp
            src/packing/wallpaper.cpp
            )

pybind11_add_module(_packing ${SOURCES})
