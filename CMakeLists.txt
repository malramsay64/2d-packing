cmake_minimum_required (VERSION 3.12)

# projectname is the same as the main-executable
project(packing)

include_directories(src/spdlog/include)

add_subdirectory(src/pybind11)

set(SOURCES src/packing/module.cpp
            src/packing/math.cpp
            src/packing/random.cpp
            src/packing/basis.cpp
            src/packing/shapes.cpp
            src/packing/wallpaper.cpp
            )

pybind11_add_module(_packing ${SOURCES})
