cmake_minimum_required (VERSION 3.12)

# projectname is the same as the main-executable
project(packing)

include_directories(src/spdlog/include)

add_subdirectory(src/pybind11)

file(GLOB SOURCES "src/packing/*.cpp")

pybind11_add_module(_packing ${SOURCES})
