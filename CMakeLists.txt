cmake_minimum_required(VERSION 3.21)
project(vleo_uq LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(PYBIND11_FINDPYTHON ON)

find_package(pybind11 CONFIG REQUIRED)
# Prefer Eigen's CMake package, but fall back to header-only include discovery
# for HPC/module environments that do not ship Eigen3Config.cmake.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  find_path(EIGEN3_INCLUDE_DIR
    NAMES Eigen/Dense
    HINTS
      ENV EIGEN3_INCLUDE_DIR
      ENV EIGEN3_ROOT
      ENV EBROOTEIGEN
    PATH_SUFFIXES
      include
      include/eigen3
      eigen3
  )
  if(NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR
      "Could not find Eigen3. Set EIGEN3_INCLUDE_DIR (or EIGEN3_ROOT/EBROOTEIGEN) "
      "to a path containing Eigen/Dense.")
  endif()
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
endif()
find_package(OpenMP)

add_subdirectory(force_moment_module_v4)

add_library(vleo_core cpp/src/propagator.cpp)
target_include_directories(vleo_core PUBLIC cpp/include)
target_link_libraries(vleo_core PUBLIC Eigen3::Eigen force_moment_module)
if(OpenMP_CXX_FOUND)
  target_link_libraries(vleo_core PUBLIC OpenMP::OpenMP_CXX)
endif()

pybind11_add_module(_vleo_uq cpp/bindings/bindings.cpp)
target_link_libraries(_vleo_uq PRIVATE vleo_core)
target_include_directories(_vleo_uq PRIVATE cpp/include force_moment_module_v4/include)

include(GNUInstallDirs)
install(TARGETS _vleo_uq
  LIBRARY DESTINATION vleo_uq
  RUNTIME DESTINATION vleo_uq
  ARCHIVE DESTINATION vleo_uq
)
