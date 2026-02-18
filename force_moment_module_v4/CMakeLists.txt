cmake_minimum_required(VERSION 3.15)
project(force_moment_module LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
enable_testing()

find_package(MPI QUIET COMPONENTS CXX)

find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  find_path(EIGEN3_INCLUDE_DIR
            NAMES Eigen/Dense
            PATH_SUFFIXES eigen3
            PATHS
              $ENV{EIGEN3_INCLUDE_DIR}
              /usr/include
              /usr/local/include
              /opt/homebrew/include
              /opt/include)
  if(NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Eigen3 not found. Install eigen3-devel or set EIGEN3_INCLUDE_DIR to the directory that contains Eigen/Dense.")
  endif()
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
                        INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
endif()

add_library(force_moment_module
  src/aerodynamics.cpp
  src/aero_adapter.cpp
  src/body_importer.cpp
  src/env_config.cpp
  src/obj_loader.cpp
  src/rotate_body.cpp
  src/shadowing.cpp
  src/shuttlecock_aero.cpp
  src/smu_stubs.cpp
  src/vleo_aerodynamics.cpp
  external/tinyobjloader/tiny_obj_loader.cc
)

target_include_directories(force_moment_module
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/external/nlohmann
  PRIVATE
    ${PROJECT_SOURCE_DIR}/external/tinyobjloader
)

target_link_libraries(force_moment_module
  PUBLIC
    Eigen3::Eigen
)

# GCC < 9 requires explicit linkage for std::filesystem.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
  target_link_libraries(force_moment_module PUBLIC stdc++fs)
endif()

add_executable(shuttlecock_demo examples/shuttlecock_demo.cpp)
target_link_libraries(shuttlecock_demo PRIVATE force_moment_module)

add_executable(aero_benchmark benchmarks/aero_benchmark.cpp)
target_link_libraries(aero_benchmark PRIVATE force_moment_module)

add_executable(shadowing_test tests/shadowing_test.cpp)
target_link_libraries(shadowing_test PRIVATE force_moment_module)
target_compile_definitions(shadowing_test PRIVATE FORCE_MOMENT_MODULE_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
add_test(NAME shadowing_test COMMAND shadowing_test)

add_executable(shuttlecock_batch examples/shuttlecock_batch.cpp)
target_link_libraries(shuttlecock_batch PRIVATE force_moment_module)

add_executable(obj_force_batch examples/obj_force_batch.cpp)
target_link_libraries(obj_force_batch PRIVATE force_moment_module)

if(MPI_FOUND)
  add_executable(force_moment_parallel_mpi examples/force_moment_parallel_mpi.cpp)
  target_link_libraries(force_moment_parallel_mpi PRIVATE force_moment_module MPI::MPI_CXX)
else()
  message(WARNING "MPI not found. The force_moment_parallel_mpi target will not be built.")
endif()
