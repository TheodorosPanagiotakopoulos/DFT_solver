cmake_minimum_required(VERSION 3.20)
project(dft-fd CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(USE_ARPACK "Enable ARPACK++ for large sparse eigenproblems" ON)
option(USE_KOKKOS "Enable Kokkos for performance portability" ON)

# Dependencies
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(Threads REQUIRED)
find_package(TBB REQUIRED)
find_package(Boost REQUIRED COMPONENTS unit_test_framework math_c99)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)
find_package(Boost REQUIRED)
# Boost.Units doesn't have a separate component; header-only
# Boost.Math is header-first; we use <boost/math/special_functions/cbrt.hpp>

if (USE_KOKKOS)
  find_package(Kokkos REQUIRED)
endif()

if (USE_ARPACK)
  # ARPACK++ is usually header-only + links to ARPACK/BLAS/LAPACK underneath.
  # Provide a hint variable if needed: -DARPACKPP_INCLUDE_DIR=/path/to/include
  if (NOT DEFINED ARPACKPP_INCLUDE_DIR)
    set(ARPACKPP_INCLUDE_DIR "" CACHE PATH "Path to ARPACK++ headers")
  endif()
  if (ARPACKPP_INCLUDE_DIR)
    include_directories(${ARPACKPP_INCLUDE_DIR})
  endif()
  find_library(ARPACK_LIB NAMES arpack arpack64 REQUIRED)
  find_package(LAPACK REQUIRED)
  find_package(BLAS REQUIRED)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(dft-fd
  src/main.cpp
)

target_link_libraries(dft-fd
  PRIVATE
    Eigen3::Eigen
    TBB::tbb
    Threads::Threads
)

if (USE_KOKKOS)
  target_link_libraries(dft-fd PRIVATE Kokkos::kokkos)
  target_compile_definitions(dft-fd PRIVATE DFT_USE_KOKKOS)
endif()

if (USE_ARPACK)
  target_compile_definitions(dft-fd PRIVATE DFT_USE_ARPACKPP)
  target_link_libraries(dft-fd PRIVATE ${ARPACK_LIB} LAPACK::LAPACK BLAS::BLAS)
endif()

# Treat warnings as errors in CI / development
if (NOT MSVC)
  target_compile_options(dft-fd PRIVATE -Wall -Wextra -Wpedantic)
else()
  target_compile_options(dft-fd PRIVATE /W4)
endif()
