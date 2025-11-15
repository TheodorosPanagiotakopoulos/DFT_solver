cmake_minimum_required(VERSION 3.14)
project(dft-fd CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(USE_ARPACK "Enable ARPACK++ for large sparse eigenproblems" ON)
option(USE_KOKKOS "Enable Kokkos for performance portability" ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(TBB REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost REQUIRED COMPONENTS math_c99)

if(USE_KOKKOS)
  find_package(Kokkos REQUIRED)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(dft-fd
  src/main.cpp
)

target_link_libraries(dft-fd
  PRIVATE Eigen3::Eigen
  PRIVATE TBB::tbb
  PRIVATE Threads::Threads
)

if(USE_KOKKOS)
  target_compile_definitions(dft-fd PRIVATE DFT_USE_KOKKOS)
  target_link_libraries(dft-fd PRIVATE Kokkos::kokkos)
endif()

if(USE_ARPACK)
  target_compile_definitions(dft-fd PRIVATE DFT_USE_ARPACKPP)
endif()

if(NOT MSVC)
  target_compile_options(dft-fd PRIVATE -Wall -Wextra -Wpedantic)
else()
  target_compile_options(dft-fd PRIVATE /W4)
endif()
.
