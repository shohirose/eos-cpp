cmake_minimum_required(VERSION 3.13)
project(eoscpp CXX)

find_package(Eigen3 REQUIRED)

add_library(eos INTERFACE)
target_compile_features(eos
  INTERFACE
    cxx_std_17
  )
target_compile_options(eos
  INTERFACE
    -Wall
    -Wextra
    -Wpedantic
  )
target_include_directories(eos
  INTERFACE
    include
  )
target_link_libraries(eos
  INTERFACE
    Eigen3::Eigen
  )

enable_testing()
add_subdirectory(test)