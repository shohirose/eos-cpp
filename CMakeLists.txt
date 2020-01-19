cmake_minimum_required(VERSION 3.13)
project(eoscpp CXX)

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

enable_testing()
add_subdirectory(test)