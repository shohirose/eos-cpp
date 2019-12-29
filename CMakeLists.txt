cmake_minimum_required(VERSION 3.8)
project(eoscpp CXX)

add_library(eos INTERFACE)
target_compile_features(eos
  INTERFACE
    cxx_std_11
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

add_subdirectory(test)