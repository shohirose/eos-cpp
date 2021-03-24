cmake_minimum_required(VERSION 3.13)

if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
endif()

project(eoscpp CXX)

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)

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
    Boost::headers
  )

enable_testing()
add_subdirectory(test)