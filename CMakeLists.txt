cmake_minimum_required(VERSION 3.14)

# Use vcpkg if available
if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
endif()

project(eoscpp CXX)

add_subdirectory(src)

option(EOSCPP_BUILD_TEST "Build unit tests" ON)

if(EOSCPP_BUILD_TEST)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
  )
  # Prevent overriding the parent project's compiler/linker
  # settings on Windows
  if(WIN32)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  endif()
  # Disable unnecessary googletest options
  set(BUILD_GMOCK OFF CACHE BOOL "Build gmock" FORCE)
  set(INSTALL_GTEST OFF CACHE BOOL "Install gtest" FORCE)
  FetchContent_MakeAvailable(googletest)

  # Find GNU Scientific Library
  find_package(GSL REQUIRED)
  
  add_subdirectory(test)
endif()