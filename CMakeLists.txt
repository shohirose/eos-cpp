cmake_minimum_required(VERSION 3.13)

# Use vcpkg if available
if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
endif()

project(eoscpp CXX)

find_package(GSL REQUIRED)
find_package(Microsoft.GSL CONFIG REQUIRED)
add_subdirectory(src)

option(EOSCPP_BUILD_TEST "Build unit tests" ON)

if(EOSCPP_BUILD_TEST)
  # Disable googletest options
  option(BUILD_GMOCK "Build gmock" OFF)
  option(INSTALL_GTEST "Install gtest" OFF)

  # Prevent overriding the parent project's compiler/linker
  # settings on Windows
  if(WIN32)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  endif()

  add_subdirectory(third-party/googletest)

  enable_testing()
  add_subdirectory(test)
endif()