cmake_minimum_required(VERSION 3.14)

# Use vcpkg if available
if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
endif()

project(eoscpp CXX)

# =============================================================================
#                                   Library
# =============================================================================
add_library(eos INTERFACE)
target_compile_features(eos INTERFACE cxx_std_17)
target_compile_options(eos
  INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic>
  )
target_include_directories(eos
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  )
target_compile_definitions(eos
  INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:NOMINMAX _USE_MATH_DEFINES>
  )

# Use precompiled headers
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16")
  target_precompile_headers(eos
    INTERFACE
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/eos/CubicEosBase.hpp>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/eos/VanDerWaalsEos.hpp>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/eos/SoaveRedlichKwongEos.hpp>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/eos/PengRobinsonEos.hpp>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/eos/IdentityCorrectionFactor.hpp>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/eos/SoaveCorrectionFactor.hpp>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/eos/MathematicalConstants.hpp>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/eos/ThermodynamicConstants.hpp>
    )
endif()
add_library(eos::eos ALIAS eos)

# Install headers
install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/include/eos
  DESTINATION include
  FILES_MATCHING PATTERN "*.hpp"
  )
# Requires CMake 3.14 or higher
include(GNUInstallDirs)
# Generate eos-config.cmake
install(TARGETS eos
  EXPORT eos-config
  )
# Install eos-config.cmake
install(EXPORT eos-config
  NAMESPACE eos::
  DESTINATION lib/cmake/eos
  )

# =============================================================================
#                                  Unit Tests
# =============================================================================
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