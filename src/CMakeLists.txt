add_library(eos
    lucas_method.cpp
    polynomial.cpp
  )
target_compile_features(eos
  PUBLIC
    cxx_std_17
  )
target_compile_options(eos
  PUBLIC
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic>
  )
target_include_directories(eos
  PUBLIC
    #$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    #$<INSTALL_INTERFACE:include>
    ${PROJECT_SOURCE_DIR}/include
  )
target_link_libraries(eos
  PUBLIC
    Microsoft.GSL::GSL
  PRIVATE
    GSL::gsl
  )
target_compile_definitions(eos
  PUBLIC
    $<$<CXX_COMPILER_ID:MSVC>:NOMINMAX _USE_MATH_DEFINES>
  )