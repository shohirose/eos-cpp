add_library(eos INTERFACE
    #lucas_method.cpp
  )
target_compile_features(eos
  INTERFACE
    cxx_std_17
  )
target_compile_options(eos
  INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic>
  )
target_include_directories(eos
  INTERFACE
    #$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    #$<INSTALL_INTERFACE:include>
    ${PROJECT_SOURCE_DIR}/include
  )
target_compile_definitions(eos
  INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:NOMINMAX _USE_MATH_DEFINES>
  )