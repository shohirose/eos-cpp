find_package(GTest REQUIRED)

add_executable(eos_test
  eos_test.cpp
  )
target_link_libraries(eos_test
  PUBLIC
    eos
    GTest::GTest
    GTest::Main
  )
add_test(
  NAME eos_test
  COMMAND eos_test
  )