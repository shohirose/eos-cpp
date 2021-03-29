find_package(GTest REQUIRED)

function(add_unit_test name)
  add_executable(${name}
    ${name}.cpp
    )
  target_link_libraries(${name}
    PUBLIC
      eos
      GTest::GTest
      GTest::Main
    )
  add_test(
    NAME ${name}
    COMMAND ${name}
    )
endfunction()

add_unit_test(cubic_eos_test)
add_unit_test(vapor_liquid_flash_test)
add_unit_test(lucas_method_test)