project(test-ccfec)

add_executable(test-ccfec
  runner.cpp
  ReedSolomonFecTest.cpp
)

if (AROS_BUILD)
    aros_add_test(test-ccfec ccfec)
else()
    target_link_libraries(test-ccfec ccfec boost_unit_test_framework)
endif()
