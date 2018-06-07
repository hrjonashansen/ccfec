project(test-ccfec)

add_executable(test-ccfec
  runner.cpp
  ReedSolomonFecTest.cpp
)

if (AROS_BUILD)
    aros_add_test(test-ccfec ccfec)
endif()
