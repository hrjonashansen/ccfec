project(test-ccfec)

add_executable(test-ccfec
  runner.cpp
  ReedSolomonFecTest.cpp
)

if (AROS_BUILD)
    aros_add_test(test-ccfec ccfec)
else()
    target_link_libraries(test-ccfec ccfec boost_unit_test_framework)

    if (CCFEC_TEST_ALL_FF_IMPLEMENTATIONS)
        add_executable(test-ccfec-online
          runner.cpp
          ReedSolomonFecTest.cpp
        )
        target_link_libraries(test-ccfec-online ccfec-online boost_unit_test_framework)

        add_executable(test-ccfec-simple
          runner.cpp
          ReedSolomonFecTest.cpp
        )
        target_link_libraries(test-ccfec-simple ccfec-simple boost_unit_test_framework)

        if (SSSE3_INSTRUCTIONS_ENABLED)
            add_executable(test-ccfec-ssse3
              runner.cpp
              ReedSolomonFecTest.cpp
            )
            target_link_libraries(test-ccfec-ssse3 ccfec-ssse3 boost_unit_test_framework)
        endif()

        if (AVX2_INSTRUCTIONS_ENABLED)
            add_executable(test-ccfec-avx2
              runner.cpp
              ReedSolomonFecTest.cpp
            )
            target_link_libraries(test-ccfec-avx2 ccfec-avx2 boost_unit_test_framework)
        endif()

        if (NEON_INSTRUCTIONS_ENABLED)
            add_executable(test-ccfec-neon
              runner.cpp
              ReedSolomonFecTest.cpp
            )
            target_link_libraries(test-ccfec-neon ccfec-neon boost_unit_test_framework)
        endif()
    endif() # CCFEC_TEST_ALL_FF_IMPLEMENTATIONS
endif()
