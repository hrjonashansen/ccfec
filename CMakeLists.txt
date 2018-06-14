project(ccfec)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")

if (NOT AROS_BUILD)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -flto -Wall -Wextra -pedantic")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -march=native -flto -Wall -Wextra -pedantic")
endif()

add_library(ccfec STATIC
    encoder.c
    decoder.c
    ff.c
    ff_basic.c
)

target_include_directories(ccfec PUBLIC "include")

add_subdirectory(test EXCLUDE_FROM_ALL)

if (NOT AROS_BUILD)
    add_subdirectory(example)
endif()
