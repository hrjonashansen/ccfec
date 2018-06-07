
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")

project(ccfec)

add_library(ccfec STATIC
    encoder.c
    decoder.c
    ff.c
    ff_online.c
)

target_include_directories(ccfec PUBLIC "include")

add_subdirectory(test EXCLUDE_FROM_ALL)
