project(example-ccfec)

add_executable(example-ccfec
  example.c
)
target_include_directories(example-ccfec PUBLIC "../include")
target_link_libraries(example-ccfec ccfec)

add_executable(cxx-example-ccfec
  cxx-example.cpp
)


target_include_directories(cxx-example-ccfec PUBLIC "../include")
target_link_libraries(cxx-example-ccfec ccfec)
