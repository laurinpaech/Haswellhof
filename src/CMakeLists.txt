cmake_minimum_required(VERSION 3.1)

project(surf)

# Defining SOURCES as all header and source project files in current directory
FILE(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Adding executable surf with SOURCES files
add_executable(surf 
    stb_image.h
    helper.h
    integral_image.h
    integral_image.cpp
    interest_point.h
    response_layer.h
    fasthessian.h
    fasthessian.cpp
    descriptor.h
    descriptor.cpp
    main.cpp
)

# Adding include directories for executable surf
target_include_directories(surf PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Adding libraries to be linked with executable surf
target_link_libraries(surf)

# Adding executable benchmark with SOURCES files
# TODO: (cjancik) if the project structure changes - delete everything in build folder and call "cmake .." again.
add_executable(benchmark 
    stb_image.h
    helper.h
    integral_image.h
    integral_image.cpp
    interest_point.h
    response_layer.h
    fasthessian.h
    fasthessian.cpp
    descriptor.h
    descriptor.cpp
    benchmark_data_to_file.h
    benchmark_data_to_file.cpp
    benchmarking.h
    benchmarking.cpp
    tsc_x86.h
    main_benchmark.cpp
)

# Adding include directories for executable benchmark
target_include_directories(benchmark PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Adding libraries to be linked with executable benchmark
target_link_libraries(benchmark)
