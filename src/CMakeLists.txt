cmake_minimum_required(VERSION 3.1)

project(surf)

# Defining SOURCES as all header and source project files in current directory
FILE(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Adding executable surf
add_executable(surf 
    base/stb_image.h
    base/helper.h
    base/helper.cpp
    base/integral_image.h
    base/integral_image.cpp
    base/interest_point.h
    base/response_layer.h
    base/fasthessian.h
    base/fasthessian.cpp
    base/descriptor.h
    base/descriptor.cpp
    base/main.cpp
)

# Adding include directories for executable surf
target_include_directories(surf 
    PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/base
)

# Adding libraries to be linked with executable surf
target_link_libraries(surf)

# Adding executable benchmark
# TODO: (cjancik) if the project structure changes - delete everything in build folder and call "cmake .." again.
add_executable(benchmark 
    base/stb_image.h
    base/helper.h
    base/helper.cpp
    base/integral_image.h
    base/integral_image.cpp
    base/interest_point.h
    base/response_layer.h
    base/fasthessian.h
    base/fasthessian.cpp
    base/descriptor.h
    base/descriptor.cpp
    benchmark/benchmark_data_to_file.h
    benchmark/benchmark_data_to_file.cpp
    benchmark/benchmarking.h
    benchmark/benchmarking.cpp
    benchmark/tsc_x86.h
    benchmark/main_benchmark.cpp
    # Add optimization header and source files here
    optimization/integral_image_opt.h
    optimization/integral_image_opt.cpp
    optimization/integral_image_simd.h
    optimization/integral_image_simd.cpp
    optimization/fasthessian_opt.h
    optimization/fasthessian_opt.cpp
    optimization/descriptor_opt.h
    optimization/descriptor_opt.cpp
    optimization/descriptors_wrapper.h
    optimization/descriptors_wrapper.cpp
)

# Adding include directories for executable benchmark
target_include_directories(benchmark 
    PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/base
    ${CMAKE_CURRENT_SOURCE_DIR}/benchmark 
    ${CMAKE_CURRENT_SOURCE_DIR}/optimization
)

# Adding libraries to be linked with executable benchmark
target_link_libraries(benchmark)

# Adding executable optimization
add_executable(optimization 
    base/stb_image.h
    base/helper.h
    base/helper.cpp
    base/integral_image.h
    base/integral_image.cpp
    base/interest_point.h
    base/response_layer.h
    base/fasthessian.h
    base/fasthessian.cpp
    base/descriptor.h
    base/descriptor.cpp
    benchmark/benchmark_data_to_file.h
    benchmark/benchmark_data_to_file.cpp
    benchmark/benchmarking.h
    benchmark/benchmarking.cpp
    benchmark/tsc_x86.h
    # Add optimization header and source files here
    optimization/integral_image_opt.h
    optimization/integral_image_opt.cpp
    optimization/integral_image_simd.h
    optimization/integral_image_simd.cpp
    optimization/fasthessian_opt.h
    optimization/fasthessian_opt.cpp
    optimization/main_opt.cpp
    optimization/descriptor_opt.h
    optimization/descriptor_opt.cpp
    optimization/descriptors_wrapper.h
    optimization/descriptors_wrapper.cpp
)

# Adding include directories for executable optimization
target_include_directories(optimization 
    PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/base 
    ${CMAKE_CURRENT_SOURCE_DIR}/benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/optimization
)

# Adding libraries to be linked with executable optimization
target_link_libraries(optimization)

# Adding executable validation
add_executable(validation 
    base/stb_image.h
    base/helper.h
    base/helper.cpp
    base/integral_image.h
    base/integral_image.cpp
    base/interest_point.h
    base/response_layer.h
    base/fasthessian.h
    base/fasthessian.cpp
    base/descriptor.h
    base/descriptor.cpp
    benchmark/benchmark_data_to_file.h
    benchmark/benchmark_data_to_file.cpp
    benchmark/benchmarking.h
    benchmark/benchmarking.cpp
    benchmark/tsc_x86.h
    # Add optimization header and source files here
    optimization/integral_image_opt.h
    optimization/integral_image_opt.cpp
    optimization/integral_image_simd.h
    optimization/integral_image_simd.cpp
    optimization/fasthessian_opt.h
    optimization/fasthessian_opt.cpp
    optimization/descriptor_opt.h
    optimization/descriptor_opt.cpp
    optimization/descriptors_wrapper.h
    optimization/descriptors_wrapper.cpp
    #autotune/descriptor_opt_gen.h
    #autotune/descriptor_opt_gen.cpp
    #autotune/fasthessian_opt_gen.h
    #autotune/fasthessian_opt_gen.cpp
    # Add validation header and source files here
    validation/validation.h
    validation/validation.cpp
    validation/validation_integral_image.h
    validation/validation_integral_image.cpp
    validation/main_validation.cpp
)

# Adding include directories for executable validation
target_include_directories(validation 
    PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/base 
    ${CMAKE_CURRENT_SOURCE_DIR}/benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/optimization
    ${CMAKE_CURRENT_SOURCE_DIR}/validation
    ${CMAKE_CURRENT_SOURCE_DIR}/autotune
)

# Adding libraries to be linked with executable validation
target_link_libraries(validation)

# Adding executable playground
add_executable(playground 
    base/stb_image.h
    base/helper.h
    base/helper.cpp
    base/integral_image.h
    base/integral_image.cpp
    base/interest_point.h
    base/response_layer.h
    base/fasthessian.h
    base/fasthessian.cpp
    base/descriptor.h
    base/descriptor.cpp
    benchmark/benchmark_data_to_file.h
    benchmark/benchmark_data_to_file.cpp
    benchmark/benchmarking.h
    benchmark/benchmarking.cpp
    benchmark/tsc_x86.h
    # Add optimization header and source files here
    optimization/integral_image_opt.h
    optimization/integral_image_opt.cpp
    optimization/integral_image_simd.h
    optimization/integral_image_simd.cpp
    optimization/fasthessian_opt.h
    optimization/fasthessian_opt.cpp
    optimization/descriptor_opt.h
    optimization/descriptor_opt.cpp
    optimization/descriptors_wrapper.h
    optimization/descriptors_wrapper.cpp
    #autotune/descriptor_opt_gen.h
    #autotune/descriptor_opt_gen.cpp
    # Add validation header and source files here
    validation/validation.h
    validation/validation.cpp
    validation/validation_integral_image.h
    validation/validation_integral_image.cpp
    # Add playground header and source files here
    playground/playground.h
    playground/playground.cpp
    playground/main_playground.cpp
)

# Adding include directories for executable playground
target_include_directories(playground 
    PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/base 
    ${CMAKE_CURRENT_SOURCE_DIR}/benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/optimization
    ${CMAKE_CURRENT_SOURCE_DIR}/validation
    ${CMAKE_CURRENT_SOURCE_DIR}/playground
)

# Adding libraries to be linked with executable playground
target_link_libraries(playground)

# Adding executable autotune
add_executable(autotune 
    base/stb_image.h
    base/helper.h
    base/helper.cpp
    base/integral_image.h
    base/integral_image.cpp
    base/interest_point.h
    base/response_layer.h
    base/fasthessian.h
    base/fasthessian.cpp
    base/descriptor.h
    base/descriptor.cpp
    benchmark/benchmark_data_to_file.h
    benchmark/benchmark_data_to_file.cpp
    benchmark/benchmarking.h
    benchmark/benchmarking.cpp
    benchmark/tsc_x86.h
    # Add optimization header and source files here
    optimization/fasthessian_opt.h
    optimization/fasthessian_opt.cpp
    optimization/descriptor_opt.h
    optimization/descriptor_opt.cpp
    optimization/descriptors_wrapper.h
    optimization/descriptors_wrapper.cpp
    optimization/integral_image_opt.h
    optimization/integral_image_opt.cpp
    # Add autotune header and source files here
    autotune/descriptor_opt_gen.h
    autotune/descriptor_opt_gen.cpp
    #autotune/fasthessian_opt_gen.h
    #autotune/fasthessian_opt_gen.cpp
    autotune/main_autotune.cpp
)

# Adding include directories for executable autotune
target_include_directories(autotune 
    PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/base 
    ${CMAKE_CURRENT_SOURCE_DIR}/benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/optimization
    ${CMAKE_CURRENT_SOURCE_DIR}/validation
    ${CMAKE_CURRENT_SOURCE_DIR}/autotune
)

# Adding libraries to be linked with executable autotune
target_link_libraries(autotune)