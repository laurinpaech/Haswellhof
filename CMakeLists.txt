cmake_minimum_required(VERSION 3.1)
project(surf)

# Setting Release as build type if not defined otherwise
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# Setting optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -DNDEBUG")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mfma -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# Setting cmake module path for potential later findXXX cmake functionality
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Potential thirdparty libraries
add_subdirectory(ext)

# Our surf implementation
add_subdirectory(src)
