cmake_minimum_required(VERSION 3.1)

project(surf)

# Defining SOURCES as all header and source project files in current directory
FILE(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Adding executable surf with SOURCES files
add_executable(${PROJECT_NAME} ${SOURCES})

# Adding include directories for executable surf
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Adding libraries to be linked with executable surf
target_link_libraries(${PROJECT_NAME})
