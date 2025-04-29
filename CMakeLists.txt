cmake_minimum_required(VERSION 3.10)

project(SmithWaterman LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SOURCES
    "${CMAKE_SOURCE_DIR}/*.cpp"
    "${CMAKE_SOURCE_DIR}/alignment/*.cpp"
)
include_directories(${CMAKE_SOURCE_DIR})

add_executable(main ${SOURCES})

target_include_directories(main PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external/xsimd/include
    ${CMAKE_SOURCE_DIR}/alignment
)

target_compile_options(main PRIVATE -O3 -march=native)