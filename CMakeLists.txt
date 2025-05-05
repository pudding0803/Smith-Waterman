cmake_minimum_required(VERSION 3.18)
project(SmithWaterman LANGUAGES CXX CUDA)

enable_language(CUDA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

file(GLOB CPP_SOURCES
    "${CMAKE_SOURCE_DIR}/*.cpp"
    "${CMAKE_SOURCE_DIR}/alignment/*.cpp"
)

file(GLOB CUDA_SOURCES
    "${CMAKE_SOURCE_DIR}/alignment/*.cu"
)

add_executable(main ${CPP_SOURCES} ${CUDA_SOURCES})

set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE CUDA)

target_include_directories(main PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external/xsimd/include
    ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    ${CMAKE_SOURCE_DIR}/alignment
)

target_compile_options(main PRIVATE
    -O3
    -march=native
)

set_target_properties(main PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES 86
)
