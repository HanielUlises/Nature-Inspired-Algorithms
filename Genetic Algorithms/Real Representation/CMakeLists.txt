cmake_minimum_required(VERSION 3.10)
project(ContrastEnhancement)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenCV REQUIRED)

include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/PM)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/SBX)

add_executable(ContrastEnhancement
    contrast_enhancement.cpp
    genetic_algo_img.cpp
    PM/polynomial_mut.cpp
    SBX/sbx.cpp
)

target_link_libraries(ContrastEnhancement ${OpenCV_LIBS})