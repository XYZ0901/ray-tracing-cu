cmake_minimum_required(VERSION 3.16)

set(CMAKE_CUDA_ARCHITECTURES 86)

project(rt-cu CUDA CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenCV REQUIRED)

message(STATUS "OpenCV library status:")
message(STATUS "    version: ${OpenCV_VERSION}")
message(STATUS "    libraries: ${OpenCV_LIBS}")
message(STATUS "    include path: ${OpenCV_INCLUDE_DIRS}")

include_directories(${OpenCV_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR}/heads)

add_executable(rt-cu main.cu)

set_target_properties(rt-cu PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(rt-cu ${OpenCV_LIBS})
