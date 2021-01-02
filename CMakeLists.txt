cmake_minimum_required(VERSION 3.18.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/usr/bin/g++-9")

project(stiff)

add_subdirectory(src)
add_subdirectory(includes)

find_package(Armadillo REQUIRED)

include_directories(
    ${ARMADILLO_INCLUDE_DIRS}
    includes
)

add_executable(stiff main.cpp)

target_link_libraries(
    stiff
    comp_lib
    ${ARMADILLO_LIBRARIES})