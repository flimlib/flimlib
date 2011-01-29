#
# CMakeLists.txt
#

cmake_minimum_required(VERSION 2.6)

project(slim-curve)

include_directories(include "src/main/c")

file(GLOB_RECURSE SLIM_CURVE_SRC "src/main/c/*.c")

add_library(slim-curve SHARED ${SLIM_CURVE_SRC})
