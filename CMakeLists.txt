cmake_minimum_required(VERSION 3.10)

project(SPECTRUM LANGUAGES C VERSION 0.1)

set (CMAKE_C_FLAGS "-O3 -ffast-math -Wall")

add_subdirectory(src)
add_subdirectory(examples)
