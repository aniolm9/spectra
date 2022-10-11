# minimum cmake version
cmake_minimum_required(VERSION 3.10)

# project name
project(SPECTRUM LANGUAGES C VERSION 0.1)

# set compiler flags, mainly optimization
set(CMAKE_C_FLAGS "-O3 -ffast-math -Wall")

# include directories
include_directories(${SPECTRUM_SOURCE_DIR}/include)

# find kissfft library
find_library(KISSFFT_LIB kissfft-float kiss_fft_float HINTS ${SPECTRUM_SOURCE_DIR}/lib)

# load CMakeLists from other subdirectories
add_subdirectory(src)
add_subdirectory(examples)
