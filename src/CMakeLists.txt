# include directories
include_directories(${SPECTRUM_SOURCE_DIR}/include)

# set sources
set(SPECTRUM_SOURCES spectrum.c windows.c)

# set libraries
set(SPECTRUM_LIBS m kissfft-float gsl gslcblas)

# add the library
add_library(spectrum SHARED ${SPECTRUM_SOURCES})
target_link_libraries(spectrum ${SPECTRUM_LIBS})
