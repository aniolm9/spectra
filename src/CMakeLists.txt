# set sources
file(GLOB_RECURSE SPECTRUM_SOURCES "*.c")

# set libraries
set(SPECTRUM_LIBS m ${KISSFFT_LIB} gsl gslcblas)

# add the library
add_library(spectrum SHARED ${SPECTRUM_SOURCES})
target_link_libraries(spectrum ${SPECTRUM_LIBS})
