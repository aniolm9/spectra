# set sources
file(GLOB_RECURSE SPECTRUM_SOURCES "*.c")

# set libraries
set(SPECTRUM_LIBS m ${KISSFFT_LIB} gsl gslcblas)

# find include directories
file(GLOB INCLUDE_SUBDIRS ${SPECTRUM_INCLUDE_DIR}/*)
list(FILTER INCLUDE_SUBDIRS EXCLUDE REGEX ".h$")

# add the library
add_library(spectrum SHARED ${SPECTRUM_SOURCES})
target_link_libraries(spectrum ${SPECTRUM_LIBS})
target_include_directories(spectrum PRIVATE ${SPECTRUM_INCLUDE_DIR} ${INCLUDE_SUBDIRS})
