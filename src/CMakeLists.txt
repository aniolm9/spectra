# set sources
file(GLOB_RECURSE SPECTRA_SOURCES "*.c")

# set libraries
set(SPECTRA_LIBS m ${KISSFFT_LIB} gsl gslcblas)

# find include directories
file(GLOB INCLUDE_SUBDIRS ${SPECTRA_INCLUDE_DIR}/*)
list(FILTER INCLUDE_SUBDIRS EXCLUDE REGEX ".h$")

# add the library
add_library(spectra SHARED ${SPECTRA_SOURCES})
target_link_libraries(spectra ${SPECTRA_LIBS})
target_include_directories(spectra PRIVATE ${SPECTRA_INCLUDE_DIR} ${INCLUDE_SUBDIRS})
