# find include directories
file(GLOB INCLUDE_SUBDIRS ${SPECTRA_INCLUDE_DIR}/*)
list(FILTER INCLUDE_SUBDIRS EXCLUDE REGEX ".h$")

# add executables
add_executable(detect_signal_aic_ed detect_signal_aic_ed.c)
target_link_libraries(detect_signal_aic_ed spectra)
target_include_directories(detect_signal_aic_ed PRIVATE ${SPECTRA_INCLUDE_DIR} ${INCLUDE_SUBDIRS})

# add custom target to build all the examples
add_custom_target(examples)
add_dependencies(examples detect_signal_aic_ed)
