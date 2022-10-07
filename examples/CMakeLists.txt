# add executables
add_executable(detect_signal_aic_ed detect_signal_aic_ed.c)
target_link_libraries(detect_signal_aic_ed spectrum)

# add custom target to build all the examples
add_custom_target(examples)
add_dependencies(examples detect_signal_aic_ed)
