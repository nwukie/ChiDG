add_subdirectory(elementary)
add_subdirectory(gaussian)
add_subdirectory(trigonometric)
add_subdirectory(polynomial)
add_subdirectory(acoustics)
add_subdirectory(convecting_vortex)
add_subdirectory(scalar_adv_diff_bl_solution)
add_subdirectory(pressure_pulse)
add_subdirectory(isentropic_vortex)

add_lib_sources(
                type_function.f90
                type_function_wrapper.f90
                mod_function.f90
                )
