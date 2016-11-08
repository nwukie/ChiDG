add_subdirectory(forward_euler)                
add_subdirectory(backward_euler)


add_lib_sources(
                type_time_integrator.f90
                mod_time_integrators.f90
                type_timer.f90
                steady.f90
                )

                
