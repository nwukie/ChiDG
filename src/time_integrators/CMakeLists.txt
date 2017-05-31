add_subdirectory(steady)
add_subdirectory(marching)
add_subdirectory(spectral)

add_lib_sources(
                type_time_integrator.f90
                type_system_assembler.f90
                mod_time_integrators.f90
                type_timer.f90
                type_time_manager.f90
                mod_time.f90
                )

                
