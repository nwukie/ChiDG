#add_subdirectory(forward_euler)                
#add_subdirectory(backward_euler)
#add_subdirectory(quasi_newton)
add_subdirectory(newton)


add_lib_sources(type_timescheme.f90
                mod_timescheme.f90
                type_timer.f90
                )

                
