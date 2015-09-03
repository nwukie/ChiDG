add_subdirectory(forward_euler)                
add_subdirectory(backward_euler)


add_lib_sources(atype_time_scheme.f90
                atype_solverdata.f90 
                mod_time_scheme.f90
                mod_solverdata.f90
                solverdata_base.f90)

                
