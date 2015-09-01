add_subdirectory(forward_euler)                
add_subdirectory(backward_euler)


add_lib_sources(atype_solver.f90
                atype_solverdata.f90 
                mod_solver.f90
                mod_solverdata.f90
                solverdata_base.f90)

                
