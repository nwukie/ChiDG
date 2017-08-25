## Add driver for executable interface
#add_exe_sources(driver.f90)    


## Add wrappers for python interface
add_python_sources(
                   f90wrap_type_chidg.f90.in
                   f90wrap_type_function.f90.in
                   f90wrap_type_equation_set.f90.in
                   f90wrap_mod_equations.f90.in
                   )


