add_subdirectory(grid_generators)
add_subdirectory(duct_acoustics)

add_lib_sources(messenger.f90 
                mod_test_utilities.f90
                mod_string.f90
                mod_condition.f90
                mod_eigenvalues.f90
                mod_eigen.f90
                mod_euler_eigenmodes.f90
                )

