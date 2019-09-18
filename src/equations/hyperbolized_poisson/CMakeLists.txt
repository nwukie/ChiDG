add_subdirectory(bc)

add_lib_sources(
                HP_boundary_average.f90
                HP_volume.f90
                HP_LaxFriedrichs.f90
                HP_bc.f90
                eqn_hyperbolized_poisson.f90
                )
