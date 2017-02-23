add_subdirectory(bc)

add_lib_sources(
                SA_boundary_average_advective_operator.f90
                SA_volume_advective_operator.f90
                SA_LaxFriedrichs_operator.f90
                SA_bc_operator.f90
                )
