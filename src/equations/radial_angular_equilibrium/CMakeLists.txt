add_subdirectory(bc)

add_lib_sources(
                rae_volume_operator.f90
                rae_boundary_average_operator.f90
                rae_laxfriedrichs_operator.f90
                rae_bc_operator.f90
                rae_volume_source.f90
                model_rae.f90
                eqn_rae.f90
                )
