add_subdirectory(bc)

add_lib_sources(
                tm_volume_operator.f90
                tm_boundary_average_operator.f90
                tm_laxfriedrichs_operator.f90
                tm_bc_operator.f90
                tm_volume_source.f90
                model_tm.f90
                eqn_tm.f90
                )
