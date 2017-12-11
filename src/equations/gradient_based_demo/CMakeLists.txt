add_subdirectory(bc)

add_lib_sources(
                graddemo_P_boundary_operator.f90
                graddemo_P_volume_operator.f90
                graddemo_P_bc_operator.f90
                eqn_graddemo.f90
                graddemo_volume_operator.f90
                graddemo_gradp_volume_operator.f90
                )
