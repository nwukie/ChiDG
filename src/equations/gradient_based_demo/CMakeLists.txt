add_subdirectory(bc)

add_lib_sources(
                eqn_graddemo.f90
                graddemo_P_boundary_operator.f90
                graddemo_P_volume_operator.f90
                graddemo_P_bc_operator.f90
                )
