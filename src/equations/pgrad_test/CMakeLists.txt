add_subdirectory(bc)

add_lib_sources(
                pgradtest_boundary_operator.f90
                pgradtest_volume_operator.f90
                pgradtest_bc_operator.f90
                eqn_pgradtest.f90
                )
