add_subdirectory(bc)

add_lib_sources(
                rac_volume_operator.f90
                rac_boundary_average_operator.f90
                rac_laxfriedrichs_operator.f90
                rac_bc_operator.f90
                rac_volume_source.f90
                model_rac.f90
                eqn_rac.f90
                )
