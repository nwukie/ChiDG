add_subdirectory(bc)

add_lib_sources(eqn_euler.f90
                euler_volume_operator.f90
                euler_boundary_average_operator.f90
                euler_laxfriedrichs_operator.f90
                euler_roe_operator.f90
                euler_properties.f90
                mod_entropy.f90
                mod_timestep.f90
                )
