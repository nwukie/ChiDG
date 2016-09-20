add_subdirectory(bc)

add_lib_sources(eqn_scalar_diffusion.f90
                SD_boundary_operator.f90
                SD_volume_operator.f90
                SD_volume_source.f90
                SD_bc_operator.f90
                SD_properties.f90
                )
