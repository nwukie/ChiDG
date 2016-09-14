add_subdirectory(bc)

add_lib_sources(eqn_linear_diffusion.f90
                LD_boundary_diffusive_operator.f90
                LD_volume_diffusive_operator.f90
                LD_volume_source.f90
                LD_bc_operator.f90
                LD_properties.f90
                )
