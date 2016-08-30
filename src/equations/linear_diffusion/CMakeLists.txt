add_subdirectory(bc)

add_lib_sources(eqn_lineardiffusion.f90
                LD_boundary_diffusive_flux.f90
                LD_volume_diffusive_flux.f90
                LD_volume_source.f90
                LD_properties.f90
                )
