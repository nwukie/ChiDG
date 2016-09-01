add_subdirectory(bc)

add_lib_sources(eqn_linear_advection.f90
                LA_boundary_average_advective_flux.f90
                LA_volume_advective_flux.f90
                LA_LaxFriedrichs_flux.f90
                LA_properties.f90
                )
