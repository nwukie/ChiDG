add_subdirectory(bc)

add_lib_sources(eqn_euler.f90
                EULER_volume_advective_flux.f90
                EULER_volume_advective_flux_cacheoptimized.f90
                EULER_boundary_average_advective_flux.f90
                EULER_LaxFriedrichs_flux.f90
                EULER_Roe_flux.f90
                EULER_Roe_flux_cacheoptimized.f90
                EULER_properties.f90
                mod_entropy.f90
                mod_timestep.f90
                )
