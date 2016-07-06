#add_subdirectory(bc)

add_lib_sources(eqn_burgersinviscid.f90
                BU_boundary_average_advective_flux.f90
                BU_volume_advective_flux.f90
                BU_LaxFriedrichs_flux.f90
                BU_properties.f90
                )
