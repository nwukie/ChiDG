add_subdirectory(bc)


add_lib_sources(
                mod_spalart_allmaras.f90
                spalart_allmaras_source.f90
                spalart_allmaras_advection_boundary_average.f90
                spalart_allmaras_laxfriedrichs.f90
                spalart_allmaras_volume_advection.f90
                spalart_allmaras_bc_advection.f90
                spalart_allmaras_boundary_diffusion.f90
                spalart_allmaras_volume_diffusion.f90
                spalart_allmaras_bc_diffusion.f90
                model_spalart_allmaras_turbulent_model_fields.f90
                )

