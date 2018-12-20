add_subdirectory(bc)


add_lib_sources(
                mod_sst.f90
                mod_modify_jacobian_sst.f90
                sst_source.f90
                sst_roe_operator.f90
                sst_advection_boundary_average.f90
                sst_laxfriedrichs.f90
                sst_volume_advection.f90
                sst_bc_advection.f90
                sst_boundary_diffusion.f90
                sst_volume_diffusion.f90
                sst_bc_diffusion.f90
                sst_artificial_viscosity_operator.f90
                sst_artificial_viscosity_bc_operator.f90
                sst_artificial_viscosity_boundary_average_operator.f90
                model_sst_turbulence_quantities.f90
                model_sst_turbulence_kinetic_energy.f90
                model_sst_blended_coefficients.f90
                model_sst_source_terms.f90
                )

