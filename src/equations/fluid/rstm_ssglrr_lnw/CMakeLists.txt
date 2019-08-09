add_subdirectory(bc)


add_lib_sources(
                mod_rstm_ssglrrw.f90
                mod_modify_jacobian_rstm.f90
                rstm_ssglrrw_source.f90
                rstm_ssglrrw_blended_source.f90
                rstm_ssglrrw_advection_boundary_average.f90
                rstm_ssglrrw_laxfriedrichs.f90
                rstm_ssglrrw_volume_advection.f90
                rstm_ssglrrw_bc_advection.f90
                rstm_ssglrrw_boundary_diffusion.f90
                rstm_ssglrrw_volume_diffusion.f90
                rstm_ssglrrw_bc_diffusion.f90
                rstm_ssglrrw_artificial_viscosity_operator.f90
                rstm_ssglrrw_artificial_viscosity_bc_operator.f90
                rstm_ssglrrw_artificial_viscosity_boundary_average_operator.f90
                model_rstm_ssglrrw_turbulence_quantities.f90
                model_rstm_ssglrrw_blended_turbulence_quantities.f90
                model_rstm_ssglrrw_production.f90
                model_rstm_ssglrrw_isotropic_dissipation.f90
                model_rstm_ssglrrw_blended_coefficients.f90
                model_rstm_ssglrrw_lrr_coefficients.f90
                model_rstm_ssglrrw_pressure_strain_correlation.f90
                model_rstm_ssglrrw_simple_diffusion.f90
                model_rstm_ssglrrw_generalized_diffusion.f90
                model_rstm_ssglrrw_artificial_viscosity.f90
                model_rstm_ssglrrw_realize_source.f90
                model_rstm_ssglrrw_realizable_stress.f90
                #model_rstm_ssglrrw_reynolds_stress.f90
                )

