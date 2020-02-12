add_test_sources(adjointx_equations_advection_diffusion
                 integration_adjointx_1D_LAD_mod.pf
                 integration_adjointx_1D_LAD_dual_mod.pf
                 )

add_test_target(adjointx_equations_advection_diffusion reg_parallel ${adjointx_equations_advection_diffusion_sources})

