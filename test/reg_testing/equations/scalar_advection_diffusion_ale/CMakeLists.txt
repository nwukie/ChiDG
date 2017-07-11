add_test_sources(equations_scalar_advection_diffusion_ale
                 integration_1D_linear_adv_diff_ale_mod.pf
                 #integration_2D_linear_advection_diffusion_mod.pf
                 )


add_test_target(equations_scalar_advection_diffusion_ale reg_parallel ${equations_scalar_advection_diffusion_ale_sources})

