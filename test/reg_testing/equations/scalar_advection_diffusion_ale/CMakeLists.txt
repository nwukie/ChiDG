add_test_sources(equations_scalar_advection_diffusion_ale
                 integration_1d_linear_adv_diff_ale_mod.pf
                 )


add_test_target(equations_scalar_advection_diffusion_ale reg_parallel ${equations_scalar_advection_diffusion_ale_sources})

