add_test_sources(equations_scalar_diffusion
                 integration_1D_linear_diffusion_mod.pf
                 integration_2D_linear_diffusion_mod.pf
                 integration_1D_nonlinear_diffusion_mod.pf
                 )


add_test_target(equations_scalar_diffusion reg_parallel ${equations_scalar_diffusion_sources})

