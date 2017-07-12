add_test_sources(equations_scalar_advection_diffusion
                 integration_1D_linear_adv_diff_mod.pf
                 #radial_equilibrium_mod.pf
                 #integration_2D_linear_advection_diffusion_mod.pf
                 )


add_test_target(equations_scalar_advection_diffusion reg_parallel ${equations_scalar_advection_diffusion_sources})

