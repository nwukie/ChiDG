add_test_sources(adjoint_equations_advection_diffusion
                 integration_adjoint_1D_LAD_mod.pf
#                 integration_adjoint_1D_LAD_ILU0_mod.pf
                 )


add_test_target(adjoint_equations_advection_diffusion reg_parallel ${adjoint_equations_advection_diffusion_sources})

