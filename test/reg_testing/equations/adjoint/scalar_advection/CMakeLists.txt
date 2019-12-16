add_test_sources(adjoint_equations_advection
                 integration_adjoint_1D_LA_mod.pf
#                 integration_adjoint_1D_LA_ILU0_mod.pf
                 )

add_test_target(adjoint_equations_advection reg_parallel ${adjoint_equations_advection_sources})

