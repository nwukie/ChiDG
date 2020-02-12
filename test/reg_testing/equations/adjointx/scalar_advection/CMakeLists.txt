add_test_sources(adjointx_equations_advection
                 integration_adjointx_1D_LA_mod.pf
                 )

add_test_target(adjointx_equations_advection reg_parallel ${adjointx_equations_advection_sources})

