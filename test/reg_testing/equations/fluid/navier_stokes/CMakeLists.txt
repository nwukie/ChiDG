add_test_sources(equations_fluid_navier_stokes
                 integration_force_mod.pf
                 )
add_test_target(equations_fluid_navier_stokes reg_parallel ${equations_fluid_navier_stokes_sources})


