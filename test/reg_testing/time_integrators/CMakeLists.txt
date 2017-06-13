add_test_sources(time_integrators
                 integration_test_DIRK_mod.pf
                 )


add_test_target(time_integrators reg_parallel ${time_integrators_sources})

