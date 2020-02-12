add_test_sources(functionals
                test_fcl_mod.pf
                test_fcl_group_mod.pf
                test_yforce_mod.pf
                test_xforce_mod.pf
                test_functional_dQ_linearization_mod.pf
                test_functional_dX_linearization__d1e1_mod.pf
                test_functional_dX_linearization__d1e2_mod.pf
                test_functional_dX_linearization__d2e2_mod.pf
                )
                
add_test_target(functionals unit_parallel ${functionals_sources})







