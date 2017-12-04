add_test_sources(preconditioners
                test_preconditioner_dummy.pf
                #test_blockjacobi.pf 
                test_ILU0_mod.pf
                test_ILU0_HB_mod.pf
                )
               
add_test_target(preconditioners unit ${preconditioners_sources})


