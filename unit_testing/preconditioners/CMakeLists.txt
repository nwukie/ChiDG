set(test_files
                #test_blockjacobi.pf 
                test_ILU0_mod.pf
                )
               
add_test_target(preconditioners unit ${test_files})


