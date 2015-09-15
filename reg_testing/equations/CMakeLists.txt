set(test_files reg
               rtest_linearadvection_mod.pf
               rtest_linearadvection_implicit_convergence.pf
               rtest_duallinearadvection_mod.pf
               )
add_test_target(equations ${test_files})
