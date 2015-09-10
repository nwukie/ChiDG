set(test_files unit
               test_linearadvection_mod.pf
               test_duallinearadvection_mod.pf
               test_euler_mod.pf
               )
add_test_target(equations ${test_files})
