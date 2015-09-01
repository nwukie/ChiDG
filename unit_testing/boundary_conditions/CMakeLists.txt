set(test_files unit
    test_bcset_mod.pf
    test_bc_mod.pf
    test_bc_periodic_mod.pf
    ) 


add_test_target(boundary_conditions ${test_files})
