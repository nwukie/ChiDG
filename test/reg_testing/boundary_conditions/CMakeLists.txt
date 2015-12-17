add_test_sources(boundary_conditions
                 rtest_periodic_boundary_connectivity.pf
                 )
add_test_target(boundary_conditions reg ${boundary_conditions_sources})
