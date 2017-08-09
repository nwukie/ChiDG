add_test_sources(polynomial
               test_legendre.pf 
               test_ordering.pf
               test__project_function_xyz__mod.pf
               test_integrate.pf
               )
                
add_test_target(polynomial unit ${polynomial_sources})
