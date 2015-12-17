add_test_sources(polynomial
               test_Legendre_1D.pf 
               test_Legendre_2D.pf
               test_Ordering.pf
               test__project_function_xyz__mod.pf
               test_integrate.pf
               )
                
add_test_target(polynomial unit ${polynomial_sources})
