set(test_files test_Legendre_1D.pf 
               test_Legendre_2D.pf
               test_Ordering.pf
               test__project_function_xyz__mod.pf)
add_test_target(polynomial ${test_files})
