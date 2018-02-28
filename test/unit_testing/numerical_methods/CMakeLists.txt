add_test_sources(numerical_methods
                test_bisect.pf
                test_linspace.pf
                test_inv.pf
                test_determinant.pf
                test_dft.pf
                )
                
add_test_target(numerical_methods unit ${numerical_methods_sources})







