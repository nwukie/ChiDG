add_test_sources(numerical_methods
                test_bisect.pf
                test_linspace.pf
                )
                
add_test_target(numerical_methods unit ${numerical_methods_sources})







