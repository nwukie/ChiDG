add_test_sources(equations
                test_equation_dummy.pf
                test_equationset.pf
                test_linearadvection_mod.pf
                test_duallinearadvection_mod.pf
                test_euler_mod.pf
                )

add_test_target(equations unit ${equations_sources})
