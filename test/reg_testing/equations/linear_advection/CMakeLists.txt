add_test_sources(equations_serial
                 rtest_linearadvection_implicit_convergence.pf
                 rtest_linearadvection_mod.pf
                 rtest_duallinearadvection_mod.pf
                 )
add_test_target(equations_serial reg_parallel ${equations_serial_sources})


