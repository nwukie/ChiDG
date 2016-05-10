add_test_sources(linear_solvers
                 test_fgmres.pf
                 test_fgmres_ref.pf
                 test_gaussseidel_standard.pf
                 test_fgmres_standard.pf
                 )


add_test_target(linear_solvers unit ${linear_solvers_sources})
