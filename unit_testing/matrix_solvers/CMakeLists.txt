set(test_files unit
               test_directsolver.pf
               test_gaussseidel.pf
               test_sor.pf
               test_gmres.pf
               test_gmres_ref.pf
               test_fgmres.pf
               test_fgmres_ref.pf
               )


add_test_target(matrix_solvers ${test_files})
