set(test_files unit
               test_directsolver.pf
               test_gaussseidel.pf
               )


add_test_target(matrix_solvers ${test_files})
