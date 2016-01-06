add_test_sources(matrixsolvers
#                 test_directsolver.pf
#                 test_gaussseidel.pf
#                 test_sor.pf
#                 test_gmres.pf
#                 test_gmres_ref.pf
                 test_fgmres.pf
                 test_fgmres_ref.pf
                 test_gaussseidel_standard.pf
                 test_fgmres_standard.pf
                 )


add_test_target(matrixsolvers unit ${matrixsolvers_sources})
