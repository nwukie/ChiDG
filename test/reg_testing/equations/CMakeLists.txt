#add_subdirectory(linear_advection)
add_subdirectory(euler)


add_test_sources(equations 
                 test_dummy.pf
#                 rtest_linearadvection_mod.pf
#                 rtest_linearadvection_implicit_convergence.pf
#                 rtest_duallinearadvection_mod.pf
                 )

add_test_target(equations reg ${equations_sources} ${test_files})




