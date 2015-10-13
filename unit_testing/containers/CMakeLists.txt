add_test_sources(containers
                test_expansion.pf 
                test_densematrix.pf 
                test_densevector.pf 
                test_densevector_operators_mod.pf
                test_blockvector_mod.pf
                test_blockvector_build.pf
                test_blockvector_distribute.pf
                test_blockvector_operators_mod.pf
                test_blockmatrix.pf
                test_blockmatrix_build.pf
                test_blockmatrix_mod.pf
                test_dict.pf
                test_operator_dot.pf
                test_operator_mv.pf
                test_rvector.pf
                test_ivector.pf
                )
               
add_test_target(containers unit ${containers_sources})



# Copy files needed to test containers
file(GLOB TEST_FILES "*.h5")
file(COPY ${TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
