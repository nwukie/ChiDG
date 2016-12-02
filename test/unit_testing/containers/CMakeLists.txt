add_test_sources(containers
                test_densematrix.pf 
                test_densevector.pf 
                test_densevector_operators_mod.pf
                test_blockvector_mod.pf
                test_blockvector_distribute.pf
                test_blockvector_operators_mod.pf
                test_blockmatrix.pf
                test_blockmatrix_mod.pf
                test_dict.pf
                test_operator_block_dot.pf
                test_rvector.pf
                test_ivector.pf
                test_mvector.pf
                test_pvector.pf
                test_chidgVector_operators_mod.pf
                test_chidgMatrix_operator_mv.pf
                test_chidgMatrix_operator_chimera_mv.pf
                test_function_status_mod.pf
                test_cache_data_field_mod.pf
                test_cache_data_mod.pf
#                test_chidg_cache_mod.pf
                )
               
add_test_target(containers unit ${containers_sources})



## Copy files needed to test containers
#file(GLOB TEST_FILES "*.h5")
#file(COPY ${TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
