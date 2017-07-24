add_test_sources(containers
                test_densematrix.pf 
                test_densevector.pf 
                test_densevector_operators_mod.pf
                test_domain_vector_mod.pf
                test_domain_vector_distribute.pf
                test_domain_vector_operators_mod.pf
                test_domain_matrix.pf
                test_domain_matrix_mod.pf
                test_operator_domain_dot.pf
                test_dict.pf
                test_rvector.pf
                test_ivector.pf
                test_mvector.pf
                test_pvector.pf
                test_chidg_vector_operators_mod.pf
                test_chidg_vector_restrict_mod.pf
                test_chidg_matrix_operator_mv.pf
                test_chidg_matrix_operator_chimera_mv.pf
                test_chidg_matrix_operator_HB_mv.pf
                test_function_status_mod.pf
                test_solverdata_auxiliary_fields.pf
                test_cache_data_field_mod.pf
                test_cache_data_mod.pf
                test_densematrix_vector.pf
#                test_chidg_cache_mod.pf
                )
               
add_test_target(containers unit ${containers_sources})

