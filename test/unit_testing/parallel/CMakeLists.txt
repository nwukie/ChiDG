add_test_sources(parallel
                test_partition_metis.pf
                test_partition_connectivity.pf
                test_initialize_communication_mod.pf
                )
                
add_test_target(parallel unit ${parallel_sources})







