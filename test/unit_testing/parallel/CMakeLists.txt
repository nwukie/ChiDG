add_test_sources(parallel
                test_partition_metis.pf
                test_partition_connectivity.pf
                )
                
add_test_target(parallel unit ${parallel_sources})







