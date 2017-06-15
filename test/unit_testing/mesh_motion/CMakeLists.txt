add_test_sources(mesh_motion
                 test_pmm_mod.pf
                 ) 


add_test_target(mesh_motion unit ${mesh_motion_sources})
