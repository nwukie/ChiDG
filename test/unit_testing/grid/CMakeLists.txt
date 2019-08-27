add_test_sources(grid
#                test_point.pf
#                test_compute_element_mappings.pf 
##                test_element__init_geometry.pf
#                test_element_linearA_mod.pf
#                test_element_linearB_mod.pf
#                test_element_linearC_mod.pf
#                test_element_linearD_cylindrical_mod.pf
#                #test_element_ale_pmm_translation_linearA_mod.pf
#                test_element_ale_pmm_dilation_linearA_mod.pf
#                test_domain_linear_mod.pf
#                test_mesh_mod.pf
#                test_chidgdata__add_domain.pf
##                test__grid_operators__mod.pf
#                test_interpolate_mod.pf
#                test_interpolate_general_mod.pf
#                test_interpolate_chimera_mod.pf
#                test_interpolate_chimera_multipledonor_mod.pf
                test_hole_cutting.pf
                )
                
add_test_target(grid unit ${grid_sources})







