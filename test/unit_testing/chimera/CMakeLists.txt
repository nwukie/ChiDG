add_test_sources(chimera
                test__detect_chimera_faces.pf
                test__find_gq_donor_mod.pf
                test__find_gq_donor__cylinder_abutting.pf
                test__detect_chimera_donors.pf
                test__compute_chimera_interpolators.pf
                )
                

add_test_target(chimera unit ${chimera_sources})

