add_test_sources(chimera
                test__detect_chimera_faces.pf
                test__compute_gq_donor_mod.pf
                test__compute_gq_donor__cylinder_abutting.pf
                test__detect_chimera_donors.pf
                test__compute_chimera_interpolators.pf
                )
                

add_test_target(chimera unit ${chimera_sources})



# Copy files needed to test Chimera
file(GLOB TEST_CHIMERA_FILES "*.h5")
file(COPY ${TEST_CHIMERA_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)


