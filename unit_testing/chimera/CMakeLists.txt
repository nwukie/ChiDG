add_test_sources(chimera
                test__detect_chimera_interfaces.pf
                )
                

add_test_target(chimera unit ${chimera_sources})



# Copy files needed to test Chimera
file(GLOB TEST_CHIMERA_FILES "*.h5")
file(COPY ${TEST_CHIMERA_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)


