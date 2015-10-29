add_test_sources(equations
                 rtest_smoothbump_rotations_mod.pf
                 )




# Copy files needed to test equation
file(GLOB TEST_FILES "*.h5")
file(COPY ${TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)



