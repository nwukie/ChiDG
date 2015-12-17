add_test_sources(spatial
                test_update_space_linearadvection.pf
                test_update_space_euler.pf
                )
                
add_test_target(spatial unit ${spatial_sources})


# Copy files needed to test io
file(GLOB TEST_FILES "*.h5")
file(COPY ${TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
