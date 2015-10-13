add_test_sources(io
                test_hdfio_grid.pf
                test_hdfio_variable_D1_E1_M1_mod.pf
                test_hdfio_variable_D1_E27_M1_mod.pf)
                
add_test_target(io unit ${io_sources})


# Copy files needed to test io
file(GLOB TEST_FILES "*.h5")
file(COPY ${TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
