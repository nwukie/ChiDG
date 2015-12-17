add_test_sources(io
                test_hdfio_grid.pf
                test_hdfio_solution_mod.pf
                test_file_utilities.pf
                )
                
add_test_target(io unit ${io_sources})


# Copy files needed to test io
file(GLOB TEST_FILES "*.h5")
file(COPY ${TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
