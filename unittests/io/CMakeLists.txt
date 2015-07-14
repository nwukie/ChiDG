set(test_files  test_hdfio_read_grid.pf test_hdfio_var_rw_mod.pf)
                
add_test_target(io ${test_files})


# Copy files needed to test io
file(GLOB TEST_FILES "*.h5")
file(COPY ${TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
