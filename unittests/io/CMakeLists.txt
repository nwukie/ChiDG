set(test_files  test_hdfio_read_grid.pf)
                
add_test_target(io ${test_files})


# Copy files needed to test io
file(GLOB TEST_FILES "*.h5")
file(COPY ${TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
