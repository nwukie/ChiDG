set(test_files test_densematrix.pf 
               test_densevector.pf 
               test_expansion.pf 
               test_blockmatrix.pf
               test_blockmatrix_mod.pf
               test_dict.pf)
               
add_test_target(containers ${test_files})



# Copy files needed to test containers
file(GLOB TEST_FILES "*.h5")
file(COPY ${TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
