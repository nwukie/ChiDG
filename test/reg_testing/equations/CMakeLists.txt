add_subdirectory(linear_advection)
add_subdirectory(euler)



# Propagate sources and update paths from previous subdirectories
add_test_sources(equations)



add_test_target(equations reg ${equations_sources} ${test_files})




