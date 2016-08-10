add_subdirectory(smoothbump)

add_test_sources(equations_parallel)
add_test_target(equations_parallel reg_parallel ${equations_parallel_sources})


