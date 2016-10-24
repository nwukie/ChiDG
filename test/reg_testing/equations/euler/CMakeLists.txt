add_subdirectory(smoothbump)
#add_subdirectory(cylinder)
#add_subdirectory(naca2412)

add_test_sources(equations_euler)
add_test_target(equations_euler reg_parallel ${equations_euler_sources})


