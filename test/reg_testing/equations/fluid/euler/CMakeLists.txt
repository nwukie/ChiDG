add_subdirectory(smoothbump)
add_subdirectory(cylinder)
#add_subdirectory(naca2412)

add_test_sources(equations_fluid_euler)
add_test_target(equations_fluid_euler reg_parallel ${equations_fluid_euler_sources})


