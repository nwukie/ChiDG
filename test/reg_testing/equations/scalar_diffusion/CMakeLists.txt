add_test_sources(equations_diffusion
                 integration_1D_linear_diffusion_mod.pf
                 integration_1D_nonlinear_diffusion_mod.pf
                 )


add_test_target(equations_diffusion reg_parallel ${equations_diffusion_sources})


# Copy files needed to test equation
file(GLOB TEST_FILES "grids/*.h5")
file(COPY ${TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)



