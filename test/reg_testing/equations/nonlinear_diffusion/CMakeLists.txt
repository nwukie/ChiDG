add_test_sources(equations_diffusion
                 integration_1D_linear_diffusion_mod.pf
                 )


add_test_target(equations_linear_diffusion reg ${equations_linear_diffusion_sources})


# Copy files needed to test equation
file(GLOB TEST_FILES "*.h5")
file(COPY ${TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)



