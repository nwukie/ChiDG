add_test_sources(equations_diffusion
#                 integration_1D_linear_diffusion_mod.pf
#                 integration_2D_linear_diffusion_mod.pf
                 integration_1D_nonlinear_diffusion_mod.pf
                 )


add_test_target(equations_diffusion reg_parallel ${equations_diffusion_sources})


# Copy files needed to test equation
file(GLOB 1D_TEST_FILES "1D_grids/*.h5")
file(GLOB 2D_TEST_FILES "2D_grids/*.h5")
file(COPY ${1D_TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
file(COPY ${2D_TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)



