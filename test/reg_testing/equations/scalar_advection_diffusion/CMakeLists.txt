add_test_sources(equations_advection_diffusion
                 integration_1D_linear_adv_diff_mod.pf
                 #integration_2D_linear_advection_diffusion_mod.pf
                 )


add_test_target(equations_advection_diffusion reg_parallel ${equations_advection_diffusion_sources})


# Copy files needed to test equation
#file(GLOB TEST_FILES "1D_grids/*.h5")
#file(GLOB TEST_FILES "2D_grids/*.h5")
#file(COPY ${TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)



