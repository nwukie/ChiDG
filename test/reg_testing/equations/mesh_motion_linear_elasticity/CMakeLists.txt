add_test_sources(equations_mesh_motion_linear_elasticity
                #integration_1D_test_linear_elasticity_mod.pf
                 integration_1D_mesh_motion_lin_elast_mod.pf
                 integration_2D_mesh_motion_lin_elast_mod.pf
                 #integration_2D_linear_linear_elasticity_mod.pf
                 #integration_1D_nonlinear_linear_elasticity_mod.pf
                 )


add_test_target(equations_mesh_motion_linear_elasticity reg_parallel ${equations_mesh_motion_linear_elasticity_sources})

