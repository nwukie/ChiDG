add_test_sources(equations_mesh_motion_diffusion
                 integration_1D_mesh_motion_diffusion_mod.pf
                 )


add_test_target(equations_mesh_motion_diffusion reg_parallel ${equations_mesh_motion_diffusion_sources})

