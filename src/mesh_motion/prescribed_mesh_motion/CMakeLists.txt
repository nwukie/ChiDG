add_subdirectory(static)
add_subdirectory(sinusoidal)
add_subdirectory(constant_motion_1d_x)
add_subdirectory(isotropic_dilation)
add_subdirectory(hpaf_ho_ws)
add_subdirectory(rigid_body_motion)
add_subdirectory(constant_velocity)

add_lib_sources(
                type_prescribed_mesh_motion.f90
                type_prescribed_mesh_motion_domain_data.f90
                type_prescribed_mesh_motion_group.f90
                type_prescribed_mesh_motion_group_wrapper.f90
                type_prescribed_mesh_motion_function.f90
                type_prescribed_mesh_motion_function_wrapper.f90
                mod_prescribed_mesh_motion_function.f90
                )
