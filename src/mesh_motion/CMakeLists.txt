add_subdirectory(prescribed_mesh_motion)
add_subdirectory(rbf_mesh_motion)
add_subdirectory(structural_models)


add_lib_sources(
    mod_mesh_motion.f90
    type_mesh_motion.f90
    type_mesh_motion_wrapper.f90
    type_mesh_motion_domain_data.f90
    type_mesh_motion_group.f90
    type_mesh_motion_group_wrapper.f90
    )
