add_subdirectory(static)
add_subdirectory(sinusoidal)

add_lib_sources(
                type_prescribed_mesh_motion.f90
                type_prescribed_mesh_motion_wrapper.f90
                mod_prescribed_mesh_motion.f90
                )
