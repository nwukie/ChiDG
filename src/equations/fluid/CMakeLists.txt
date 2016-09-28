add_subdirectory(euler)
add_subdirectory(viscous)
add_subdirectory(bc)

add_lib_sources(
                eqn_euler.f90
                eqn_navier_stokes.f90
                )

