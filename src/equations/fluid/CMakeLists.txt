add_subdirectory(euler)
add_subdirectory(viscous)
add_subdirectory(bc)
add_subdirectory(spalart_allmaras)

add_lib_sources(
                eqn_euler.f90
                eqn_navier_stokes.f90
                type_fluid_pseudo_timestep.f90
                type_ideal_gas.f90
                )

