add_subdirectory(euler)
add_subdirectory(viscous)
add_subdirectory(bc)
add_subdirectory(spalart_allmaras)

add_lib_sources(
                eqn_euler.f90
                eqn_navier_stokes.f90
                eqn_navier_stokes_av.f90
                eqn_laminar_navier_stokes.f90
                type_fluid_pseudo_timestep.f90
                type_ideal_gas.f90
                type_sutherlands_law.f90
                type_constant_viscosity.f90
                type_stokes_hypothesis.f90
                type_reynolds_analogy.f90
                type_zero_turbulent_model_fields.f90
                type_fluid_wave_speed.f90
                )

