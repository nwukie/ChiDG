add_subdirectory(euler)
add_subdirectory(viscous)
add_subdirectory(bc)
add_subdirectory(spalart_allmaras)
add_subdirectory(geometric_conservation)
#add_subdirectory(rans_lowcache)

add_lib_sources(
                eqn_euler.f90
                eqn_navier_stokes.f90
                eqn_rans.f90
                mod_fluid.f90
                mod_force.f90
                mod_entropy.f90
                type_fluid_pseudo_timestep.f90
                type_ideal_gas.f90
                type_shear_stress.f90
                type_temperature_gradient.f90
                type_sutherlands_law.f90
                type_constant_viscosity.f90
                type_stokes_hypothesis.f90
                type_reynolds_analogy.f90
                type_zero_turbulent_model_fields.f90
                type_model_vorticity.f90
                type_fluid_wave_speed.f90
                )

