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
                model_ideal_gas.f90
                model_shear_stress.f90
                model_temperature_gradient.f90
                model_sutherlands_law.f90
                model_constant_viscosity.f90
                model_stokes_hypothesis.f90
                model_reynolds_analogy.f90
                model_zero_turbulent_model_fields.f90
                model_vorticity.f90
                model_fluid_wave_speed.f90
                )

