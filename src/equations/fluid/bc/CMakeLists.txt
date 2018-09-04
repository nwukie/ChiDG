add_subdirectory(nonlocal)

add_lib_sources(
                bc_state_fluid_averaging.f90
                bc_state_inlet_total.f90
                bc_state_inlet_total_profile.f90
                bc_state_inlet_characteristic.f90
                bc_state_outlet_constant_pressure.f90
                bc_state_outlet_pressure_profile.f90
                bc_state_outlet_linear_pressure.f90
                bc_state_outlet_auxiliary_equations.f90
                bc_state_outlet_neumann_pressure_fd.f90
                bc_state_outlet_neumann_pressure_localdg.f90
                bc_state_outlet_neumann_pressure_localdg_new.f90
                bc_state_outlet_neumann_pressure_localdg_new_point.f90
                bc_state_outlet_neumann_pressure_globaldg.f90
                bc_state_outlet_neumann_LODI_localdg.f90
                bc_state_wall.f90
                bc_state_moving_wall.f90
                bc_state_fluid_extrapolate.f90
                bc_state_momentum_inlet.f90
                bc_state_symmetry.f90
                bc_state_slipwall.f90
                bc_state_farfield.f90
                bc_state_farfield_perturbation.f90
                bc_state_outlet_steady_1dchar.f90
                bc_state_outlet_3dgiles.f90
                bc_state_outlet_3dgiles_innerproduct.f90
#                bc_state_inlet_giles_quasi3d_steady.f90
#                bc_state_outlet_giles_quasi3d_steady.f90
#                bc_state_outlet_giles_quasi3d_steady_extrapolate_average.f90
#                bc_state_outlet_giles_quasi3d_unsteady_HB.f90
#                bc_state_inlet_giles_quasi3d_unsteady_HB.f90
#                bc_state_turbo_interface_steady.f90
#                bc_giles_HB_base.f90
                bc_state_outlet_characteristic_quasi3d_unsteady_HB.f90
                )
