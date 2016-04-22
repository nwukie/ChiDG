add_subdirectory(bc)

add_lib_sources(eqn_primitive_linearized_euler_axisymmetric.f90
                PRIMLINEULERAXI_volume_advective_flux_real.f90
                PRIMLINEULERAXI_volume_advective_flux_imag.f90
                PRIMLINEULERAXI_boundary_average_advective_flux_real.f90
                PRIMLINEULERAXI_boundary_average_advective_flux_imag.f90
                PRIMLINEULERAXI_LaxFriedrichs_flux_real.f90
                PRIMLINEULERAXI_LaxFriedrichs_flux_imag.f90
                PRIMLINEULERAXI_volume_advective_source_real.f90
                PRIMLINEULERAXI_volume_advective_source_imag.f90
                PRIMLINEULERAXI_properties.f90
                PRIMLINEULERAXI_volume_advective_sourceterms_real.f90
                mod_primitive_linearized_euler_axisymmetric.f90
                )
