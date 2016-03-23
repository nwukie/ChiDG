add_subdirectory(bc)

add_lib_sources(eqn_primitive_linearized_euler.f90
                PRIMLINEULER_volume_advective_flux_real.f90
                PRIMLINEULER_volume_advective_flux_imag.f90
                PRIMLINEULER_boundary_average_advective_flux_real.f90
                PRIMLINEULER_boundary_average_advective_flux_imag.f90
                PRIMLINEULER_LaxFriedrichs_flux_real.f90
                PRIMLINEULER_LaxFriedrichs_flux_imag.f90
                PRIMLINEULER_volume_advective_source_real.f90
                PRIMLINEULER_volume_advective_source_imag.f90
                PRIMLINEULER_properties.f90
                mod_primitive_linearized_euler.f90
                )
