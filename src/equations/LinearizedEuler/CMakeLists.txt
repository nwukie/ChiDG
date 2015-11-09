add_subdirectory(bc)

add_lib_sources(eqn_linearized_euler.f90
                LINEULER_volume_advective_flux_real.f90
                LINEULER_volume_advective_flux_imag.f90
                LINEULER_boundary_average_advective_flux_real.f90
                LINEULER_boundary_average_advective_flux_imag.f90
                LINEULER_LaxFriedrichs_flux_real.f90
                LINEULER_LaxFriedrichs_flux_imag.f90
                LINEULER_volume_advective_source_real.f90
                LINEULER_volume_advective_source_imag.f90
                LINEULER_properties.f90
                mod_linearized_euler.f90
                LINEULER_volume_advective_sourceterms_real.f90
                )
