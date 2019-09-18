add_subdirectory(bc)

add_lib_sources(eqn_primitive_linearized_euler.f90
                mod_primitive_linearized_euler.f90
                PRIMLINEULER_volume_advection.f90
                PRIMLINEULER_boundary_average.f90
                PRIMLINEULER_bc.f90
                PRIMLINEULER_LaxFriedrichs.f90
                PRIMLINEULER_temporal_source.f90
                PRIMLINEULER_axial_source.f90
                PRIMLINEULER_circumferential_source.f90
                PRIMLINEULER_equation_source.f90
                PRIMLINEULER_divergence_source.f90
                )
