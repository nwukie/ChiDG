add_subdirectory(bc)

add_lib_sources(
                artificial_viscosity_boundary_average_operator.f90
                artificial_viscosity_volume_operator.f90
                artificial_viscosity_bc_operator.f90
                artificial_viscosity_source.f90
                type_artificial_viscosity_resolution_sensor.f90
                type_artificial_viscosity_jump_sensor.f90
                )
