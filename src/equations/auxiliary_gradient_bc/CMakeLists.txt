add_subdirectory(bc)

add_lib_sources(
                auxiliary_gradient_boundary_operator.f90
                auxiliary_gradient_volume_operator.f90
                auxiliary_gradient_bc_operator.f90
                )
