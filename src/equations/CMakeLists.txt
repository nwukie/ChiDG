add_subdirectory(Scalar)
add_subdirectory(LinearAdvection)
add_subdirectory(DualLinearAdvection)
add_subdirectory(Euler)


add_lib_sources(type_equation.f90 
                atype_equationset.f90 
                mod_equations.f90 
                atype_flux.f90
                atype_volume_flux.f90
                atype_boundary_flux.f90
                type_boundary_flux_wrapper.f90
                type_volume_flux_wrapper.f90
                type_properties.f90
                )

