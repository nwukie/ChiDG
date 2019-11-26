add_subdirectory(surface_functionals)
add_subdirectory(volume_functionals)

add_lib_sources(
                type_evaluator.f90
                type_functional_wrapper.f90
                type_functional_group.f90
                mod_update_functionals.f90
                mod_functional.f90
                mod_functional_operators.f90
                )
