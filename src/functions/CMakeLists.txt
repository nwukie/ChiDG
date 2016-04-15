add_subdirectory(elementary)
add_subdirectory(gaussian)
add_subdirectory(trigonometric)
add_subdirectory(polynomial)

add_lib_sources(
                type_function.f90
                type_function_wrapper.f90
                mod_function.f90
                )
