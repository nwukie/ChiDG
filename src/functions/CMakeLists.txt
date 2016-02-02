add_subdirectory(elementary)
add_subdirectory(gaussian)
#add_subdirectory(isentropic_vortex)
#add_subdirectory(sod_shock_tube)
#add_subdirectory(roe_check)

add_lib_sources(type_function.f90
                mod_function.f90)
