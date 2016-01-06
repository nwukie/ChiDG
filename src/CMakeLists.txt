add_subdirectory(chidg)
add_subdirectory(constants)
add_subdirectory(grid)
add_subdirectory(containers)
add_subdirectory(tools)
add_subdirectory(io)
add_subdirectory(polynomial)
add_subdirectory(quadrature)    
add_subdirectory(timescheme)
add_subdirectory(equations)
add_subdirectory(spatial)
add_subdirectory(functions)
add_subdirectory(boundary_conditions)
add_subdirectory(matrixsolvers)
add_subdirectory(preconditioners)
add_subdirectory(materials)
add_subdirectory(chimera)
add_subdirectory(mpi)
#add_subdirectory(partitioners)

# Propagate lib_srcs list to parent directory
add_lib_sources() 


# Add main executable and propagate to parent
add_exe_sources(driver.f90)    
#add_exe_sources(dummy_driver.f90)    

    

