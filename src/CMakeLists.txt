add_subdirectory(chidg)
add_subdirectory(constants)
add_subdirectory(grid)
add_subdirectory(containers)
add_subdirectory(tools)
add_subdirectory(io)
add_subdirectory(polynomial)
add_subdirectory(quadrature)    
add_subdirectory(solvers)
add_subdirectory(equations)
add_subdirectory(spatial)
add_subdirectory(functions)
add_subdirectory(boundary_conditions)
add_subdirectory(matrixsolvers)

# Propagate lib_srcs list to parent directory
add_lib_sources() 


# Add main executable and propagate to parent
add_exe_sources(driver.f90)    
    

