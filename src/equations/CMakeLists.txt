add_subdirectory(linear_advection)
#add_subdirectory(LinearDiffusion)
add_subdirectory(dual_linear_advection)
add_subdirectory(euler)
#add_subdirectory(PrimitiveLinearizedEuler)



#add_subdirectory(Acoustics)



add_lib_sources(type_equation.f90 
                type_equation_set.f90 
                type_equation_builder.f90
                type_equation_builder_wrapper.f90
                mod_equations.f90 
#                type_flux.f90
                mod_operators.f90
                type_operator.f90
                type_operator_wrapper.f90
                type_bc_operator.f90
                type_bc_operator_wrapper.f90
                type_bc_patch.f90
                type_bc.f90
                type_properties.f90
                mock_equation_set.f90
                type_equationset_function_data.f90
                )

