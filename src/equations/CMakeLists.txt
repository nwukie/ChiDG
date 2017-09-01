add_subdirectory(scalar)
add_subdirectory(fluid)
add_subdirectory(dual_linear_advection)
add_subdirectory(wall_distance)
add_subdirectory(mesh_motion)
add_subdirectory(artificial_viscosity)
add_subdirectory(HB_test_cases)


#add_subdirectory(PrimitiveLinearizedEuler)
#add_subdirectory(Acoustics)



add_lib_sources(type_field.f90 
                type_equation_set.f90 
                type_equation_set_wrapper.f90
                type_equation_builder.f90
                type_equation_builder_wrapper.f90
                mod_equations.f90 
                mod_operators.f90
                mod_models.f90
                type_operator.f90
                type_operator_wrapper.f90
                type_properties.f90
                type_equationset_function_data.f90
                type_pseudo_timestep.f90
                mock_equation_set.f90
                type_model.f90
                type_model_wrapper.f90
                mock_model.f90
                )

