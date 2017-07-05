add_subdirectory(Euler)
add_subdirectory(scalar_advection)
add_subdirectory(scalar_advection_ale)
add_subdirectory(DualLinearAdvection)
#add_subdirectory(Acoustics)

add_test_sources(equations
                test_equationset.pf
                test_equationset__add_operator.pf
                test_model.f90
                )

add_test_target(equations unit ${equations_sources})
