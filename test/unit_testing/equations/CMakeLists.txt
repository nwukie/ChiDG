add_subdirectory(Euler)
add_subdirectory(LinearAdvection)
add_subdirectory(DualLinearAdvection)
#add_subdirectory(Acoustics)

add_test_sources(equations
                test_equation_dummy.pf
                test_equationset.pf
                test_equationset__add_equation.pf
                test_equationset__add_operator.pf
                )

add_test_target(equations unit ${equations_sources})
