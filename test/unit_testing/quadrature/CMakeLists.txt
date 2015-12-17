add_test_sources(quadrature
                 test_GaussLegendre.pf
                 )

add_test_target(quadrature unit ${quadrature_sources})
