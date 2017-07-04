add_test_sources(quadrature
                 test_GaussLegendre.pf
                 test_reference_element.pf
                 )

add_test_target(quadrature unit ${quadrature_sources})
