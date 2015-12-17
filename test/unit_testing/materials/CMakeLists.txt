add_test_sources(materials
                 test_perfectgas.pf
                 )


add_test_target(materials unit ${materials_sources})
