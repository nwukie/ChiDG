add_subdirectory(periodic)

add_lib_sources(
                type_bc.f90 
                type_bcset.f90 
                type_bcwrapper.f90 
                #mock_bc.f90 
                mod_bc.f90 
                type_bcproperty.f90
                type_bcproperty_set.f90
                type_bcdata.f90
                bc_empty.f90
                )

