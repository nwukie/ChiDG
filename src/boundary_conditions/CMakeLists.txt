add_subdirectory(periodic)

add_lib_sources(
                mod_bc.f90 
                type_bc_state.f90
                type_bc_state_group.f90
                type_bc_state_wrapper.f90
                type_bc_patch.f90
                type_bc_patch_group.f90
                type_domain_patch_data.f90
                type_bc_element_coupling.f90
                type_element_coupling_data.f90
                type_bcproperty.f90
                type_bcproperty_set.f90
                bc_empty.f90
                )

