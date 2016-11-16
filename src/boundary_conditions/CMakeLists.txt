add_subdirectory(periodic)

add_lib_sources(
                type_bc.f90
                type_bc_state.f90
                type_bc_state_wrapper.f90
                type_bc_patch.f90
                type_bcset.f90 
                mod_bc.f90 
                type_bcproperty.f90
                type_bcproperty_set.f90
                type_bc_patch_data.f90
                type_bc_group.f90
                bc_empty.f90
                type_bcset_coupling.f90
                type_bc_coupling.f90
                )

