#add_subdirectory(adjoint_edit)
#add_subdirectory(functional_edit)

add_lib_sources( 
                chidg_edit.f90 
                mod_chidg_edit_domaininfo.f90
                mod_chidg_edit_boundaryconditions.f90
                mod_chidg_edit_matrixsolver.f90
                mod_chidg_edit_timescheme.f90
                mod_chidg_edit_meshmotion.f90
                mod_chidg_edit_printoverview.f90
                mod_chidg_edit_adjoint.f90
                )
