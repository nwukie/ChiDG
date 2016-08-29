module type_volume_flux
    use mod_kinds,          only: rk, ik
    use type_flux,          only: flux_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
!    use type_mesh,          only: mesh_t
!    use type_solverdata,    only: solverdata_t
!    use type_element_info,  only: element_info_t
!    use type_function_info, only: function_info_t
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(flux_t), abstract, public :: volume_flux_t

    contains

        procedure(compute_interface), deferred :: compute

    end type volume_flux_t
    !*******************************************************************************




    abstract interface
        !subroutine compute_interface(self,mesh,sdata,prop,elem_info,function_info)
        subroutine compute_interface(self,worker,prop)
            use mod_kinds,  only: ik
            import volume_flux_t
            import chidg_worker_t
            import properties_t
!            import mesh_t
!            import solverdata_t
!            import element_info_t
!            import function_info_t

            class(volume_flux_t),   intent(in)      :: self
            type(chidg_worker_t),   intent(inout)   :: worker
            class(properties_t),    intent(inout)   :: prop
!            type(mesh_t),           intent(in)      :: mesh(:)
!            type(solverdata_t),     intent(inout)   :: sdata
!            type(element_info_t),   intent(in)      :: elem_info
!            type(function_info_t),  intent(in)      :: function_info
        end subroutine
    end interface




contains








end module type_volume_flux
