module type_boundary_flux
    use mod_kinds,          only: rk, ik
    use type_flux,          only: flux_t
    use type_mesh,          only: mesh_t
    use type_solverdata,    only: solverdata_t
    use type_properties,    only: properties_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------
    type, extends(flux_t), abstract, public :: boundary_flux_t

    contains

        procedure(compute_interface), deferred :: compute

    end type boundary_flux_t
    !**********************************************************************************




    abstract interface
        subroutine compute_interface(self,mesh,sdata,prop,face_info,function_info)
            use mod_kinds,  only: ik
            import boundary_flux_t
            import mesh_t
            import solverdata_t
            import properties_t
            import face_info_t
            import function_info_t

            class(boundary_flux_t), intent(in)      :: self
            type(mesh_t),           intent(in)      :: mesh(:)
            type(solverdata_t),     intent(inout)   :: sdata
            class(properties_t),    intent(inout)   :: prop
            type(face_info_t),      intent(in)      :: face_info
            type(function_info_t),  intent(in)      :: function_info
        end subroutine
    end interface




contains








end module type_boundary_flux
