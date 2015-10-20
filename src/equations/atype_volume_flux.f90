module atype_volume_flux
    use mod_kinds,          only: rk, ik
    use atype_flux,         only: flux_t
    use type_mesh,          only: mesh_t
    use type_solverdata,    only: solverdata_t
    use type_properties,    only: properties_t
    implicit none



    type, extends(flux_t), abstract, public :: volume_flux_t

    contains
        procedure(compute_interface), deferred :: compute
    end type volume_flux_t




    abstract interface
        subroutine compute_interface(self,mesh,sdata,prop,idom,ielem,iblk)
            use mod_kinds,  only: ik
            import volume_flux_t
            import mesh_t
            import solverdata_t
            import properties_t

            class(volume_flux_t),   intent(in)      :: self
            type(mesh_t),           intent(in)      :: mesh(:)
            type(solverdata_t),     intent(inout)   :: sdata
            class(properties_t),    intent(inout)   :: prop
            integer(ik),            intent(in)      :: idom
            integer(ik),            intent(in)      :: ielem
            integer(ik),            intent(in)      :: iblk
        end subroutine
    end interface




contains








end module atype_volume_flux
