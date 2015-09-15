module atype_boundary_flux
    use mod_kinds,          only: rk, ik
    use atype_flux,         only: flux_t
    use type_mesh,          only: mesh_t
    use atype_solverdata,   only: solverdata_t
    use type_properties,    only: properties_t



    type, extends(flux_t), abstract, public :: boundary_flux_t

    contains
        procedure(compute_interface), deferred :: compute
    end type boundary_flux_t




    abstract interface
        subroutine compute_interface(self,mesh,sdata,ielem,iface,iblk,prop)
            use mod_kinds,  only: ik
            import boundary_flux_t
            import mesh_t
            import solverdata_t
            import properties_t

            class(boundary_flux_t), intent(in)      :: self
            class(mesh_t),          intent(in)      :: mesh
            class(solverdata_t),    intent(inout)   :: sdata
            integer(ik),            intent(in)      :: ielem, iface, iblk
            class(properties_t),    intent(inout)   :: prop
        end subroutine
    end interface




contains








end module atype_boundary_flux
