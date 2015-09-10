module atype_boundary_flux



    type, public, abstract  :: boundary_flux_t

    contains
        procedure   :: compute
    end type boundary_flux_t



contains

    subroutine compute(self,mesh,sdata,ielem,iface,iblk)
        class(boundary_flux_t), intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iface, iblk


    end subroutine



end module atype_boundary_flux
