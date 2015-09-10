module type_volume_flux
    use mod_kinds,          only: rk, ik
    use type_mesh,          only: mesh_t
    use atype_solverdata,   only: solverdata_t



    type, public :: volume_flux_t

    contains
        procedure :: compute
    end type volume_flux_t





contains



    subroutine compute(self,mesh,sdata,ielem,iblk)
        class(volume_flux_t),   intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iblk
    end subroutine



end module type_volume_flux
