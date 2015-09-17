module mock_bc
    use mod_kinds,          only: rk, ik
    use atype_bc,           only: bc_t
    use type_mesh,          only: mesh_t
    use atype_solverdata,   only: solverdata_t
    use type_properties,    only: properties_t
    implicit none
    private

    !> Mock boundary condition for testing abstract boundary condition
    !!
    !!  @author Nathan A. Wukie
    !--------------------------------------------------------
    type, extends(bc_t), public :: bc_m

    contains
        procedure   :: compute
        final       :: destructor
    end type bc_m

contains

    !> Specialized boundary condition compute procedure.
    subroutine compute(self,mesh,sdata,ielem,iface,iblk,prop)
            class(bc_m),            intent(inout)   :: self
            type(mesh_t),           intent(in)      :: mesh
            class(solverdata_t),    intent(inout)   :: sdata
            integer(ik),            intent(in)      :: ielem
            integer(ik),            intent(in)      :: iface
            integer(ik),            intent(in)      :: iblk
            class(properties_t),    intent(inout)   :: prop





    end subroutine
    
    subroutine destructor(self)
        type(bc_m), intent(in) :: self
    end subroutine

end module mock_bc
