module atype_bc
    use mod_kinds,      only: rk, ik
    use type_mesh,      only: mesh_t
    implicit none
    private


    !> Abstract base-type for boundary conditions
    !!
    !!
    !!
    !-------------------------------------------------
    type, public :: bc_t
        private

        integer(ik), allocatable :: ielems(:)    !> Indices of elements associated with boundary condition

        integer :: a
        integer :: field_name
    contains
        procedure :: init
        final :: destructor
    end type bc_t

contains

    subroutine init(self,mesh)
        class(bc_t),    intent(in)  :: self
        type(mesh_t),   intent(in)  :: mesh


    end subroutine
    
    subroutine destructor(self)
        type(bc_t), intent(in) :: self
    end subroutine

end module atype_bc
