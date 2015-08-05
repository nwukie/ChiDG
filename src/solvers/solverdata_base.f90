module solverdata_base
    use mod_kinds,          only: rk,ik
    use atype_solverdata,   only: solverdata_t
    use type_mesh,          only: mesh_t
    implicit none
    private


    type, extends(solverdata_t), public :: base_d

    contains
        procedure   :: init

        final :: destructor
    end type base_d

contains


    !> Solver initialization
    subroutine  init(self,mesh)
        class(base_d),   intent(inout)   :: self
        type(mesh_t),   intent(in)      :: mesh

        call self%init_base(mesh)   !> Initialize base components

    end subroutine









    subroutine destructor(self)
        type(base_d),      intent(in) :: self
    end subroutine

end module solverdata_base
