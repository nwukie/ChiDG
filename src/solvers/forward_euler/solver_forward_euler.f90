module solver_forward_euler
    use mod_kinds,      only: rk,ik
    use atype_solver,   only: solver_t
    use type_mesh,      only: mesh_t
    implicit none
    private


    type, extends(solver_t), public :: forward_euler_s

    contains
        procedure   :: init
        procedure   :: solve


        final :: destructor
    end type forward_euler_s

contains


    !> Solver initialization
    subroutine  init(self,mesh)
        class(forward_euler_s),     intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh

        call self%init_base(mesh)   !> Initialize base components
    end subroutine






    !> Solve method
    subroutine solve(self,mesh)
        class(forward_euler_s),     intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh

    end subroutine






    
    subroutine destructor(self)
        type(forward_euler_s),      intent(in) :: self
    end subroutine

end module solver_forward_euler
