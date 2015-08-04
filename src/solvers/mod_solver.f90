module mod_solver
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use atype_solver,   only: solver_t


    ! Import solvers

    implicit none


    ! Instantiate solvers


    logical :: initialized = .false.



contains


    subroutine CreateSolver(solverString,solver)
        character(*),                       intent(in)      :: solverString
        class(solver_t),    allocatable,    intent(inout)   :: solver



        select case (trim(solverString))
            case default
                call signal(3,'CreateSolver -- equation string not recognized')
        end select




    end subroutine

end module mod_solver
