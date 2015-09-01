module mod_solver
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use atype_solver,       only: solver_t




    ! Import solverdata types
    use solver_forward_euler,   only: forward_euler_s
    use solver_backward_euler,  only: backward_euler_s
    implicit none



    ! Instantiate solver types for sourcing
    type(forward_euler_s)   :: FORWARD_EULER
    type(backward_euler_s)  :: BACKWARD_EULER




    logical :: initialized = .false.



contains





    subroutine create_solver(solverString,solver)
        character(*),                   intent(in)      :: solverString
        class(solver_t), allocatable,   intent(inout)   :: solver

        select case (trim(solverString))
            case ('ForwardEuler','forwardeuler','Forward_Euler','forward_euler','FE','fe')
                allocate(solver, source=FORWARD_EULER)


            case ('BackwardEuler','backwardeuler','Backward_Euler','backward_euler','BE','be')
                allocate(solver, source=BACKWARD_EULER)



            case default
                call signal(FATAL,'create_solver -- solver string not recognized')
        end select



        ! Make sure the solver was allocated
        if (.not. allocated(solver)) call signal(FATAL,"create_solver: solver was not allocated. Check that the desired solver was registered and instantiated in the mod_solver module")


    end subroutine








end module mod_solver
