module mod_nonlinear_solver
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_dict,              only: dict_t




    ! Import solverdata types
    use newton,                         only: newton_t
    use quasi_newton,                   only: quasi_newton_t
    implicit none



    ! Instantiate solver types for sourcing
    type(newton_t)                      :: NEWTON
    type(quasi_newton_t)                :: QUASI_NEWTON



    logical :: initialized = .false.



contains




    !> Create a concrete timescheme
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/15/2016
    !!
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------------------------------
    subroutine create_nonlinear_solver(string,instance,options)
        character(*),                           intent(in)      :: string
        class(nonlinear_solver_t), allocatable, intent(inout)   :: instance
        type(dict_t), optional,                 intent(inout)   :: options

        select case (trim(string))

            case ('newton','Newton')
                allocate(instance, source=NEWTON)

            case ('quasi_newton','Quasi_Newton')
                allocate(instance, source=QUASI_NEWTON)

            case default
                call chidg_signal(FATAL,'create_nonlinear_solver -- solver string not recognized')

        end select





        !
        ! Call options initialization if present
        !
        if (present(options)) then
            call instance%set(options)
        end if



        !
        ! Make sure the solver was allocated
        !
        if (.not. allocated(instance)) call chidg_signal(FATAL,"create_nonlinear_solver: solver was not allocated. Check that the desired solver was registered and instantiated in the mod_nonlinear_solver module")


    end subroutine create_nonlinear_solver
    !***********************************************************************************************************************








end module mod_nonlinear_solver
