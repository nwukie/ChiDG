module mod_equations
    use mod_kinds,                      only: rk,ik
    use atype_equationset,              only: equationset_t

    ! Import Equations
    use eqn_scalar,                     only: scalar_e
    implicit none

    ! Instantiate Equations
    type(scalar_e)  :: scalar


    logical :: uninitialized = .true.


contains


    subroutine initialize_equations()


        if (uninitialized) then
            ! List of equations to initialize
            call scalar%init()

        end if

        uninitialized = .false.
    end subroutine






    subroutine AssignEquationSet(eqnstring,eqnset)
        character(*),                      intent(in)      :: eqnstring
        class(equationset_t), allocatable, intent(inout)   :: eqnset



        select case (trim(eqnstring))
            case ('scalar')
                allocate(eqnset, source=scalar)

            case default
                stop "Error: AssignEquationSet -- equation string not recognized"
        end select

    end subroutine



end module mod_equations
