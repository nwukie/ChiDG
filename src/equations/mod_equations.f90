module mod_equations
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use type_equationset,               only: equationset_t

    ! Import Equations
    use eqn_scalar,                     only: scalar_e
    use eqn_linearadvection,            only: linearadvection_e
    use eqn_duallinearadvection,        only: duallinearadvection_e
    use eqn_euler,                      only: euler_e
    use eqn_linearized_euler,           only: linearized_euler_e
    implicit none

    ! Instantiate Equations
    type(scalar_e)              :: SCALAR
    type(linearadvection_e)     :: LINEARADVECTION
    type(duallinearadvection_e) :: DUALLINEARADVECTION
    type(euler_e)               :: EULER
    type(linearized_euler_e)    :: LINEULER

    logical :: uninitialized = .true.


contains


    subroutine initialize_equations()


        if (uninitialized) then
            ! List of equations to initialize
            call SCALAR%init()
            call LINEARADVECTION%init()
            call DUALLINEARADVECTION%init()
            call EULER%init()
            call LINEULER%init()

        end if

        uninitialized = .false.
    end subroutine





    !>  EquationSet Factory
    !!      - procedure for allocating a concrete instance of an equationset_t
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in] eqnstring    Character string for the equation set name
    !!  @param[in] eqnset       Allocatable equationset_t class to be instantiated
    !-------------------------------------------------------------------------------------
    subroutine create_equationset(eqnstring,eqnset)
        character(*),                      intent(in)      :: eqnstring
        class(equationset_t), allocatable, intent(inout)   :: eqnset



        select case (trim(eqnstring))
            case ('scalar','Scalar','SCALAR')
                allocate(eqnset, source=SCALAR)

            case ('linearadvection','LinearAdvection','la','LA')
                allocate(eqnset, source=LINEARADVECTION)

            case ('duallinearadvection','DualLinearAdvection','dla','DLA')
                allocate(eqnset, source=DUALLINEARADVECTION)

            case('euler','Euler','EULER')
                allocate(eqnset, source=EULER)

            case('linearizedeuler','linearized_euler','LinearizedEuler','Linearized_Euler')
                allocate(eqnset, source=LINEULER)

            case default
                call chidg_signal(FATAL,'create_equationset -- equation string not recognized')
        end select

    end subroutine



end module mod_equations
