!>
!!
!!  @author Nathan A. Wukie
!!  @date   2/8/2016
!!
!!
!!
!--------------------------------------------------------------------------------------------------
module mod_equations
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use type_equationset,               only: equationset_t
    use type_evector,                   only: evector_t

    !
    ! Import Equations
    !
    use eqn_scalar,                     only: scalar_e
    use eqn_linearadvection,            only: linearadvection_e
    use eqn_duallinearadvection,        only: duallinearadvection_e
    use eqn_euler,                      only: euler_e
    use eqn_linearized_euler,           only: linearized_euler_e
    implicit none



    !
    ! Vector of registered equations.
    !
    type(evector_t)             :: registered_equations
    logical                     :: initialized = .false.


contains





    !>  Register equations in a module vector. This is called from chidg%init('env').
    !!
    !!  This allows the available equations to be queried in the same way that they 
    !!  are registered for allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine register_equations()
        integer :: neqns, ieqn

        !
        ! Instantiate Equations
        !
        type(scalar_e)              :: SCALAR
        type(linearadvection_e)     :: LINEARADVECTION
        type(duallinearadvection_e) :: DUALLINEARADVECTION
        type(euler_e)               :: EULER
        type(linearized_euler_e)    :: LINEARIZED_EULER


        if ( .not. initialized ) then
            !
            ! Register in global vector
            !
            call registered_equations%push_back(SCALAR)
            call registered_equations%push_back(LINEARADVECTION)
            call registered_equations%push_back(DUALLINEARADVECTION)
            call registered_equations%push_back(EULER)
            call registered_equations%push_back(LINEARIZED_EULER)





            !
            ! Initialize each equation in set. Doesn't need editing
            !
            neqns = registered_equations%size()
            do ieqn = 1,neqns
                call registered_equations%data(ieqn)%item%init()
            end do

            !
            ! Confirm initialization
            !
            initialized = .true.

        end if

    end subroutine register_equations
    !************************************************************************************







    !>  EquationSet Factory
    !!      - procedure for allocating a concrete instance of an equationset_t
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!  @param[in] eqnstring    Character string for the equation set name
    !!  @param[in] eqnset       Allocatable equationset_t class to be instantiated
    !!
    !-------------------------------------------------------------------------------------
    subroutine create_equationset(eqnstring,eqnset)
        character(*),                      intent(in)      :: eqnstring
        class(equationset_t), allocatable, intent(inout)   :: eqnset

        integer :: ierr, eindex

        !
        ! Find equation set in 'available_equations' vector
        !
        eindex = registered_equations%index_by_name(eqnstring)



        !
        ! Check equationset was found in 'available_equations'
        !
        if (eindex == 0) call chidg_signal_one(FATAL,"create_equationset: equation string not recognized", trim(eqnstring))



        !
        ! Allocate conrete equationset_t instance
        !
        allocate(eqnset, source=registered_equations%data(eindex)%item, stat=ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"create_equationset: error allocating equationset from global vector.")



        !
        ! Check equation was allocated
        !
        if ( .not. allocated(eqnset) ) call chidg_signal(FATAL,"create_equationset: error allocating conrete equation set.")



    end subroutine create_equationset
    !*************************************************************************************








    !>  This is really a utilitity for 'chidg edit' to dynamically list the avalable 
    !!  equation sets.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine list_equations()
        integer :: neqns, ieqn
        character(len=:),   allocatable :: ename
        
        neqns = registered_equations%size()

        do ieqn = 1,neqns

            ename = registered_equations%data(ieqn)%item%get_name()


            call write_line(trim(ename))
        end do ! ieqn

    end subroutine list_equations
    !**************************************************************************************


















end module mod_equations
