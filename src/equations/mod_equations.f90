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
    use mod_kinds,              only: rk,ik
    use mod_string,             only: string_t
    use type_equation_set,      only: equation_set_t
    use type_equation_builder,  only: equation_builder_t
    use type_evector,           only: evector_t

    !
    ! Import Equations
    !
!    use eqn_linearadvection,                         only: linearadvection_e
!    use eqn_lineardiffusion,                         only: lineardiffusion_e
!    use eqn_duallinearadvection,                     only: duallinearadvection_e
!    use eqn_euler,                                   only: euler_e
!    use eqn_primitive_linearized_euler,              only: primitive_linearized_euler_e

    use eqn_euler,  only: euler 
    implicit none



    !
    ! Vector of registered equations.
    !
    type(evector_t)             :: registered_equation_builders
    logical                     :: initialized = .false.




!    !>
!    !!
!    !!
!    !!
!    !!
!    !---------------------------------------------------------------------------------------------
!    type, public :: equation_builder_t
!
!        type(string_t)  :: name
!
!    contains
!
!        procedure(init_interface),  deferred    :: init
!        procedure(build_interface), deferred    :: build
!
!        procedure                               :: set_name
!        procedure                               :: get_name
!
!    end type equation_builder_t
!    !*********************************************************************************************
!    abstract interface
!        subroutine init_interface(self)
!            import equation_builder_t
!
!            class(equation_builder_t),  intent(inout)   :: self
!        end subroutine
!    end interface
!
!    abstract interface
!        function build_interface(self,blueprint)
!            import equation_builder_t
!
!            class(equation_builder_t),  intent(in)  :: self
!            character(len=*),           intent(in)  :: blueprint
!        end function
!    end interface
!    !*********************************************************************************************

contains


!    !>
!    !!
!    !!
!    !!
!    !!
!    !---------------------------------------------------------------------------------------------
!    subroutine set_name(self,string)
!        class(equation_builder_t),  intent(inout)   :: self
!        character(len=*),           intent(in)      :: string
!
!        self%name = string
!
!    end subroutine set_name
!    !*********************************************************************************************
!    
!
!    
!
!    !>
!    !!
!    !!
!    !!
!    !!
!    !!
!    !---------------------------------------------------------------------------------------------
!    function get_name(self) result(builder_name)
!        class(equation_builder_t),  intent(in)  :: self
!
!        character(len=:), allocatable   :: builder_name
!
!        builder_name = self%name%str
!
!    end function get_name
!    !*********************************************************************************************













    !>  Register equation builders in a module vector. This is called from chidg%init('env').
    !!
    !!  This allows the available equations to be queried in the same way that they 
    !!  are registered for allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine register_equation_builders()
        integer :: neqns, ieqn

        !
        ! Instantiate Equations
        !
        type(euler) :: euler_builder


        !
        ! Register if needed
        !
        if ( .not. initialized ) then

            ! Register in global vector
            call registered_equation_builders%push_back(euler_builder)

            ! Confirm initialization
            initialized = .true.

        end if

    end subroutine register_equation_builders
    !*************************************************************************************







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
    function build_equation_set(eqnstring,blueprint) result(eqnset)
        character(len=*),   intent(in)      :: eqnstring
        character(len=*),   intent(in)      :: blueprint

        integer                 :: ierr, bindex
        type(equation_set_t)    :: eqnset

        !
        ! Find equation set in 'available_equations' vector
        !
        bindex = registered_equations%index_by_name(eqnstring)


        !
        ! Check equationset was found in 'available_equations'
        !
        if (bindex == 0) call chidg_signal_one(FATAL,"create_equationset: equation string not recognized", trim(eqnstring))


        !
        ! Get equation set builder
        !
        builder = registered_equation_builders%at(bindex)


        !
        ! Build equation set
        !
        eqnset = builder%build(blueprint)


    end function build_equation_set
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
        
        neqns = registered_equation_builders%size()

        do ieqn = 1,neqns

            ename = registered_equation_builders%data(ieqn)%get_name()
            call write_line(trim(ename))

        end do ! ieqn

    end subroutine list_equations
    !**************************************************************************************


















end module mod_equations
