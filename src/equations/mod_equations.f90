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
    use eqn_linear_advection,       only: linear_advection
    use eqn_linear_diffusion,       only: linear_diffusion
    use eqn_dual_linear_advection,  only: dual_linear_advection
    use eqn_euler,                  only: euler 
    implicit none



    !
    ! Vector of registered equations.
    !
    type(evector_t)             :: registered_equation_builders
    logical                     :: initialized = .false.



contains




    !>  Register equation builders in a module vector. This is called from chidg%init('env').
    !!
    !!  This allows the available equations to be queried in the same way that they 
    !!  are registered for allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/2/2016
    !!  @note   Modified to take equation_builder's instead of equation_set's
    !!
    !-----------------------------------------------------------------------------------
    subroutine register_equation_builders()
        integer :: ibuild

        !
        ! Instantiate Equations
        !
        type(euler)                 :: euler_builder
        type(linear_advection)      :: linear_advection_builder
        type(linear_diffusion)      :: linear_diffusion_builder
        type(dual_linear_advection) :: dual_linear_advection_builder


        !
        ! Register if needed
        !
        if ( .not. initialized ) then

            ! Register in global vector
            call registered_equation_builders%push_back(euler_builder)
            call registered_equation_builders%push_back(linear_advection_builder)
            call registered_equation_builders%push_back(linear_diffusion_builder)
            call registered_equation_builders%push_back(dual_linear_advection_builder)




            ! Initialize each builder
            do ibuild = 1,registered_equation_builders%size()
                call registered_equation_builders%data(ibuild)%bld%init()
            end do

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

        integer                                 :: ierr, bindex
        type(equation_set_t)                    :: eqnset
        class(equation_builder_t), allocatable  :: builder

        !
        ! Find equation set in 'available_equations' vector
        !
        bindex = registered_equation_builders%index_by_name(eqnstring)


        !
        ! Check equationset was found in 'available_equations'
        !
        if (bindex == 0) call chidg_signal_one(FATAL,"build_equation_set: equation string not recognized", trim(eqnstring))


        !
        ! Get equation set builder
        !
        !builder = registered_equation_builders%at(bindex)
        allocate(builder, source=registered_equation_builders%at(bindex), stat=ierr)
        if (ierr /= 0) call AllocationError


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

            ename = registered_equation_builders%data(ieqn)%bld%get_name()
            call write_line(trim(ename))

        end do ! ieqn

    end subroutine list_equations
    !**************************************************************************************


















end module mod_equations
