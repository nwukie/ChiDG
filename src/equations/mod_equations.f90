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
    use eqn_scalar_advection,               only: scalar_advection
    use eqn_scalar_diffusion,               only: scalar_diffusion
    use eqn_dual_linear_advection,          only: dual_linear_advection
    use eqn_euler,                          only: euler 
    use eqn_navier_stokes,                  only: navier_stokes
    use eqn_rans,                           only: rans
    use eqn_wall_distance,                  only: wall_distance
    use eqn_mesh_motion_diffusion,          only: mesh_motion_diffusion 
    use eqn_mesh_motion_linear_elasticity,  only: mesh_motion_linear_elasticity
    use eqn_test_case_linear_advection,     only: test_case_linear_advection
    use eqn_test_case_poisson_equation,     only: test_case_poisson_equation
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/19/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    type, public :: equation_set_factory_t

        type(evector_t) :: equation_sets

    contains
        
        procedure   :: register
        procedure   :: produce  => produce_equation_set
        procedure   :: list     => list_equation_sets
        procedure   :: has      => has_equation_set

    end type equation_set_factory_t
    !**************************************************************************************




    type(equation_set_factory_t),   target  :: equation_set_factory
    logical                                 :: initialized = .false.



contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/19/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine register(self,eqn_set)
        class(equation_set_factory_t),  intent(inout)   :: self
        type(equation_set_t),           intent(in)      :: eqn_set

        call self%equation_sets%push_back(eqn_set)

    end subroutine register
    !***************************************************************************************










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
        type(scalar_advection)              :: scalar_advection_builder
        type(scalar_diffusion)              :: scalar_diffusion_builder
        type(dual_linear_advection)         :: dual_linear_advection_builder
        type(euler)                         :: euler_builder
        type(navier_stokes)                 :: navier_stokes_builder
        type(rans)                          :: rans_builder
        type(wall_distance)                 :: wall_distance_builder
        type(mesh_motion_diffusion)         :: mesh_motion_diffusion_builder
        type(mesh_motion_linear_elasticity) :: mesh_motion_linear_elasticity_builder
        type(test_case_linear_advection)    :: test_case_linear_advection_builder
        type(test_case_poisson_equation)    :: test_case_poisson_equation_builder


        !
        ! Register if needed
        !
        if ( .not. initialized ) then

            ! Register in global vector
            call equation_set_factory%register(scalar_advection_builder%build('default'))
            call equation_set_factory%register(scalar_diffusion_builder%build('default'))
            call equation_set_factory%register(dual_linear_advection_builder%build('default'))
            call equation_set_factory%register(euler_builder%build('default'))
            call equation_set_factory%register(navier_stokes_builder%build('default'))
            call equation_set_factory%register(rans_builder%build('default'))
            call equation_set_factory%register(wall_distance_builder%build('default'))
            call equation_set_factory%register(mesh_motion_diffusion_builder%build('default'))
            call equation_set_factory%register(mesh_motion_linear_elasticity_builder%build('default'))
            call equation_set_factory%register(test_case_linear_advection_builder%build('default'))
            call equation_set_factory%register(test_case_poisson_equation_builder%build('default'))



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
    function produce_equation_set(self,eqnstring,blueprint) result(eqnset)
        class(equation_set_factory_t),  intent(inout)   :: self
        character(*),                   intent(in)      :: eqnstring
        character(*),                   intent(in)      :: blueprint

        character(:),               allocatable :: user_msg, dev_msg
        integer                                 :: ierr, bindex
        type(equation_set_t)                    :: eqnset

        !
        ! Find equation set in 'available_equations' vector
        !
        bindex = self%equation_sets%index_by_name(eqnstring)


        !
        ! Check equationset was found in 'available_equations'
        !
        user_msg = "We can't seem to find an equation set that matches the string that was passed &
                    into the equation set factory. &
                    Maybe check that the equation set strings that were set for the domains &
                    are all valid."
        dev_msg = "Check that the equation set builder is registered properly in &
                   the equation set builder factory: src/equations/mod_equations.f90. When an equation &
                   set builder is defined, is needs to be registered in the factory by calling &
                   'call equation_set_factory%register(builder)', where 'builder' is the &
                   object that knows how to build the equation set. This could be done in  &
                   'mod_equations.register_equation_builders'. The 'register_equation_builders' routine &
                   gets called on startup and loads all the default equation builders into the factory &
                   so the library knows what equation sets it can build."
        if (bindex == 0) call chidg_signal_two(OOPS,user_msg,trim(eqnstring),dev_msg=dev_msg)


        !
        ! Get equation set builder
        !
        eqnset = self%equation_sets%at(bindex)



    end function produce_equation_set
    !*************************************************************************************








    !>  This is really a utilitity for 'chidg edit' to dynamically list the avalable 
    !!  equation sets.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine list_equation_sets(self)
        class(equation_set_factory_t),  intent(in)   :: self

        integer :: neqns, ieqn
        character(:),   allocatable :: ename
        
        neqns = self%equation_sets%size()

        do ieqn = 1,neqns

            ename = self%equation_sets%data(ieqn)%get_name()
            call write_line(trim(ename))

        end do ! ieqn

    end subroutine list_equation_sets
    !**************************************************************************************





    !>  Check if a given equation builder is registered in the factory.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2017
    !!
    !--------------------------------------------------------------------------------------
    function has_equation_set(self,equation_string) result(equation_status)
        class(equation_set_factory_t),  intent(in)  :: self
        character(*),                   intent(in)  :: equation_string



        integer                     :: neqns, ieqn
        character(:),   allocatable :: ename
        logical                     :: equation_status
        
        equation_status = .false.
        neqns = self%equation_sets%size()
        do ieqn = 1,neqns

            ename = self%equation_sets%data(ieqn)%get_name()

            equation_status = (trim(equation_string) == trim(ename))
            if (equation_status) exit

        end do ! ieqn


    end function has_equation_set
    !**************************************************************************************









end module mod_equations
