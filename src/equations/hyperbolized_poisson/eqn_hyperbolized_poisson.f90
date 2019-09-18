module eqn_hyperbolized_poisson
#include <messenger.h>
    use mod_constants,         only: ZERO, ONE, TWO, THREE
    use type_equation_set,     only: equation_set_t
    use type_equation_builder, only: equation_builder_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_operator,          only: operator_t
    use type_model,             only: model_t
    use mod_operators,          only: operator_factory
    use mod_models,             only: model_factory
    use DNAD_D
    use ieee_arithmetic
    implicit none



    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   11/09/2018
    !!
    !------------------------------------------------------------------------------------------
    type, extends(equation_builder_t), public :: hyperbolized_poisson

    contains

        procedure   :: init
        procedure   :: build

    end type hyperbolized_poisson
    !******************************************************************************************




    !>  A source term in the p-Poisson equation for computing Wall Distance.
    !!
    !!  S = 1.0
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   11/09/2018
    !!
    !------------------------------------------------------------------------------------------
    type, extends(operator_t), public :: hyperbolized_poisson_source_t


    contains

        procedure   :: init    => init_source
        procedure   :: compute => compute_source

    end type hyperbolized_poisson_source_t
    !******************************************************************************************






contains

    !-----------------------------------------------------------------------------------------
    !
    !                Methods implementing the Wall Distance Source term
    !
    !-----------------------------------------------------------------------------------------


    !>  Initialize Wall Distance source term.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   11/09/2018
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init_source(self)
        class(hyperbolized_poisson_source_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Hyperbolized Poisson Source")

        ! Set operator type
        call self%set_operator_type("Volume Advective Source")

        ! Set operator equations
        call self%add_primary_field("u")
        call self%add_primary_field("p")

    end subroutine init_source
    !******************************************************************************************





    !>  Implement Wall Distance source term.
    !!
    !!  S = 1.0
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   11/09/2018
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_source(self,worker,prop)
        class(hyperbolized_poisson_source_t),      intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:)   :: source_u, source_p, source_q, source_r, p, q, r

        ! Interpolate solution to quadrature nodes to initialize derivatives
        source_u = worker%get_field('u','value','element')
        source_u = ONE

        
        p = worker%get_field('p','value','element')
        source_p = -p
        q = worker%get_field('q','value','element')
        source_q = -q
        r = worker%get_field('r','value','element')
        source_r = -r


        ! Integrate volume flux
        call worker%integrate_volume_source('u',source_u)
        call worker%integrate_volume_source('p',source_p)
        call worker%integrate_volume_source('q',source_q)
        call worker%integrate_volume_source('r',source_r)


    end subroutine compute_source
    !*******************************************************************************************








    !-----------------------------------------------------------------------------------------
    subroutine init(self)
        class(hyperbolized_poisson),   intent(inout)  :: self

        call self%set_name('Hyperbolized Poisson')

    end subroutine init
    !*****************************************************************************************





    !-----------------------------------------------------------------------------------------
    function build(self,blueprint) result(hyperbolized_poisson_eqn)
        class(hyperbolized_poisson),    intent(in)  :: self
        character(*),                   intent(in)  :: blueprint

        character(:),   allocatable         :: user_msg
        type(equation_set_t)                :: hyperbolized_poisson_eqn
        type(hyperbolized_poisson_source_t) :: hyperbolized_poisson_source
        
        ! Set equationset name.
        call hyperbolized_poisson_eqn%set_name("Hyperbolized Poisson")
        
        ! Register p-Poisson model and wall distance source term
        call operator_factory%register(hyperbolized_poisson_source)

        ! Add spatial operators
        select case (trim(blueprint))

            case('default')

                ! Add the operators for the standard scalar diffusion equation:
                call hyperbolized_poisson_eqn%add_operator("Hyperbolized Poisson Boundary Average Operator")
                call hyperbolized_poisson_eqn%add_operator("Hyperbolized Poisson Volume Operator")
                call hyperbolized_poisson_eqn%add_operator("Hyperbolized Poisson BC Operator")
                call hyperbolized_poisson_eqn%add_operator("Hyperbolized Poisson LaxFriedrichs Operator")
                call hyperbolized_poisson_eqn%add_operator("Hyperbolized Poisson Source")

            case default
                user_msg = "build_wall_distance: I didn't recognize the construction &
                            parameter that was passed to build the equation set."
                call chidg_signal_one(FATAL,user_msg, blueprint)

        end select


    end function build
    !*****************************************************************************************






end module eqn_hyperbolized_poisson
