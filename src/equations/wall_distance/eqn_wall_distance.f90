!>  This module implements the models, source terms, and equation set for
!!  computing an approximate wall-distance field using by solving a p-Poisson
!!  equation within the DG-Chimera framework.
!!
!!  OVERVIEW:
!!  ---------
!!  This approach extends prior work that used regular Poisson equations to compute
!!  an approximate wall distance field. The prior work using Poisson equations is published 
!!  largely by Paul Tucker. The Poisson-based approach is not as accurate as some
!!  other approaches, but is easy to implement.
!!      Tucker, et al. "Transport Equation Based Wall Distance Computations Aimed at
!!                      Flows With Time-Dependent Geometry," NASA/TM-2003-212680
!!
!!  Tucker investigates some other methods based on Hamilton-Jacobi equations as well.
!!  Hamilton-Jacobi equations are not directly discretizable in a DG framework.
!!
!!  The approach taken here follows the work by Belyaev and Fayolle.
!!      Belyaev, et al. "On Variational and PDE-based Distance Function Approximations,"
!!                       COMPUTER GRAPHICS, Vol. 34, No. 8, 2015, pp. 104-118.
!!
!!  Belyaev, et al. extended the Poisson approach by using a p-Laplacian operator
!!  to create a p-Poisson equation. This has the property, that as 'p' goes to infinity,
!!  the solution of this equation satisfies the distance field. The benefit is that
!!  the approach is easy to implement, and can be made as accurate as one wants by
!!  increasing 'p' in the governing equation.
!!
!!  
!!  EQUATIONS:
!!  ----------
!!  
!!  Scalar equation: 
!!      - working variable(u)
!!      - parameter(p)
!!  
!!  div( |grad(u)|**(p-2) * grad(u) ) = 1
!!
!!
!!  Here, |grad(u)|**(p-2) is some norm of the gradient of 'u'. For the p-Laplacian it is:
!!  
!!      |grad(u)|**(p-2) = (dudx**2  +  dudy**2  +  dudz**2)**(p-2/2)
!!
!!  This can be thought of as just a nonlinear diffusion coefficient in our original
!!  scalar diffusion equation:
!!      div( mu(u) * grad(u) ) = 1
!!
!!  To implement this, we can reuse the scalar diffusion operators that have already 
!!  been implemented:
!!      - Scalar Diffusion Boundary Average Operator
!!      - Scalar Diffusion Volume Operator
!!      - Scalar Diffusion BC Operator
!!
!!  We need to implement the nonlinear diffusion coefficient model for mu(u):
!!      - Implement p-laplacian diffusion coefficient model. Used by 
!!        scalar diffusion operators.
!!
!!  We need to implement the source term:
!!      - Implement new operator (S=1)
!!
!!
!!  In the end, we build a new equation set by composing the standard scalar
!!  diffusion operators with the p-laplacian model and a unit source term.
!!
!!
!----------------------------------------------------------------------------------------
module eqn_wall_distance
#include <messenger.h>
    use mod_constants,         only: ZERO, ONE, TWO, THREE
    use type_equation_set,     only: equation_set_t
    use type_equation_builder, only: equation_builder_t
    use type_scalar,           only: scalar_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_operator,          only: operator_t
    use type_model,             only: model_t
    use mod_operators,          only: operator_factory
    use mod_models,             only: model_factory
    use DNAD_D
    implicit none



    !>  A builder for an approximate Wall Distance equation set based on a p-Poisson equation.
    !!
    !!  Build using the standard scalar diffusion operators:
    !!      Scalar Diffusion Boundary Average Operator
    !!      Scalar Diffusion Volume Operator
    !!      Scalar Diffusion BC Operator
    !!
    !!  Use the p-Laplace model implemented in this file:
    !!      p-Laplace
    !!
    !!  Use the Source term(S = 1.0) implemented in this file:
    !!      Wall Distance Source
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/16/2016
    !!
    !------------------------------------------------------------------------------------------
    type, extends(equation_builder_t), public :: wall_distance

    contains

        procedure   :: init
        procedure   :: build

    end type wall_distance
    !******************************************************************************************



    !>  A new diffusion coefficient model for the scalar equations. This facilitates
    !!  the implementation of a p-Laplace operator, which is like a normal
    !!  laplace operator, but with a nonlinear diffusion coefficient:
    !!
    !!      dif(mu * grad(u))
    !!
    !!  where
    !!
    !!      mu = (dudx**2 + dudy**2 + dudz**2)**((p-2)/2)
    !!
    !!  and 'p' is the p-Laplace parameter.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/30/2016
    !!
    !------------------------------------------------------------------------------------------
    type, extends(model_t), public :: p_laplace_model

    contains

        procedure   :: init    => init_model
        procedure   :: compute => compute_model

    end type p_laplace_model
    !******************************************************************************************






    !>  A source term in the p-Poisson equation for computing Wall Distance.
    !!
    !!  S = 1.0
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !------------------------------------------------------------------------------------------
    type, extends(operator_t), public :: wall_distance_source_t


    contains

        procedure   :: init    => init_source
        procedure   :: compute => compute_source

    end type wall_distance_source_t
    !******************************************************************************************








    !>  Parameter 'p' in the p-Poisson equation
    !!
    !!   p=2 :: linear Poisson equation
    !!   p>2 :: nonlinear p-Poisson equation
    !!
    !!  Default: p=2
    !!
    !! Procedures:
    !!   set_p_poisson_parameter
    !!   get_p_poisson_parameter
    !!
    !------------------------------------------
    real(rk)    :: p = 2._rk



contains



    !---------------------------------------------------------------------------------------
    !
    !                   Scalar Diffusion Coefficient Model : p-Laplacian
    !
    !---------------------------------------------------------------------------------------


    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init_model(self)   
        class(p_laplace_model), intent(inout)   :: self

        call self%set_name('p-Laplace')

        call self%add_model_field('Scalar Diffusion Coefficient')

    end subroutine init_model
    !***************************************************************************************






    !>  Compute a p-Laplace diffusion coefficient.
    !!
    !!  In the scalar diffusion equation: 
    !!      dif(mu * grad(u)) = S
    !!
    !!  The p-Laplace equation is given by defining the diffusion coefficient as:
    !!      mu = (dudx**2 + dudy**2 + dudz**2)**((p-2)/2)
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute_model(self,worker)
        class(p_laplace_model),     intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            u, dudx, dudy, dudz, mag2, mu

        !
        ! Interpolate solution to quadrature nodes
        !
        u    = worker%get_primary_field_general('u', 'value')
        dudx = worker%get_primary_field_general('u', 'ddx')
        dudy = worker%get_primary_field_general('u', 'ddy')
        dudz = worker%get_primary_field_general('u', 'ddz')



        !
        ! Compute magnitude of gradient
        !
        mag2 = dudx*dudx + dudy*dudy + dudz*dudz

        
        !
        ! Compute p-Laplace diffusion coefficient
        !
        ! mu = (dudx**2 + dudy**2 + dudz**2)**((p-2)/2)
        !
        if (abs(p-2._rk) > 1.e-8_rk) then
            mu = mag2**((p-TWO)/TWO)
        else
            mu = mag2
            mu = ONE
        end if



        call worker%store_model_field('Scalar Diffusion Coefficient','value',mu)


    end subroutine compute_model
    !***************************************************************************************








    !-----------------------------------------------------------------------------------------
    !
    !                Methods implementing the Wall Distance Source term
    !
    !-----------------------------------------------------------------------------------------


    !>  Initialize Wall Distance source term.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init_source(self)
        class(wall_distance_source_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Wall Distance Source")

        ! Set operator type
        call self%set_operator_type("Volume Advective Source")

        ! Set operator equations
        call self%add_primary_field("u")

    end subroutine init_source
    !******************************************************************************************





    !>  Implement Wall Distance source term.
    !!
    !!  S = 1.0
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_source(self,worker,prop)
        class(wall_distance_source_t),      intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:)   :: source
        real(rk),   allocatable, dimension(:)   :: x, y


        !
        ! Interpolate solution to quadrature nodes to initialize derivatives
        !
        source = worker%get_primary_field_element('u','value')
        source = ONE


        !
        ! Integrate volume flux
        !
        call worker%integrate_volume('u',source)


    end subroutine compute_source
    !*******************************************************************************************








    !-----------------------------------------------------------------------------------------
    !
    !                Methods implementing the Wall Distance equation builder.
    !
    !-----------------------------------------------------------------------------------------

    !>  Wall Distance builder initialization.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init(self)
        class(wall_distance),   intent(inout)  :: self

        call self%set_name('Wall Distance : p-Poisson')

    end subroutine init
    !*****************************************************************************************





    !>  Build a p-Poisson equation.
    !!
    !!  This substitutes the standard Laplace operator in the Poisson equation 
    !!  with a p-Laplace operator, creating a p-Poisson equation.
    !!
    !!  A nice property of the p-Poisson equation is that as 'p' goes to infinity,
    !!  the solution approaches the distance field, satisfying the Eikonal equation.
    !!  This allows one to obtain approximate distance fields, selecting the fidelity
    !!  of the approximation that is required.
    !!
    !!  The practical implementation of this is to create a new nonlinear diffusion 
    !!  coefficient model that will then be used in the regular scalar diffusion equation.
    !!  A source term is also implemented so that the equation is an approximation of the
    !!  distance field.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/16/2016
    !!
    !-----------------------------------------------------------------------------------------
    function build(self,blueprint) result(wall_distance_eqn)
        class(wall_distance),   intent(in)  :: self
        character(*),           intent(in)  :: blueprint

        character(:),   allocatable     :: user_msg
        type(equation_set_t)            :: wall_distance_eqn
        type(p_laplace_model)           :: p_laplace
        type(wall_distance_source_t)    :: wall_distance_source
        


        !
        ! Register p-Poisson model and wall distance source term
        !
        call model_factory%register(p_laplace)
        call operator_factory%register(wall_distance_source)


        !
        ! Set equationset name.
        !
        call wall_distance_eqn%set_name("Wall Distance : p-Poisson")


        !
        ! Add spatial operators
        !
        select case (trim(blueprint))

            case('default')
                ! Add the operators for the standard scalar diffusion equation:
                !
                !   div(mu*grad(u)) = 0
                !
                call wall_distance_eqn%add_operator("Scalar Diffusion Boundary Average Operator")
                call wall_distance_eqn%add_operator("Scalar Diffusion Volume Operator")
                call wall_distance_eqn%add_operator("Scalar Diffusion BC Operator")


                ! Add a definition for mu
                call wall_distance_eqn%add_model('p-Laplace')

                ! Add a source term (S=1)
                call wall_distance_eqn%add_operator("Wall Distance Source")

            case default
                user_msg = "build_wall_distance: I didn't recognize the construction &
                            parameter that was passed to build the equation set."
                call chidg_signal_one(FATAL,user_msg, blueprint)

        end select


    end function build
    !*****************************************************************************************













    !>  Set the 'p' parameter in the p-Poisson equation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/16/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine set_p_poisson_parameter(p_in)
        real(rk),   intent(in)  :: p_in

        p = p_in

    end subroutine set_p_poisson_parameter
    !*****************************************************************************************



    !>  Get the 'p' parameter in the p-Poisson equation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/16/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_p_poisson_parameter() result(p_out)
        real(rk)    :: p_out

        p_out = p

    end function get_p_poisson_parameter
    !*****************************************************************************************










end module eqn_wall_distance
