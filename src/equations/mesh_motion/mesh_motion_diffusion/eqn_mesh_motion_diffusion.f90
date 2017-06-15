!>  Numerical mesh motion with a diffusion-based model.
!!
!----------------------------------------------------------------------------------------
module eqn_mesh_motion_diffusion
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



    !>  A builder for an approximate Wall Distance equation set based on a p-Poisson equation.
    !!
    !!  Build using the standard scalar diffusion operators:
    !!      Mesh Motion Diffusion Boundary Average Operator
    !!      Mesh Motion Diffusion Volume Operator
    !!      Mesh Motion Diffusion BC Operator
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
    type, extends(equation_builder_t), public :: mesh_motion_diffusion

    contains

        procedure   :: init
        procedure   :: build

    end type mesh_motion_diffusion
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
    type, extends(model_t), public :: variable_diffusivity_model

    contains

        procedure   :: init    => init_model
        procedure   :: compute => compute_model

    end type variable_diffusivity_model
    !******************************************************************************************












contains



    !---------------------------------------------------------------------------------------
    !
    !                   Mesh Motion Diffusion Coefficient Model : p-Laplacian
    !
    !---------------------------------------------------------------------------------------


    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init_model(self)   
        class(variable_diffusivity_model), intent(inout)   :: self

        call self%set_name('Variable Diffusivity')
        call self%set_dependency('f(Q-)')
        !call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Mesh Motion Diffusion Coefficient')

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
        class(variable_diffusivity_model),     intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            u, grad1_u, grad2_u, grad3_u, sumsqr, mu

        mu = ONE


        call worker%store_model_field('Mesh Motion Diffusion Coefficient','value',mu)


    end subroutine compute_model
    !***************************************************************************************








    !-----------------------------------------------------------------------------------------
    !
    !                Methods implementing the Wall Distance Source term
    !
    !-----------------------------------------------------------------------------------------






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
        class(mesh_motion_diffusion),   intent(inout)  :: self

        call self%set_name('Mesh Motion Diffusion')

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
    function build(self,blueprint) result(mesh_motion_diffusion_eqn)
        class(mesh_motion_diffusion),   intent(in)  :: self
        character(*),           intent(in)  :: blueprint

        character(:),   allocatable     :: user_msg
        type(equation_set_t)            :: mesh_motion_diffusion_eqn
        type(variable_diffusivity_model)           :: variable_diffusivity
        


        !
        ! Register p-Poisson model and wall distance source term
        !
        call model_factory%register(variable_diffusivity)


        !
        ! Set equationset name.
        !
        call mesh_motion_diffusion_eqn%set_name("Mesh Motion Diffusion")


        !
        ! Add spatial operators
        !
        select case (trim(blueprint))

            case('default')
                ! Add the operators for the standard scalar diffusion equation:
                !
                !   div(mu*grad(u)) = 0
                !
                call mesh_motion_diffusion_eqn%add_operator("Mesh Motion Diffusion Boundary Average Operator")
                call mesh_motion_diffusion_eqn%add_operator("Mesh Motion Diffusion Volume Operator")
                call mesh_motion_diffusion_eqn%add_operator("Mesh Motion Diffusion BC Operator")


                ! Add a definition for mu
                call mesh_motion_diffusion_eqn%add_model('constant')

            case default
                user_msg = "build_mesh_motion_diffusion: I didn't recognize the construction &
                            parameter that was passed to build the equation set."
                call chidg_signal_one(FATAL,user_msg, blueprint)

        end select


    end function build
    !*****************************************************************************************



















end module eqn_mesh_motion_diffusion
