!>  Numerical mesh motion with a linear_elasticity-based model.
!!
!----------------------------------------------------------------------------------------
module eqn_mesh_motion_linear_elasticity
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
    !!  Build using the standard scalar linear_elasticity operators:
    !!      Mesh Motion linear_elasticity Boundary Average Operator
    !!      Mesh Motion linear_elasticity Volume Operator
    !!      Mesh Motion linear_elasticity BC Operator
    !!
    !!
    !!  @author Eric Wolf
    !!  @date   5/16/2017 
    !!
    !------------------------------------------------------------------------------------------
    type, extends(equation_builder_t), public :: mesh_motion_linear_elasticity

    contains

        procedure   :: init
        procedure   :: build

    end type mesh_motion_linear_elasticity
    !******************************************************************************************



    !>  A new linear_elasticity coefficient model for the scalar equations. This facilitates
    !!  the implementation of a p-Laplace operator, which is like a normal
    !!  laplace operator, but with a nonlinear linear_elasticity coefficient:
    !!
    !!      dif(mu * grad(u))
    !!
    !!  where
    !!
    !!      mu = (dudx**2 + dudy**2 + dudz**2)**((p-2)/2)
    !!
    !!  and 'p' is the p-Laplace parameter.
    !!
    !!  @author Eric Wolf
    !!  @date   5/16/2017
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
    !                   Mesh Motion linear_elasticity Coefficient Model : p-Laplacian
    !
    !---------------------------------------------------------------------------------------


    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Eric Wolf
    !!  @date   5/16/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine init_model(self)   
        class(variable_diffusivity_model), intent(inout)   :: self

        call self%set_name('Linear Elasticity Model')
        call self%set_dependency('f(Q-)')
        !call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Mesh Motion Linear Elasticity Modulus')
        call self%add_model_field('Mesh Motion Linear Elasticity Poisson Ratio')

    end subroutine init_model
    !***************************************************************************************






    !>  Compute a p-Laplace linear_elasticity coefficient.
    !!
    !!
    !!
    !!  @author Eric Wolf
    !!  @date   5/16/2017
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute_model(self,worker)
        class(variable_diffusivity_model),     intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            u, grad1_u, grad2_u, grad3_u, sumsqr, mu

        mu = ONE


        call worker%store_model_field('Mesh Motion Linear Elasticity Modulus','value',mu)
        mu = 0.3_rk
        call worker%store_model_field('Mesh Motion Linear Elasticity Poisson Ratio','value',mu)


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
    !!  @author Eric Wolf
    !!  @date   5/16/2017
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init(self)
        class(mesh_motion_linear_elasticity),   intent(inout)  :: self

        call self%set_name('Mesh Motion Linear Elasticity')

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
    !!  The practical implementation of this is to create a new nonlinear linear_elasticity 
    !!  coefficient model that will then be used in the regular scalar linear_elasticity equation.
    !!  A source term is also implemented so that the equation is an approximation of the
    !!  distance field.
    !!
    !!
    !!  @author Eric Wolf
    !!  @date  5/16/2017 
    !!
    !-----------------------------------------------------------------------------------------
    function build(self,blueprint) result(mesh_motion_linear_elasticity_eqn)
        class(mesh_motion_linear_elasticity),   intent(in)  :: self
        character(*),           intent(in)  :: blueprint

        character(:),   allocatable     :: user_msg
        type(equation_set_t)            :: mesh_motion_linear_elasticity_eqn
        type(variable_diffusivity_model)           :: variable_diffusivity
        


        !
        ! Register p-Poisson model and wall distance source term
        !
        call model_factory%register(variable_diffusivity)


        !
        ! Set equationset name.
        !
        call mesh_motion_linear_elasticity_eqn%set_name("Mesh Motion Linear Elasticity")


        !
        ! Add spatial operators
        !
        select case (trim(blueprint))

            case('default')
                ! Add the operators for the standard scalar linear_elasticity equation:
                !
                !   div(mu*grad(u)) = 0
                !
                call mesh_motion_linear_elasticity_eqn%add_operator("Mesh Motion Linear Elasticity Boundary Average Operator")
                call mesh_motion_linear_elasticity_eqn%add_operator("Mesh Motion Linear Elasticity Volume Operator")
                call mesh_motion_linear_elasticity_eqn%add_operator("Mesh Motion Linear Elasticity BC Operator")


                ! Add a definition for mu
                call mesh_motion_linear_elasticity_eqn%add_model('constant')

            case default
                user_msg = "build_mesh_motion_linear_elasticity: I didn't recognize the construction &
                            parameter that was passed to build the equation set."
                call chidg_signal_one(FATAL,user_msg, blueprint)

        end select


    end function build
    !*****************************************************************************************



















end module eqn_mesh_motion_linear_elasticity
