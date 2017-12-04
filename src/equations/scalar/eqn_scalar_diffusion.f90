module eqn_scalar_diffusion
#include <messenger.h>
    use mod_constants,                  only: ZERO, ONE
    use type_equation_set,              only: equation_set_t
    use type_equation_builder,          only: equation_builder_t
    use type_chidg_worker,              only: chidg_worker_t
    use type_model,                     only: model_t
    use type_operator,                  only: operator_t
    use type_properties,                only: properties_t
    use mod_models,                     only: model_factory
    use mod_operators,                  only: operator_factory
    use DNAD_D
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/19/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------------
    type, extends(equation_builder_t), public :: scalar_diffusion


    contains

        procedure   :: init
        procedure   :: build

    end type scalar_diffusion
    !******************************************************************************************************


    !>  Default model for scalar diffusion coefficient.
    !!
    !!  Provides definition for diffusion coefficient, k
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/28/2017
    !!
    !------------------------------------------------------------------------------------------
    type, extends(model_t), public :: default_diffusion_model

    contains

        procedure   :: init    => init_model
        procedure   :: compute => compute_model

    end type default_diffusion_model
    !******************************************************************************************


    !>  A custom source term for the current test case
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!  S(x) = 4*pi*pi*sin(2*pi*x)
    !!
    !---------------------------------------------------------------------------------------
    type, extends(operator_t), public :: pressure_gradient_bc_source_t


    contains

        procedure   :: init    => init_source
        procedure   :: compute => compute_source

    end type pressure_gradient_bc_source_t
    !***************************************************************************************






contains


    !---------------------------------------------------------------------------------------
    !
    !                             Implement Coefficient Model
    !
    !---------------------------------------------------------------------------------------

    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init_model(self)   
        class(default_diffusion_model), intent(inout)   :: self

        call self%set_name('Default Diffusion Coefficient Model')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Scalar Diffusion Coefficient')

    end subroutine init_model
    !***************************************************************************************




    !>  Compute a default advection velocity vector.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/6/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute_model(self,worker)
        class(default_diffusion_model),     intent(in)      :: self
        type(chidg_worker_t),               intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: u, k

        !
        ! Interpolate solution to quadrature nodes
        !
        u = worker%get_field('u', 'value')


        ! Initialize derivatives
        k = u
        k = -1.0

        call worker%store_model_field('Scalar Diffusion Coefficient','value',k)

    end subroutine compute_model
    !***************************************************************************************





    !-------------------------------------------------------------------------------
    !                           Volume Source Methods
    !-------------------------------------------------------------------------------

    !>  Initialize the new volume source operator.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/31/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init_source(self)
        class(pressure_gradient_bc_source_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Pressure Gradient Source")

        ! Set operator type
        call self%set_operator_type("Volume Diffusive Flux")

        ! Set operator equations
        call self%add_primary_field("u")

    end subroutine init_source
    !********************************************************************************



    !>  Implement the volume source definition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/31/2017
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine compute_source(self,worker,prop)
        class(pressure_gradient_bc_source_t),   intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:)   :: source
        real(rk),   allocatable, dimension(:)   :: r


        !
        ! Interpolate solution to quadrature nodes
        !
        source = worker%get_field('u','grad1','element')

        r = worker%coordinate('1','volume')

        source = (200._rk/r  +  200_rk)


        !
        ! Integrate volume flux
        !
        call worker%integrate_volume_source('u',source)


    end subroutine compute_source
    !***************************************************************************************






















    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(scalar_diffusion),   intent(inout)  :: self

        call self%set_name('Scalar Diffusion')

    end subroutine init
    !*********************************************************************************************

    



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/19/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------------
    function build(self,blueprint) result(scalar_diffusion_eqn)
        class(scalar_diffusion),    intent(in)  :: self
        character(*),               intent(in)  :: blueprint

        character(:),       allocatable     :: user_msg
        type(equation_set_t)                :: scalar_diffusion_eqn

        !
        ! Instantiate/Register advection model, bc source
        !
        type(default_diffusion_model)           :: diffusion_model
        type(pressure_gradient_bc_source_t)     :: pressure_gradient_source
        call model_factory%register(diffusion_model)
        call operator_factory%register(pressure_gradient_source)


        !
        ! Set equationset name.
        !
        call scalar_diffusion_eqn%set_name("Scalar Diffusion")



        !
        ! Add spatial operators
        !
        select case (trim(blueprint))
        
            case('default')
                call scalar_diffusion_eqn%add_operator("Scalar Diffusion Boundary Average Operator")
                call scalar_diffusion_eqn%add_operator("Scalar Diffusion Volume Operator")
                call scalar_diffusion_eqn%add_operator("Scalar Diffusion BC Operator")

                call scalar_diffusion_eqn%add_model("Default Diffusion Coefficient Model")

                call scalar_diffusion_eqn%add_operator("Pressure Gradient Source")


            case default
                user_msg = "build scalar diffusion: I didn't recognize the construction &
                            parameter that was passed to build the equation set."
                call chidg_signal_one(FATAL, user_msg, blueprint)


        end select


    end function build
    !*********************************************************************************************************





end module eqn_scalar_diffusion
