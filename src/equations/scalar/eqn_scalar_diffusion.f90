module eqn_scalar_diffusion
#include <messenger.h>
    use mod_constants,                  only: ZERO, ONE
    use type_equation_set,              only: equation_set_t
    use type_equation_builder,          only: equation_builder_t
    use type_chidg_worker,              only: chidg_worker_t
    use type_model,                     only: model_t
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
        !k = -1.0_rk
        k = 0.1_rk

        call worker%store_model_field('Scalar Diffusion Coefficient','value',k)

    end subroutine compute_model
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
        call model_factory%register(diffusion_model)


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

            case default
                user_msg = "build scalar diffusion: I didn't recognize the construction &
                            parameter that was passed to build the equation set."
                call chidg_signal_one(FATAL, user_msg, blueprint)


        end select


    end function build
    !*********************************************************************************************************





end module eqn_scalar_diffusion
