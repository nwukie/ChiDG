module eqn_test_case_poisson_equation
#include<messenger.h>
    use mod_constants
    use type_equation_set,              only: equation_set_t
    use type_equation_builder,          only: equation_builder_t
    use type_linear_coefficient_model,  only: linear_coefficient_model_t
    use type_chidg_worker,              only: chidg_worker_t
    use type_model,                     only: model_t
    use mod_models,                     only: model_factory
    use DNAD_D
    implicit none



    !>  Compose a linear diffusion equation
    !!
    !!  @author Mayank Sharma
    !!  @date   1/2/2017
    !!
    !!  \f$ \alpha \nabla u \f$
    !!
    !----------------------------------------------------------------------------------------------------
    type, extends(equation_builder_t), public   :: test_case_poisson_equation


    contains
        procedure   :: init
        procedure   :: build


    end type test_case_poisson_equation
    !****************************************************************************************************



    !>  Model for the diffusivity
    !!
    !!  Provides definition for \f$ \alpha \f$
    !!
    !!  @author Mayank Sharma
    !!  @date   1/2/2017
    !!
    !----------------------------------------------------------------------------------------------------
    type, extends(model_t), public  :: test_case_poisson_equation_model


    contains
        procedure   :: init    => init_model
        procedure   :: compute => compute_model


    end type test_case_poisson_equation_model
    !****************************************************************************************************



contains



    !>  Initialize the model with a name and the model fields it is correspoding to
    !!
    !!  @author Mayank Sharma
    !!  @date   1/2/2017
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine init_model(self)
        class(test_case_poisson_equation_model), intent(inout)  :: self

        call self%set_name('Test Case Poisson Equation Model')

        call self%add_model_field('Scalar Diffusion Coefficient')


    end subroutine init_model
    !****************************************************************************************************



    !>  Compute the diffusivity for the Poisson equation
    !!
    !!  @author Mayank Sharma   
    !!  @date   1/2/2017
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine compute_model(self,worker)
        class(test_case_poisson_equation_model), intent(in)     :: self
        type(chidg_worker_t),                    intent(inout)  :: worker

        type(AD_D), dimension(:), allocatable   :: u, alpha

        !
        ! Interpolate solution to quadrature nodes
        !
        u = worker%get_primary_field_general('u', 'value')

        !
        ! Initialize derivatives
        !
        alpha = u

        !
        ! Set default value = 0.1
        !
        alpha = TENTH

        call worker%store_model_field('Scalar Diffusion Coefficient', 'value', alpha)


    end subroutine compute_model
    !****************************************************************************************************



    !----------------------------------------------------------------------------------------------------
    !
    !                                       Implement Builder
    !
    !----------------------------------------------------------------------------------------------------

    !>
    !!
    !!  @author Mayank Sharma
    !!  @date   2/1/2017
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine init(self)
        class(test_case_poisson_equation), intent(inout)    :: self

        call self%set_name("Test Case Poisson Equation")


    end subroutine init
    !****************************************************************************************************



    !>
    !!
    !!  @author Mayank Sharma
    !!  @date   2/1/2017
    !!
    !----------------------------------------------------------------------------------------------------
    function build(self,blueprint) result(poisson_equation_eqn)
        class(test_case_poisson_equation), intent(in)   :: self
        character(len = *),                intent(in)   :: blueprint

        type(equation_set_t)                    :: poisson_equation_eqn
        type(test_case_poisson_equation_model)  :: poisson_equation_model

        !
        ! Register the Poisson equation model
        !
        call model_factory%register(poisson_equation_model)

        !
        ! Set equationset name
        !
        call poisson_equation_eqn%set_name("Test Case Poisson Equation")

        !
        ! Add spatial operator
        !
        select case (trim(blueprint))

            case('default')
                call poisson_equation_eqn%add_operator('Scalar Diffusion Boundary Average Operator')
                call poisson_equation_eqn%add_operator('Scalar Diffusion Volume Operator')
                call poisson_equation_eqn%add_operator('Scalar Diffusion BC Operator')

                call poisson_equation_eqn%add_model("Test Case Poisson Equation Model")

            case default
                call chidg_signal_one(FATAL, "build_poisson_equation: didn't recognize the construction &
                                              parameter that was passed to build the equation set", blueprint)

        end select


    end function build
    !****************************************************************************************************




















end module eqn_test_case_poisson_equation
