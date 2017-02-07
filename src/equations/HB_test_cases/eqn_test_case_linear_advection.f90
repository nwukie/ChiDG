module eqn_test_case_linear_advection
#include<messenger.h>
    use mod_constants,                  only: ZERO,ONE
    use type_equation_set,              only: equation_set_t
    use type_equation_builder,          only: equation_builder_t
    use type_linear_coefficient_model,  only: linear_coefficient_model_t
    use type_chidg_worker,              only: chidg_worker_t
    use type_model,                     only: model_t   
    use mod_models,                     only: model_factory
    use DNAD_D
    implicit none



    !>  Compose a linear advection equation
    !!
    !!  @author Mayank Sharma
    !!  @date   1/2/2017
    !!
    !!`\f$ \nabla \cdot (\vec{c} u) \f$
    !!
    !---------------------------------------------------------------------------------------------
    type, extends(equation_builder_t), public   :: test_case_linear_advection


    contains
        procedure   :: init
        procedure   :: build


    end type test_case_linear_advection
    !*********************************************************************************************



    !>  Model for the linear advection velocity
    !!
    !!  Provides definition for \vec{c}
    !!
    !!  @author Mayank Sharma
    !!  @date   31/1/2017
    !!
    !---------------------------------------------------------------------------------------------
    type, extends(model_t), public  :: test_case_linear_advection_model


    contains
        procedure   :: init    => init_model
        procedure   :: compute => compute_model


    end type test_case_linear_advection_model
    !*********************************************************************************************



contains



    !>  Initialize model with a name and the model fields it is contributing to
    !!
    !!  @author Mayank Sharma
    !!  @date   1/2/2017
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init_model(self)
        class(test_case_linear_advection_model), intent(inout)  :: self

        call self%set_name('Test Case Linear Advection Model')

        call self%add_model_field('Linear X-Advection Velocity')
        call self%add_model_field('Linear Y-Advection Velocity')
        call self%add_model_field('Linear Z-Advection Velocity')


    end subroutine init_model
    !*********************************************************************************************



    !>  Compute a linear advection velocity vector
    !!
    !!  @author Mayank Sharma
    !!  @sate   1/2/2017
    !!
    !---------------------------------------------------------------------------------------------
    subroutine compute_model(self,worker)
        class(test_case_linear_advection_model), intent(in)     :: self
        type(chidg_worker_t),                    intent(inout)  :: worker

        type(AD_D), dimension(:), allocatable   :: u, cx, cy, cz

        !
        ! Interpolate solution to quadrature nodes
        !
        u = worker%get_primary_field_general('u', 'value')

        !
        ! Initialize derivatives
        !
        cx = u
        cy = u
        cz = u

        !
        ! Set default values (1,0,0)
        !
        cx = ONE
        cy = ZERO
        cz = ZERO

        call worker%store_model_field('Linear X-Advection Velocity', 'value', cx)
        call worker%store_model_field('Linear Y-Advection Velocity', 'value', cy)
        call worker%store_model_field('Linear Z-Advection Velocity', 'value', cz)

    
    end subroutine compute_model
    !*********************************************************************************************



    !---------------------------------------------------------------------------------------------
    !
    !                                   Implement Builder
    !
    !---------------------------------------------------------------------------------------------

    !>
    !!
    !!  @author Mayank Sharma
    !!  @date   1/2/2017
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(test_case_linear_advection), intent(inout)    :: self

        call self%set_name('Test Case Linear Advection')


    end subroutine init
    !*********************************************************************************************



    !>
    !!
    !!  @author Mayank Sharma
    !!  @date   1/2/2017
    !!
    !---------------------------------------------------------------------------------------------
    function build(self,blueprint) result(linear_advection_eqn)
        class(test_case_linear_advection), intent(in)   :: self
        character(len = *),                intent(in)   :: blueprint

        type(equation_set_t)                    :: linear_advection_eqn
        type(test_case_linear_advection_model)  :: linear_advection_model

        !
        ! Register the linear advection model
        !
        call model_factory%register(linear_advection_model)

        !
        ! Set equationset name
        !
        call linear_advection_eqn%set_name("Test Case Linear Advection")

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))

            case('default')
                call linear_advection_eqn%add_operator('Scalar Advection Boundary Average Operator')
                call linear_advection_eqn%add_operator('Scalar Advection Volume Operator')
                call linear_advection_eqn%add_operator('Scalar Advection LaxFriedrichs Operator')
                call linear_advection_eqn%add_operator('Scalar Advection BC Operator')

                call linear_advection_eqn%add_model('Test Case Linear Advection Model')

            case default
                call chidg_signal_one(FATAL, "build_linear_adveciton: didn't recognize the construction &
                                              parameter that was passed to build the equation set", blueprint)

        end select


    end function build
    !*********************************************************************************************




















end module eqn_test_case_linear_advection
