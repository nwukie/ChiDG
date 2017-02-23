module eqn_scalar_advection
#include <messenger.h>
    use mod_constants,                  only: ZERO, ONE
    use type_equation_set,              only: equation_set_t
    use type_equation_builder,          only: equation_builder_t
    use type_linear_coefficient_model,  only: linear_coefficient_model_t
    use type_chidg_worker,              only: chidg_worker_t
    use type_model,                     only: model_t
    use mod_models,                     only: model_factory
    use DNAD_D
    implicit none



    !>  Compose a scalar advection equation.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  div( vec(c)*u ) = 0
    !!
    !-----------------------------------------------------------------------------------------------------
    type, extends(equation_builder_t), public :: scalar_advection

    contains

        procedure   :: init
        procedure   :: build

    end type scalar_advection
    !******************************************************************************************************





    !>  Default model for scalar advection velocity.
    !!
    !!  Provides definition for vec(c)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/6/2016
    !!
    !------------------------------------------------------------------------------------------
    type, extends(model_t), public :: default_advection_model

    contains

        procedure   :: init    => init_model
        procedure   :: compute => compute_model

    end type default_advection_model
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
        class(default_advection_model), intent(inout)   :: self

        call self%set_name('Default Advection Velocity Model')
        call self%set_dependency('Q-')

        call self%add_model_field('Scalar Advection Velocity-1')
        call self%add_model_field('Scalar Advection Velocity-2')
        call self%add_model_field('Scalar Advection Velocity-3')

    end subroutine init_model
    !***************************************************************************************






    !>  Compute a default advection velocity vector.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/6/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute_model(self,worker)
        class(default_advection_model),     intent(in)      :: self
        type(chidg_worker_t),               intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            u, c1, c2, c3

        !
        ! Interpolate solution to quadrature nodes
        !
        u = worker%get_primary_field_general('u', 'value')


        ! Initialize derivatives
        c1 = u
        c2 = u
        c3 = u

        ! Set default values (1,0,0)
        c1 = ZERO
        c2 = ONE
        c3 = ZERO



        call worker%store_model_field('Scalar Advection Velocity-1','value',c1)
        call worker%store_model_field('Scalar Advection Velocity-2','value',c2)
        call worker%store_model_field('Scalar Advection Velocity-3','value',c3)


    end subroutine compute_model
    !***************************************************************************************











    !----------------------------------------------------------------------------------------
    !
    !                                   Implement Builder
    !
    !----------------------------------------------------------------------------------------

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !----------------------------------------------------------------------------------------
    subroutine init(self)
        class(scalar_advection),   intent(inout)  :: self

        call self%set_name('Scalar Advection')

    end subroutine init
    !****************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    function build(self,blueprint) result(scalar_advection_eqn)
        class(scalar_advection),    intent(in)  :: self
        character(*),               intent(in)  :: blueprint

        type(equation_set_t)            :: scalar_advection_eqn
        type(default_advection_model)   :: advection_model

        !
        ! Register advection model
        !
        call model_factory%register(advection_model)
        

        !
        ! Set equationset name.
        !
        call scalar_advection_eqn%set_name("Scalar Advection")


        !
        ! Add spatial operators
        !
        select case (trim(blueprint))

            case('default')
                call scalar_advection_eqn%add_operator('Scalar Advection Boundary Average Operator')
                call scalar_advection_eqn%add_operator('Scalar Advection Volume Operator')
                call scalar_advection_eqn%add_operator('Scalar Advection LaxFriedrichs Operator')
                call scalar_advection_eqn%add_operator('Scalar Advection BC Operator')


                call scalar_advection_eqn%add_model('Default Advection Velocity Model')

            case default
                call chidg_signal_one(FATAL, "build_scalar_advection: I didn't recognize the construction &
                                              parameter that was passed to build the equation set.", blueprint)

        end select


    end function build
    !****************************************************************************************





end module eqn_scalar_advection
