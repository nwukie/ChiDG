module eqn_scalar_advection_diffusion_ale
#include <messenger.h>
    use type_equation_set,              only: equation_set_t
    use type_equation_builder,          only: equation_builder_t
    use DNAD_D
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------------
    type, extends(equation_builder_t), public :: scalar_advection_diffusion_ale

    contains

        procedure   :: init
        procedure   :: build

    end type scalar_advection_diffusion_ale
    !******************************************************************************************************






contains






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(scalar_advection_diffusion_ale),   intent(inout)  :: self

        call self%set_name('Scalar Advection Diffusion ALE')

    end subroutine init
    !*********************************************************************************************





    !>
    !!
    !!
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------------
    function build(self,blueprint) result(scalar_advection_diffusion_ale_eqn)
        class(scalar_advection_diffusion_ale),  intent(in)  :: self
        character(*),                       intent(in)  :: blueprint

        character(:),       allocatable     :: user_msg
        type(equation_set_t)                :: scalar_advection_diffusion_ale_eqn
        

        !
        ! Set equationset name.
        !
        call scalar_advection_diffusion_ale_eqn%set_name("Scalar Advection Diffusion ALE")


        !
        ! Add spatial operators
        !
        select case (trim(blueprint))

            case('default')
                call scalar_advection_diffusion_ale_eqn%add_operator('Scalar Advection ALE Boundary Average Operator')
                call scalar_advection_diffusion_ale_eqn%add_operator('Scalar Advection ALE Volume Operator')
                call scalar_advection_diffusion_ale_eqn%add_operator('Scalar Advection ALE LaxFriedrichs Operator')
                call scalar_advection_diffusion_ale_eqn%add_operator('Scalar Advection ALE BC Operator')
                call scalar_advection_diffusion_ale_eqn%add_operator("Scalar Diffusion ALE Boundary Average Operator")
                call scalar_advection_diffusion_ale_eqn%add_operator("Scalar Diffusion ALE Volume Operator")
                call scalar_advection_diffusion_ale_eqn%add_operator("Scalar Diffusion ALE BC Operator")


            case default
                user_msg = "build_scalar_advection_diffusion_ale: I didn't recogvize the &
                            construction parameter that was passed to build the equation &
                            set."
                call chidg_signal_one(FATAL, user_msg, blueprint)

        end select


    end function build
    !*********************************************************************************************************





end module eqn_scalar_advection_diffusion_ale
