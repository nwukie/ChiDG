module eqn_scalar_diffusion
#include <messenger.h>
    use type_equation_set,              only: equation_set_t
    use type_equation_builder,          only: equation_builder_t
    use type_linear_coefficient_model,  only: linear_coefficient_model_t
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






contains

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
        class(scalar_diffusion), intent(in)  :: self
        character(len=*),        intent(in)  :: blueprint

        type(equation_set_t)                :: scalar_diffusion_eqn
        type(linear_coefficient_model_t)    :: linear_coefficient_model
        

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

                call scalar_diffusion_eqn%prop%add_scalar(linear_coefficient_model)

            case default
                call chidg_signal_one(FATAL, "build scalar diffusion: I didn't recognize the construction parameter &
                                              that was passed to build the equation set.", blueprint)


        end select


    end function build
    !*********************************************************************************************************





end module eqn_scalar_diffusion
