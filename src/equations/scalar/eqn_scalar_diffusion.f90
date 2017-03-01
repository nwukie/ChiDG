module eqn_scalar_diffusion
#include <messenger.h>
    use type_equation_set,              only: equation_set_t
    use type_equation_builder,          only: equation_builder_t
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
        class(scalar_diffusion),    intent(in)  :: self
        character(*),               intent(in)  :: blueprint

        character(:),       allocatable     :: user_msg
        type(equation_set_t)                :: scalar_diffusion_eqn
        

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


            case default
                user_msg = "build scalar diffusion: I didn't recognize the construction &
                            parameter that was passed to build the equation set."
                call chidg_signal_one(FATAL, user_msg, blueprint)


        end select


    end function build
    !*********************************************************************************************************





end module eqn_scalar_diffusion
