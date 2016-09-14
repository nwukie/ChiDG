module eqn_linear_diffusion
#include <messenger.h>
    use type_equation_set,      only: equation_set_t
    use type_equation_builder,  only: equation_builder_t
    implicit none



    !>
    !!
    !!
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------------
    type, extends(equation_builder_t), public :: linear_diffusion


    contains

        procedure   :: init
        procedure   :: build

    end type linear_diffusion
    !******************************************************************************************************






contains

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(linear_diffusion),   intent(inout)  :: self

        call self%set_name('Linear Diffusion')

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
    function build(self,blueprint) result(linear_diffusion_eqn)
        class(linear_diffusion), intent(in)  :: self
        character(len=*),        intent(in)  :: blueprint

        type(equation_set_t)    :: linear_diffusion_eqn
        

        !
        ! Set equationset name.
        !
        call linear_diffusion_eqn%set_name("Linear Diffusion")



        !
        ! Add spatial operators
        !
        select case (trim(blueprint))
        
            case('default')
                call linear_diffusion_eqn%add_operator("Linear Diffusion Boundary Average Flux")
                call linear_diffusion_eqn%add_operator("Linear Diffusion Volume Flux")
                call linear_diffusion_eqn%add_operator("Linear Diffusion Volume Source")
                call linear_diffusion_eqn%add_operator("Linear Diffusion BC Flux")

            case default
                call chidg_signal_one(FATAL, "build lineardiffusion: I didn't recognize the construction parameter &
                                              that was passed to build the equation set.", blueprint)


        end select


    end function build
    !*********************************************************************************************************





end module eqn_linear_diffusion
