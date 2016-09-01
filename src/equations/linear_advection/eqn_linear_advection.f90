module eqn_linear_advection
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
    type, extends(equation_builder_t), public :: linear_advection

    contains

        procedure   :: init
        procedure   :: build

    end type linear_advection
    !******************************************************************************************************






contains

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(linear_advection),   intent(inout)  :: self

        call self%set_name('LinearAdvection')

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
    function build(self,blueprint) result(linear_advection_eqn)
        class(linear_advection),    intent(in)  :: self
        character(len=*),           intent(in)  :: blueprint

        type(equation_set_t)    :: linear_advection_eqn
        

        !
        ! Set equationset name.
        !
        call linear_advection_eqn%set_name("LinearAdvection")


        !
        ! Add spatial operators
        !
        select case (trim(blueprint))

            case('default')
                call linear_advection_eqn%add_operator('LA Volume Flux')
                call linear_advection_eqn%add_operator('LA Boundary Average Flux')
                call linear_advection_eqn%add_operator('LA LaxFriedrichs Flux')

            case default
                call chidg_signal_one(FATAL, "build_linear_advection: I didn't recognize the construction &
                                              parameter that was passed to build the equation set.", blueprint)

        end select


    end function build
    !*********************************************************************************************************





end module eqn_linear_advection
