module eqn_scalar_advection
#include <messenger.h>
    use type_equation_set,      only: equation_set_t
    use type_equation_builder,  only: equation_builder_t
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------------
    type, extends(equation_builder_t), public :: scalar_advection

    contains

        procedure   :: init
        procedure   :: build

    end type scalar_advection
    !******************************************************************************************************






contains

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(scalar_advection),   intent(inout)  :: self

        call self%set_name('Scalar Advection')

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
    function build(self,blueprint) result(scalar_advection_eqn)
        class(scalar_advection),    intent(in)  :: self
        character(len=*),           intent(in)  :: blueprint

        type(equation_set_t)    :: scalar_advection_eqn
        

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

            case default
                call chidg_signal_one(FATAL, "build_scalar_advection: I didn't recognize the construction &
                                              parameter that was passed to build the equation set.", blueprint)

        end select


    end function build
    !*********************************************************************************************************





end module eqn_scalar_advection
