module eqn_dual_linear_advection
#include <messenger.h>
    use type_equation_set,      only: equation_set_t
    use type_equation_builder,  only: equation_builder_t
    implicit none


    !> This equation set exists really just to test equationsets with 
    !! more than one equation. The idea is just to compute the linear
    !! advecdtion solution twice at the same time. The equations are 
    !! independent of each other. So, we can verify, for example,
    !! the volume flux jacobians for each equation. They should be the
    !! same as for the single LinearAdvection equation set
    !!
    !!  @author Nathan A. Wukie
    !!  
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!  @note   Reformulated into an equation builder
    !!
    !-------------------------------------------------------------------------------------------------
    type, extends(equation_builder_t), public :: dual_linear_advection


    contains

        procedure   :: init
        procedure   :: build

    end type dual_linear_advection
    !*************************************************************************************************





contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(dual_linear_advection),   intent(inout)   :: self

        call self%set_name('DualLinearAdvection')

    end subroutine init
    !*********************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!  @note   Reformulated into an equation builder
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    function build(self,blueprint) result(dual_linear_advection_eqn)
        class(dual_linear_advection),   intent(in)  :: self
        character(len=*),               intent(in)  :: blueprint

        type(equation_set_t)    :: dual_linear_advection_eqn

        !
        ! Set equationset name.
        !
        call dual_linear_advection_eqn%set_name("DualLinearAdvection")


        !
        ! Add spatial operators
        !
        select case (trim(blueprint))

            case('default')
                call dual_linear_advection_eqn%add_operator('DLA Volume Flux')
                call dual_linear_advection_eqn%add_operator('DLA Boundary Average Flux')
                call dual_linear_advection_eqn%add_operator('DLA LaxFriedrichs Flux')

            case default
                call chidg_signal_one(FATAL, "build_dual_linear_advection: I didn't recognize the construction &
                                              parameter that was passed to build the equation set.", blueprint)

        end select


    end function build
    !*************************************************************************************************








end module eqn_dual_linear_advection
