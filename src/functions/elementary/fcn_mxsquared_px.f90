module fcn_mxsquared_px
    use mod_kinds,      only: rk,ik
    use type_function,  only: function_t
    use type_point,     only: point_t
    use mod_constants,   only: ONE
    implicit none
    private



    !>  x-squared function.
    !!
    !!  \f$     f(t,\vec{x}) = x^2  \f$
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-------------------------------------------------------------------------
    type, extends(function_t), public :: mxsquared_px_f

    contains

        procedure   :: init
        procedure   :: compute

    end type mxsquared_px_f
    !**************************************************************************



contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(mxsquared_px_f),  intent(inout)   :: self

        !
        ! Set function name
        !
        call self%set_name("-x_squared+x")


        !
        ! Set function options to default settings
        !


    end subroutine init
    !*************************************************************************




    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!
    !--------------------------------------------------------------------------
    impure elemental function compute(self,time,coord) result(val)
        class(mxsquared_px_f),  intent(inout)  :: self
        real(rk),           intent(in)  :: time
        type(point_t),      intent(in)  :: coord

        real(rk)                        :: val

        ! f(x) = -x**2+x
        val = -(ONE) * coord%c1_  *  coord%c1_ + coord%c1_

    end function compute
    !**************************************************************************


end module fcn_mxsquared_px
