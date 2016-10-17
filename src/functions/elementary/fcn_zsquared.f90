module fcn_zsquared
    use mod_kinds,      only: rk,ik
    use type_function,  only: function_t
    use type_point,     only: point_t
    implicit none
    private



    !> z-squared function.
    !!
    !!  \f$     f(t,\vec{x}) = z^2      \f$
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !---------------------------------------------------------------------------------
    type, extends(function_t), public :: zsquared_f
    

    contains

        procedure   :: init
        procedure   :: compute

    end type zsquared_f
    !**********************************************************************************



contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(zsquared_f),  intent(inout)   :: self

        !
        ! Set function name
        !
        call self%set_name("z_squared")


        !
        ! Set function options to default settings
        !


    end subroutine init
    !*************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!
    !---------------------------------------------------------------------------------
    impure elemental function compute(self,time,coord) result(val)
        class(zsquared_f),  intent(inout)  :: self
        real(rk),           intent(in)  :: time
        type(point_t),      intent(in)  :: coord

        real(rk)                        :: val

        ! f(x) = z**2
        val = coord%c3_  *  coord%c3_

    end function compute
    !**********************************************************************************


end module fcn_zsquared
