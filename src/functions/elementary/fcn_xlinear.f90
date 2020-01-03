module fcn_xlinear
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO
    use type_function,  only: function_t
    use type_point_ad,  only: point_ad_t
    use DNAD_D
    implicit none
    private


    !> X-linear function.
    !!
    !!  \f$     f(t,\vec{x}) = const + \vec{x}*slope    \f$
    !!
    !!  @author Matteo Ugolotti 
    !!  @date   10/28/2018
    !!
    !-------------------------------------------------------------------
    type, extends(function_t), public :: xlinear_f


    contains

        procedure   :: init
        procedure   :: compute

    end type xlinear_f
    !********************************************************************



contains



    !>  Add function name and options
    !!
    !!  @author Matteo Ugolotti 
    !!  @date   10/28/2018
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(xlinear_f),  intent(inout)   :: self

        !
        ! Set function name
        !
        call self%set_name("x_linear")


        !
        ! Set function options to default settings
        !
        call self%add_option('val',1._rk)
        call self%add_option('slope',1._rk)


    end subroutine init
    !*************************************************************************






    !> Function method to return function value.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/31/2018
    !!
    !--------------------------------------------------------------------
    impure elemental function compute(self,time,coord) result(val)
        class(xlinear_f),  intent(inout)   :: self
        real(rk),           intent(in)      :: time
        type(point_ad_t),   intent(in)      :: coord

        type(AD_D)    :: val,a,b,x

        ! Get x coordiantes
        x = coord%c1_

        ! Define consant term
        a = x
        a = self%get_option_value('val')

        ! Define slope
        b = x
        b = self%get_option_value('slope')
        

        ! f(x,y,z) = const
        ! derivatives are automatically set to ZERO in this case
        val = a + b*x

    end function compute
    !********************************************************************






end module fcn_xlinear
