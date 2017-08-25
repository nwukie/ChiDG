module fcn_constant
    use mod_kinds,      only: rk,ik
    use type_function,  only: function_t
    use type_point,     only: point_t
    implicit none
    private


    !> Constant function.
    !!
    !!  \f$     f(t,\vec{x}) = const    \f$
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-------------------------------------------------------------------
    type, extends(function_t), public :: constant_f


    contains

        procedure   :: init
        procedure   :: compute

    end type constant_f
    !********************************************************************



contains



    !>  Add function name and options
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(constant_f),  intent(inout)   :: self

        !
        ! Set function name
        !
        call self%set_name("constant")


        !
        ! Set function options to default settings
        !
        call self%add_option('val',1._rk)


    end subroutine init
    !*************************************************************************








    !> Function method to return function value.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !--------------------------------------------------------------------
    impure elemental function compute(self,time,coord) result(val)
        class(constant_f),  intent(inout)   :: self
        real(rk),           intent(in)      :: time
        type(point_t),      intent(in)      :: coord

        real(rk)    :: val

        ! f(x,y,z) = const
        val = self%get_option_value('val')

    end function compute
    !********************************************************************






end module fcn_constant
