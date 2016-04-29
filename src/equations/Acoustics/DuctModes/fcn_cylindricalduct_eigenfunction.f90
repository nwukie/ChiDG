module fcn_cylindricalduct_eigenfunction
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: HALF, ZERO
    use type_point,         only: point_t
    use type_function,      only: function_t
    implicit none



    !>  Function for computing the eigenfunction for cylindrical duct modes.
    !!  This could be passed to a root-finding routine to find the zeros,
    !!  which are the cylindrical duct mode eigenvalues.
    !!  
    !!  \f$    f(alpha_{mn}) = J'_m (alpha_{mn})    \f$
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/14/2016
    !!
    !!
    !------------------------------------------------------------------------
    type, extends(function_t) :: cylindricalduct_eigenfunction_f

    contains
        procedure :: init
        procedure :: compute
    end type cylindricalduct_eigenfunction_f
    !************************************************************************







contains



    !>  Initialization for eigenfunction class
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/14/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(cylindricalduct_eigenfunction_f),  intent(inout)   :: self

        ! Set function name
        call self%add_name("cylindricalduct_eigenfunction")

        ! Add option
        !
        ! The option 'm' here corresponds to m the azimuthal mode order
        !
        call self%add_option('m', ZERO)

    end subroutine init
    !*************************************************************************

    


    !>  Function implementation for cylindrical duct eigenfunction
    !!
    !!  \f$ f(alpha_{mn}) = J'_m (alpha_{mn})    \f$
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/14/2016
    !!
    !--------------------------------------------------------------------------
    impure elemental function compute(self,time,coord) result(val)
        class(cylindricalduct_eigenfunction_f), intent(inout)   :: self
        real(rk),                               intent(in)      :: time
        type(point_t),                          intent(in)      :: coord

        real(rk)    :: m
        real(rk)    :: x
        real(rk)    :: val

        ! Get x
        x = coord%c1_

        ! Get m
        m = self%get_option_value('m')


        ! We need the (derivative of the bessel function Jm) = Jm'
        !
        !   Jm'(x) = -J1(x)                        for   m=0
        !   Jm'(x) = 1/2( Jm-1(x)  -  Jm+1(x) )    for m/=0
        !
        if ( int(m) == 0 ) then
            val = -bessel_j1(x)
        else
            val = HALF*( bessel_jn(int(m-1),x) - bessel_jn(int(m+1),x) )
        end if

    end function compute
    !**************************************************************************










end module fcn_cylindricalduct_eigenfunction
