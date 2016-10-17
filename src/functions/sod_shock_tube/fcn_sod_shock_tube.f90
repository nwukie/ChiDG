module fcn_sod_shock_tube
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_point,     only: point_t
    use type_function,  only: function_t
    implicit none
    private




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------------
    type, extends(function_t), public :: sod_shock_tube_f
        private


        integer(ik) :: ivar
        real(rk)    :: gam, beta

    contains

        procedure   :: init
        procedure   :: compute

    end type sod_shock_tube_f
    !********************************************************************************



contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(sod_shock_tube_f),  intent(inout) :: self

        !
        ! Set function name
        !
        call self%set_name("sod shock tube  ::   ")


    end subroutine init
    !*************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!
    !----------------------------------------------------------------------------------
    elemental function compute(self,time,coord) result(val)
        class(sod_shock_tube_f),    intent(in)  :: self
        real(rk),                   intent(in)  :: time
        type(point_t),              intent(in)  :: coord

        real(rk)                                :: val

        real(rk)    :: x,   y,   z, &
                       u, v, w, &
                       gam, beta, rho, p

        x = coord%c1_
        y = coord%c2_
        z = coord%c3_

        gam = self%gam
        beta = self%beta


        if (x <= 400._rk) then
            rho = ONE
            u   = ZERO
            v   = ZERO
            w   = ZERO
            p   = ONE
        else if (x > 400._rk) then
            rho = 0.125_rk
            u   = ZERO
            v   = ZERO
            w   = ZERO
            p   = 0.1_rk
        end if



        select case (self%ivar)
            ! RHO
            case (1)
                val = rho

            ! RHO-U
            case (2)
                val = rho*u

            ! RHO-V
            case (3)
                val = rho*v

            ! RHO-W
            case (4)
                val = rho*w

            ! RHO-E
            case (5)
                val = p/(gam-ONE)  +  HALF*rho*(u*u + v*v + w*w)

        end select

    end function compute
    !******************************************************************************************


end module fcn_sod_shock_tube
