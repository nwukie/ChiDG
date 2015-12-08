module fcn_isentropic_vortex
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_point,     only: point_t
    use atype_function, only: function_t
    implicit none
    private

    type, extends(function_t), public :: isentropic_vortex_f
        private


        ! constants in the gaussian function
        !
        ! f(x) = a exp(- (x-b)**2 / 2c**2)
        !
        real(rk)    :: rho = ONE
        real(rk)    :: p   = ONE


        real(rk)    :: gam  = 1.4_rk
        real(rk)    :: beta = FIVE
        real(rk)    :: xo
        real(rk)    :: yo
        real(rk)    :: zo

        real(rk)    :: uinf = ONE
        real(rk)    :: vinf = ONE
        real(rk)    :: winf = ZERO

        integer(ik) :: ivar

    contains
        procedure   :: order
        procedure   :: calc
        procedure   :: set
    end type isentropic_vortex_f



contains

    subroutine set(self,valstring,val)
        class(isentropic_vortex_f), intent(inout)   :: self
        character(*),               intent(in)      :: valstring
        real(rk),                   intent(in)      :: val


        select case (valstring)

            ! Function settings
            case('uinf','Uinf')
                self%uinf = val
            case('vinf','Vinf')
                self%vinf = val
            case('winf','Winf')
                self%winf = val
            case('beta','Beta')
                self%beta = val

            
            case('xo','XO')
                self%xo = val
            case('yo','YO')
                self%yo = val
            case('zo','ZO')
                self%zo = val




            ! Variables available
            case('var','Var','variable','Variable')
                self%ivar = NINT(val)


            case default
                call chidg_signal(FATAL,'gaussian_f%set: Invalid option string')
        end select


    end subroutine


    function order(self)
        class(isentropic_vortex_f), intent(in)   :: self
        integer(ik)                     :: order

        order = 3

    end function



    elemental function calc(self,pt)
        class(isentropic_vortex_f),  intent(in)  :: self
        type(point_t),               intent(in)  :: pt
        real(rk)                                 :: calc

        real(rk)    :: x,   y,   z, &
                       du, dv, u, v, w, &
                       gam, beta, r, T, rho, p

        x = pt%c1_
        y = pt%c2_
        z = pt%c3_

        gam = self%gam
        beta = self%beta

        r = sqrt((x - self%xo)**TWO + (y - self%yo)**TWO)
        T = ONE - ((gam - ONE)*(beta**TWO)/(EIGHT*gam*PI**TWO))*exp(ONE - r**TWO)
        rho = T**(ONE/(gam-ONE))
        p   = rho*T

        du = (beta/(TWO*PI))*exp((ONE-r**TWO)/TWO)*(-(y-self%yo))
        dv = (beta/(TWO*PI))*exp((ONE-r**TWO)/TWO)*((x-self%xo))
        u  = self%uinf + du
        v  = self%vinf + dv
        w  = ZERO





        select case (self%ivar)
            ! RHO
            case (1)
                calc = rho

            ! RHO-U
            case (2)
                calc = rho*u

            ! RHO-V
            case (3)
                calc = rho*v

            ! RHO-W
            case (4)
                calc = rho*w

            ! RHO-E
            case (5)
                calc = p/(gam-ONE)  +  HALF*rho*(u*u + v*v + w*w)

        end select

    end function


end module fcn_isentropic_vortex
