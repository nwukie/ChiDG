module fcn_roe_check
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_point,     only: point_t
    use atype_function, only: function_t
    implicit none
    private

    type, extends(function_t), public :: roe_check_f
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
    end type roe_check_f



contains

    subroutine set(self,valstring,val)
        class(roe_check_f),         intent(inout)   :: self
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
        class(roe_check_f), intent(in)   :: self
        integer(ik)                     :: order

        order = 3

    end function



    elemental function calc(self,pt)
        class(roe_check_f),  intent(in)  :: self
        type(point_t),       intent(in)  :: pt
        real(rk)                         :: calc

        real(rk)    :: x, y, z

        logical :: interior 

        x = pt%c1_
        y = pt%c2_
        z = pt%c3_

    

        interior = (x > 0.3333_rk) .and. (x < 0.6666_rk) .and. (y > 0.3333_rk) .and. (y < 0.6666_rk)

        select case (self%ivar)
            ! RHO
            case (1)
                if ( interior ) then
                    calc = 1.2_rk
                else
                    calc = 1.1_rk
                end if

            ! RHO-U
            case (2)
                if ( interior ) then
                    calc = 50._rk
                else
                    calc = 50._rk
                end if

            ! RHO-V
            case (3)
                if ( interior ) then
                    calc = 0._rk
                else
                    calc = 0._rk
                end if

            ! RHO-W
            case (4)
                if ( interior ) then
                    calc = 0._rk
                else
                    calc = 0._rk
                end if

            ! RHO-E
            case (5)
                if ( interior ) then
                    calc = 260000._rk
                else
                    calc = 260000._rk
                end if

        end select

    end function


end module fcn_roe_check
