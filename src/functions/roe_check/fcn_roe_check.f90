module fcn_roe_check
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_point,     only: point_t
    use type_function,  only: function_t
    implicit none
    private








    !>  @TODO: BROKEN Needs fixed
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(function_t), public :: roe_check_f
        private


        integer(ik) :: ivar

    contains

        procedure   :: init
        procedure   :: compute

    end type roe_check_f
    !*******************************************************************************



contains

    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(roe_check_f),  intent(inout)  :: self

        !
        ! Set function name
        !
        self%name = "roe check  ::  "



    end subroutine init
    !*************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    elemental function compute(self,time,coord) result(val)
        class(roe_check_f),     intent(in)  :: self
        real(rk),               intent(in)  :: time
        type(point_t),          intent(in)  :: coord

        real(rk)                            :: val

        real(rk)    :: x, y, z

        logical :: interior 

        x = coord%c1_
        y = coord%c2_
        z = coord%c3_

    

        interior = (x > 0.3333_rk) .and. (x < 0.6666_rk) .and. (y > 0.3333_rk) .and. (y < 0.6666_rk)

        select case (self%ivar)
            ! RHO
            case (1)
                if ( interior ) then
                    val = 1.2_rk
                else
                    val = 1.1_rk
                end if

            ! RHO-U
            case (2)
                if ( interior ) then
                    val = 50._rk
                else
                    val = 50._rk
                end if

            ! RHO-V
            case (3)
                if ( interior ) then
                    val = 0._rk
                else
                    val = 0._rk
                end if

            ! RHO-W
            case (4)
                if ( interior ) then
                    val = 0._rk
                else
                    val = 0._rk
                end if

            ! RHO-E
            case (5)
                if ( interior ) then
                    val = 260000._rk
                else
                    val = 260000._rk
                end if

        end select

    end function compute
    !********************************************************************************************


end module fcn_roe_check
