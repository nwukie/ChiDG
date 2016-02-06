module mod_function
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use type_function,  only: function_t

    ! IMPORT FUNCTIONS
    use fcn_xsquared,           only: xsquared_f
    use fcn_ysquared,           only: ysquared_f
    use fcn_zsquared,           only: zsquared_f
    use fcn_xyz,                only: xyz_f
    use fcn_gaussian,           only: gaussian_f
    use fcn_constant,           only: constant_f
!    use fcn_isentropic_vortex,  only: isentropic_vortex_f
!    use fcn_sod_shock_tube,     only: sod_shock_tube_f
!    use fcn_roe_check,          only: roe_check_f
    implicit none

    ! LIST OF FUNCTION OBJECTS
    type(xsquared_f)            :: xsquared
    type(ysquared_f)            :: ysquared
    type(zsquared_f)            :: zsquared
    type(xyz_f)                 :: xyz
    type(gaussian_f)            :: gaussian
    type(constant_f)            :: constant
!    type(isentropic_vortex_f)   :: isentropic_vortex
!    type(sod_shock_tube_f)      :: sod_shock_tube
!    type(roe_check_f)           :: roe_check


contains


    !> Factory method for allocating concrete functions
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!  @param[inout]   fcn     Incoming function base class to be allocated
    !!  @param[in]      str     Function string for selecting concrete type
    !!
    !-------------------------------------------------------------------------------------
    subroutine create_function(fcn,str)
        class(function_t), allocatable, intent(inout)   :: fcn
        character(*),                   intent(in)      :: str


        !
        ! Check if function is already allocated. If so, reset.
        !
        if ( allocated(fcn) ) then
            deallocate(fcn)
        end if



        !
        ! Allocate function to concrete type
        !
        select case (trim(str))
            case ('xsquared')
                allocate(fcn,source=xsquared)

            case ('ysquared')
                allocate(fcn,source=ysquared)

            case ('zsquared')
                allocate(fcn,source=zsquared)

            case ('xyz')
                allocate(fcn,source=xyz)

            case ('gaussian')
                allocate(fcn,source=gaussian)

            case ('constant')
                allocate(fcn,source=constant)

!            case ('isentropic_vortex','isentropic vortex')
!                allocate(fcn,source=isentropic_vortex)
!
!            case ('sod','sod_shock_tube')
!                allocate(fcn,source=sod_shock_tube)
!
!            case ('roe_check')
!                allocate(fcn,source=roe_check)

            case default
                stop "Error: assign_function -- function string not recognized"

        end select


        !
        ! Call function initialization for options
        !
        call fcn%init()


    end subroutine create_function
    !**************************************************************************************




end module mod_function
