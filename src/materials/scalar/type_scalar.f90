module type_scalar
#include <messenger.h>
    use mod_kinds,      only: rk
    use type_material,  only: material_t
    use DNAD_D
    implicit none
    private





    !>  Abstract scalar class
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/15/2016
    !!
    !--------------------------------------------------------------------------
    type, public, extends(material_t), abstract :: scalar_t


    contains

        procedure   :: compute_mu
        procedure   :: compute_cx
        procedure   :: compute_cy
        procedure   :: compute_cz

    end type scalar_t
    !**************************************************************************


contains




!    !>
!    !!
!    !!  @author Nathan A. Wukie (AFRL)
!    !!  @date   8/30/2016
!    !!
!    !!
!    !------------------------------------------------------------------------------------------
!    subroutine init(self)
!        class(fluid_t), intent(inout)   :: self
!
!
!        call self%contribute_equation("Density")
!        call self%contribute_equation("X-Momentum")
!        call self%contribute_equation("Y-Momentum")
!        call self%contribute_equation("Z-Momentum")
!        call self%contribute_equation("Energy")
!        call self%contribute_equation("Turbulence Viscosity")
!        call self%contribute_equation("Turbulence Kinetic Energy")
!        call self%contribute_equation("Turbulence Dissipation Rate")
!
!
!    end subroutine init
!    !******************************************************************************************

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/21/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    impure elemental function compute_mu(self,u) result(val)
        class(scalar_t),    intent(in)  :: self
        type(AD_D),         intent(in)  :: u

        type(AD_D)                      :: val
        character(len=:),   allocatable :: msg

        msg = "scalar%compute_mu: There wasn't any compute_mu routine implemented in the &
                                  allocated scalar model. Make sure this routine gets implemented &
                                  if there are diffusive scalar operators being computed"
        call chidg_signal(FATAL,msg)

    end function compute_mu
    !*******************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/21/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    impure elemental function compute_cx(self,u) result(val)
        class(scalar_t),    intent(in)  :: self
        type(AD_D),         intent(in)  :: u

        type(AD_D)                      :: val
        character(len=:),   allocatable :: msg

        msg = "scalar%compute_cx: There wasn't any compute_cx routine implemented in the &
                                  allocated scalar model. Make sure this routine gets implemented &
                                  if there are diffusive scalar operators being computed"
        call chidg_signal(FATAL,msg)

    end function compute_cx
    !*******************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/21/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    impure elemental function compute_cy(self,u) result(val)
        class(scalar_t),    intent(in)  :: self
        type(AD_D),         intent(in)  :: u

        type(AD_D)                      :: val
        character(len=:),   allocatable :: msg

        msg = "scalar%compute_cy: There wasn't any compute_cy routine implemented in the &
                                  allocated scalar model. Make sure this routine gets implemented &
                                  if there are diffusive scalar operators being computed"
        call chidg_signal(FATAL,msg)

    end function compute_cy
    !*******************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/21/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    impure elemental function compute_cz(self,u) result(val)
        class(scalar_t),    intent(in)  :: self
        type(AD_D),         intent(in)  :: u

        type(AD_D)                      :: val
        character(len=:),   allocatable :: msg

        msg = "scalar%compute_cz: There wasn't any compute_cz routine implemented in the &
                                  allocated scalar model. Make sure this routine gets implemented &
                                  if there are diffusive scalar operators being computed"
        call chidg_signal(FATAL,msg)


    end function compute_cz
    !*******************************************************************************







end module type_scalar
