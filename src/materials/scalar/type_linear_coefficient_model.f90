module type_linear_coefficient_model
#include <messenger.h>
    use mod_constants,  only: ONE, ZERO
    use type_scalar,    only: scalar_t
    use DNAD_D
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/21/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(scalar_t) :: linear_coefficient_model_t

    contains

        procedure   :: compute_mu
        procedure   :: compute_cx
        procedure   :: compute_cy
        procedure   :: compute_cz

    end type linear_coefficient_model_t
    !******************************************************************************







contains




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/21/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    impure elemental function compute_mu(self,u,dudx,dudy,dudz) result(val)
        class(linear_coefficient_model_t),  intent(in)  :: self
        type(AD_D),                         intent(in)  :: u
        type(AD_D),         intent(in), optional    :: dudx
        type(AD_D),         intent(in), optional    :: dudy
        type(AD_D),         intent(in), optional    :: dudz

        type(AD_D) :: val

        ! Allocate vals
        val =u

        ! Set to constant
        val = ONE

    end function compute_mu
    !*******************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/21/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    impure elemental function compute_cx(self,u,dudx,dudy,dudz) result(val)
        class(linear_coefficient_model_t),  intent(in)  :: self
        type(AD_D),                         intent(in)  :: u
        type(AD_D),         intent(in), optional    :: dudx
        type(AD_D),         intent(in), optional    :: dudy
        type(AD_D),         intent(in), optional    :: dudz

        type(AD_D)  :: val

        ! Allocate vals
        val =u

        ! Set to constant
        val = ONE

    end function compute_cx
    !*******************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/21/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    impure elemental function compute_cy(self,u,dudx,dudy,dudz) result(val)
        class(linear_coefficient_model_t),  intent(in)  :: self
        type(AD_D),                         intent(in)  :: u
        type(AD_D),         intent(in), optional    :: dudx
        type(AD_D),         intent(in), optional    :: dudy
        type(AD_D),         intent(in), optional    :: dudz

        type(AD_D) :: val

        ! Allocate vals
        val =u

        ! Set to constant
        val = ZERO

    end function compute_cy
    !*******************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/21/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    impure elemental function compute_cz(self,u,dudx,dudy,dudz) result(val)
        class(linear_coefficient_model_t),  intent(in)  :: self
        type(AD_D),                         intent(in)  :: u
        type(AD_D),         intent(in), optional    :: dudx
        type(AD_D),         intent(in), optional    :: dudy
        type(AD_D),         intent(in), optional    :: dudz

        type(AD_D)  :: val

        ! Allocate vals
        val =u

        ! Set to constant
        val = ZERO

    end function compute_cz
    !*******************************************************************************


end module type_linear_coefficient_model
