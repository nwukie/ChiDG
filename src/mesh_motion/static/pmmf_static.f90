module pmmf_static
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_point,     only: point_t
    use type_prescribed_mesh_motion_function,  only: prescribed_mesh_motion_function_t
    implicit none
    private








    !>  Static mesh (no mesh motion).
    !!
    !!  @author Eric Wolf
    !!  @date  3/15/2017 
    !!
    !!
    !-------------------------------------------------------------------------------
    type, extends(prescribed_mesh_motion_function_t), public :: static_pmmf
        private

        
    contains

        procedure   :: init
        procedure   :: compute_pos
        procedure   :: compute_vel

    end type static_pmmf
    !********************************************************************************



contains




    !>
    !!
    !!  @author Eric Wolf
    !!  @date  3/15/2017 
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(static_pmmf),  intent(inout)  :: self

        !
        ! Set function name
        !
        !self%name = "isentropic vortex  ::  weeeee!"
        call self%set_name("static")


        !
        ! Set function options to default settings
        !
!        call self%dict%set('uinf',2._rk/sqrt(5._rk))
!        call self%dict%set('vinf',1._rk/sqrt(5._rk))
!        call self%dict%set('winf',0._rk)
!        call self%dict%set('Minf',1._rk*sqrt(
!        call self%dict%set('beta',0._rk)
!        call self%dict%set('xo',1._rk)
!        call self%dict%set('yo',1._rk)
!        call self%dict%set('zo',1._rk)
!

    end subroutine init
    !*************************************************************************



    !>
    !!
    !!  @author Eric Wolf
    !!  @date  3/15/2017 
    !!
    !!
    !-----------------------------------------------------------------------------------------
    impure  elemental function compute_pos(self,time,coord) result(val)
        class(static_pmmf),     intent(inout)  :: self
        real(rk),                       intent(in)  :: time
        type(point_t),                  intent(in)  :: coord

        integer(ik)                                 :: ivar
        type(point_t)                                    :: val

        val = coord
        
    end function compute_pos
    !**********************************************************************************






    !>
    !!
    !!  @author Eric Wolf
    !!  @date  3/15/2017 
    !!
    !!
    !-----------------------------------------------------------------------------------------
    impure  elemental function compute_vel(self,time,coord) result(val)
        class(static_pmmf),     intent(inout)  :: self
        real(rk),                       intent(in)  :: time
        type(point_t),                  intent(in)  :: coord

        integer(ik)                                 :: ivar
        type(point_t)                                   :: val

        val%c1_ = ZERO
        val%c2_ = ZERO
        val%c3_ = ZERO
        
    end function compute_vel
    !**********************************************************************************


end module pmmf_static
