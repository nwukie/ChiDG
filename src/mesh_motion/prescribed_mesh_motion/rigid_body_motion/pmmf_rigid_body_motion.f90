module pmmf_rigid_body_motion
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use mod_rigid_body_motion,      only: t0, t1, rigid_body_motion_disp_old, rigid_body_motion_disp_new, rigid_body_motion_vel
    use type_prescribed_mesh_motion_function,  only: prescribed_mesh_motion_function_t
    implicit none
    private








    !>  rigid_body_motion mesh (no mesh motion).
    !!
    !!  @author Eric Wolf
    !!  @date  3/15/2017 
    !!
    !!
    !-------------------------------------------------------------------------------
    type, extends(prescribed_mesh_motion_function_t), public :: rigid_body_motion_pmmf
        private

        
    contains

        procedure   :: init
        procedure   :: compute_pos
        procedure   :: compute_vel

    end type rigid_body_motion_pmmf
    !********************************************************************************



contains




    !>
    !!
    !!  @author Eric Wolf
    !!  @date  3/15/2017 
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(rigid_body_motion_pmmf),  intent(inout)  :: self

        !
        ! Set function name
        !
        !self%name = "isentropic vortex  ::  weeeee!"
        call self%set_name("rigid_body_motion")




    end subroutine init
    !*************************************************************************



    !>
    !!
    !!  @author Eric Wolf
    !!  @date  3/15/2017 
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function compute_pos(self,time,node) result(val)
        class(rigid_body_motion_pmmf),     intent(inout)   :: self
        real(rk),               intent(in)      :: time
        real(rk),               intent(in)      :: node(3)

        integer(ik)                             :: ivar
        real(rk)                                :: val(3)

        val = node + (t1-time)/(t1-t0)*rigid_body_motion_disp_new + (time-t0)/(t1-t0)*rigid_body_motion_disp_new
        

        
    end function compute_pos
    !**********************************************************************************






    !>
    !!
    !!  @author Eric Wolf
    !!  @date  3/15/2017 
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function compute_vel(self,time,node) result(val)
        class(rigid_body_motion_pmmf),     intent(inout)   :: self
        real(rk),               intent(in)      :: time
        real(rk),               intent(in)      :: node(3)

        integer(ik)                             :: ivar
        real(rk)                                :: val(3)

        val = rigid_body_motion_vel
        
    end function compute_vel
    !**********************************************************************************


end module pmmf_rigid_body_motion
