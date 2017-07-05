module pmmf_static
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
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
        class(static_pmmf),     intent(inout)   :: self
        real(rk),               intent(in)      :: time
        real(rk),               intent(in)      :: node(3)

        integer(ik)                             :: ivar
        real(rk)                                :: val(3)

        val = node
        
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
        class(static_pmmf),     intent(inout)   :: self
        real(rk),               intent(in)      :: time
        real(rk),               intent(in)      :: node(3)

        integer(ik)                             :: ivar
        real(rk)                                :: val(3)

        val(1) = ZERO
        val(2) = ZERO
        val(3) = ZERO
        
    end function compute_vel
    !**********************************************************************************


end module pmmf_static
