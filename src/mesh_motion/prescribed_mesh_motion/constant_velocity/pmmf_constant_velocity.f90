module pmmf_constant_velocity
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_prescribed_mesh_motion_function,  only: prescribed_mesh_motion_function_t
    implicit none
    private





    !>  Constant mesh velocity (no mesh displacement).
    !!
    !!  @author Nathan A. Wukie
    !!  @date  8/15/2018
    !!
    !!
    !-------------------------------------------------------------------------------
    type, extends(prescribed_mesh_motion_function_t), public :: constant_velocity_pmmf
        private
        
    contains

        procedure   :: init
        procedure   :: compute_pos
        procedure   :: compute_vel

    end type constant_velocity_pmmf
    !********************************************************************************



contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date  8/15/2018 
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(constant_velocity_pmmf),  intent(inout)  :: self

        ! Set function name
        call self%set_name("Constant Velocity")

    end subroutine init
    !*************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date  8/15/2018 
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function compute_pos(self,time,node) result(val)
        class(constant_velocity_pmmf),  intent(inout)   :: self
        real(rk),                       intent(in)      :: time
        real(rk),                       intent(in)      :: node(3)

        real(rk) :: val(3)

        val = node
        
    end function compute_pos
    !**********************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date  8/15/2018 
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function compute_vel(self,time,node) result(val)
        class(constant_velocity_pmmf),  intent(inout)   :: self
        real(rk),                       intent(in)      :: time
        real(rk),                       intent(in)      :: node(3)

        real(rk)    :: val(3)

        val = [ZERO, 10._rk, ZERO]
        
    end function compute_vel
    !**********************************************************************************


end module pmmf_constant_velocity
