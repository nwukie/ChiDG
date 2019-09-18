module pmmf_constant_motion_1d_x
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_prescribed_mesh_motion_function,  only: prescribed_mesh_motion_function_t
    implicit none
    private








    !>  constant_motion_1d_x mesh motion. 
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017 
    !!
    !!
    !-------------------------------------------------------------------------------
    type, extends(prescribed_mesh_motion_function_t), public :: constant_motion_1d_x_pmmf
        private

        
    contains

        procedure   :: init
        procedure   :: compute_pos
        procedure   :: compute_vel

    end type constant_motion_1d_x_pmmf
    !********************************************************************************



contains




    !>
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017 
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(constant_motion_1d_x_pmmf),  intent(inout)  :: self

        !
        ! Set function name
        !
        call self%set_name("constant_motion_1d_x")


        !
        ! Set function options to default settings
        !
        call self%add_option('grid_advection_velocity',ONE)
        
    end subroutine init
    !*************************************************************************



    !>
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017 
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function compute_pos(self,time,node) result(val)
        class(constant_motion_1d_x_pmmf),     intent(inout)   :: self
        real(rk),                   intent(in)      :: time
        real(rk),                   intent(in)      :: node(3)

        integer(ik)                                 :: ivar
        real(rk)                                    :: val(3)

        real(rk)                                    :: u_grid 

        u_grid = self%get_option_value('grid_advection_velocity')
        val(1) = node(1) + time*u_grid 
       
        val(2) = node(2)
        val(3) = node(3)
        
    end function compute_pos
    !**********************************************************************************






    !>
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017 
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function compute_vel(self,time,node) result(val)
        class(constant_motion_1d_x_pmmf),     intent(inout)   :: self
        real(rk),                   intent(in)      :: time
        real(rk),                   intent(in)      :: node(3)

        integer(ik)                                 :: ivar
        real(rk)                                    :: val(3)
        real(rk)                                    :: u_grid 

        u_grid = self%get_option_value('grid_advection_velocity')

        val(1) = u_grid 
        val(2) = ZERO
        val(3) = ZERO

        
    end function compute_vel
    !**********************************************************************************


end module pmmf_constant_motion_1d_x
