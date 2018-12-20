!
! This PMMF is used in the convecting vortex regression test.
!
module pmmf_sinusoidal_convecting_vortex
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FOUR, FIVE, EIGHT, PI
    use type_prescribed_mesh_motion_function,  only: prescribed_mesh_motion_function_t
    implicit none
    private








    !>  sinusoidal_convecting_vortex mesh motion. 
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017 
    !!
    !!
    !-------------------------------------------------------------------------------
    type, extends(prescribed_mesh_motion_function_t), public :: sinusoidal_convecting_vortex_pmmf
        private

        
    contains

        procedure   :: init
        procedure   :: compute_pos
        procedure   :: compute_vel

    end type sinusoidal_convecting_vortex_pmmf
    !********************************************************************************



contains




    !>
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017 
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(sinusoidal_convecting_vortex_pmmf),  intent(inout)  :: self

        !
        ! Set function name
        !
        call self%set_name("sinusoidal_convecting_vortex")


        !
        ! Set function options to default settings
        !
        call self%add_option('L_X',20._rk)
        call self%add_option('L_Y',15._rk)
        call self%add_option('L_Z',1._rk)
        call self%add_option('GRID_MODE_X',1._rk)
        call self%add_option('GRID_MODE_Y',1._rk)
        call self%add_option('GRID_MODE_Z',1._rk)
        call self%add_option('GRID_FREQ_X',TWO*PI)
        call self%add_option('GRID_FREQ_Y',FOUR*PI)
        call self%add_option('GRID_FREQ_Z',1._rk)
        call self%add_option('GRID_AMP_X',2._rk)
        call self%add_option('GRID_AMP_Y',1.5_rk)
        call self%add_option('GRID_AMP_Z',0._rk)


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
        class(sinusoidal_convecting_vortex_pmmf),     intent(inout)   :: self
        real(rk),                   intent(in)      :: time
        real(rk),                   intent(in)      :: node(3)

        integer(ik)                                 :: ivar
        real(rk)                                    :: val(3)

        real(rk)                                    :: L_X, L_Y, L_Z, &
                                                        GRID_MODE_X, GRID_MODE_Y, GRID_MODE_Z, &
                                                        GRID_FREQ_X, GRID_FREQ_Y, GRID_FREQ_Z, &
                                                        GRID_AMP_X, GRID_AMP_Y, GRID_AMP_Z
        L_X = self%get_option_value('L_X')
        L_Y = self%get_option_value('L_Y')
        L_Z = self%get_option_value('L_Z')
        GRID_MODE_X = self%get_option_value('GRID_MODE_X')
        GRID_MODE_Y = self%get_option_value('GRID_MODE_Y')
        GRID_MODE_Z = self%get_option_value('GRID_MODE_Z')
        GRID_FREQ_X = self%get_option_value('GRID_FREQ_X')
        GRID_FREQ_Y = self%get_option_value('GRID_FREQ_Y')
        GRID_FREQ_Z = self%get_option_value('GRID_FREQ_Z')
        GRID_AMP_X = self%get_option_value('GRID_AMP_X')
        GRID_AMP_Y = self%get_option_value('GRID_AMP_Y')
        GRID_AMP_Z = self%get_option_value('GRID_AMP_Z')

        val(1) = node(1) + &
            GRID_AMP_X* &
            sin(GRID_MODE_X*TWO*PI*node(1)/L_X)* &
            sin(GRID_MODE_Y*TWO*PI*node(2)/L_Y)* &
            sin(GRID_FREQ_X*time)
        
        val(2) = node(2) + &
            GRID_AMP_Y* &
            sin(GRID_MODE_X*TWO*PI*node(1)/L_X)* &
            sin(GRID_MODE_Y*TWO*PI*node(2)/L_Y)* &
            sin(GRID_FREQ_Y*time)

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
        class(sinusoidal_convecting_vortex_pmmf),     intent(inout)   :: self
        real(rk),                   intent(in)      :: time
        real(rk),                   intent(in)      :: node(3)

        integer(ik)                                 :: ivar
        real(rk)                                    :: val(3)
        real(rk)                                    :: L_X, L_Y, L_Z, &
                                                        GRID_MODE_X, GRID_MODE_Y, GRID_MODE_Z, &
                                                        GRID_FREQ_X, GRID_FREQ_Y, GRID_FREQ_Z, &
                                                        GRID_AMP_X, GRID_AMP_Y, GRID_AMP_Z
        L_X = self%get_option_value('L_X')
        L_Y = self%get_option_value('L_Y')
        L_Z = self%get_option_value('L_Z')
        GRID_MODE_X = self%get_option_value('GRID_MODE_X')
        GRID_MODE_Y = self%get_option_value('GRID_MODE_Y')
        GRID_MODE_Z = self%get_option_value('GRID_MODE_Z')
        GRID_FREQ_X = self%get_option_value('GRID_FREQ_X')
        GRID_FREQ_Y = self%get_option_value('GRID_FREQ_Y')
        GRID_FREQ_Z = self%get_option_value('GRID_FREQ_Z')
        GRID_AMP_X = self%get_option_value('GRID_AMP_X')
        GRID_AMP_Y = self%get_option_value('GRID_AMP_Y')
        GRID_AMP_Z = self%get_option_value('GRID_AMP_Z')
        
        val(1) =  &
            GRID_AMP_X* &
            sin(GRID_MODE_X*TWO*PI*node(1)/L_X)* &
            sin(GRID_MODE_Y*TWO*PI*node(2)/L_Y)* &
            GRID_FREQ_X*cos(GRID_FREQ_X*time)
        
        val(2) =  &
            GRID_AMP_Y* &
            sin(GRID_MODE_X*TWO*PI*node(1)/L_X)* &
            sin(GRID_MODE_Y*TWO*PI*node(2)/L_Y)* &
            GRID_FREQ_Y*cos(GRID_FREQ_Y*time)

        
        val(3) = ZERO
    end function compute_vel
    !**********************************************************************************


end module pmmf_sinusoidal_convecting_vortex
