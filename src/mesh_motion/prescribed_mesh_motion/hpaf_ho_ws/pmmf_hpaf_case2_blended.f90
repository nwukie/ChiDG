module pmmf_hpaf_case2_blended
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FOUR, FIVE, EIGHT, PI
    use type_prescribed_mesh_motion_function,  only: prescribed_mesh_motion_function_t
    use mod_blending_functions,                 only: transition_polynomial_order5
    implicit none
    private








    !>  hpaf_case2_blended mesh motion. 
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017 
    !!
    !!
    !-------------------------------------------------------------------------------
    type, extends(prescribed_mesh_motion_function_t), public :: hpaf_case2_blended_pmmf
        private

        
    contains

        procedure   :: init
        procedure   :: compute_pos
        procedure   :: compute_vel

    end type hpaf_case2_blended_pmmf
    !********************************************************************************



contains




    !>
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017 
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(hpaf_case2_blended_pmmf),  intent(inout)  :: self

        !
        ! Set function name
        !
        call self%set_name("hpaf_case2_blended")


        !
        ! Set function options to default settings
        !
        call self%add_option('blend_rigid_radius', 50.0_rk) !Within this radius, the mesh is rigidly displaced
        call self%add_option('blend_transition_length', 30.0_rk) !Length scale where the mesh is blended between RB and fixed
    
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
        class(hpaf_case2_blended_pmmf),     intent(inout)  :: self
        real(rk),                       intent(in)  :: time
        real(rk),                  intent(in)  :: node(3)

        integer(ik)                                 :: ivar
        real(rk)                                    :: val(3)

        real(rk)                                    :: b1, b2, b3, &
                                                        A2, A3, x0, y0, x_ale, y_ale, &
                                                        xc, yc, height, theta, &
                                                        blend_rigid_radius, blend_transition_length, &
                                                        blend_dist, blend_val


        blend_rigid_radius = self%get_option('blend_rigid_radius')
        blend_transition_lenght = self%get_option('blend_transition_length')


        blend_dist = sqrt(val(1)**TWO+val(2)**TWO)-blend_rigid_radius

        if (blend_dist < ZERO) then
            blend_val = ZERO
        else if (blend_val > ONE) then
            blend_val = ONE
        else
            blend_val =  transition_polynomial_order5(blend_dist)
        end if

        b1 = time**TWO*(time**TWO-4._rk*time+4._rk)
        b2 = time**TWO*(3._rk-time)/4._rk
        b3 = time**THREE*(-8._rk*time**THREE+51._rk*time**TWO-111._rk*time+84._rk)/16._rk

!        !Case 1
!        height = b2
!        theta = ZERO
!
        !Case 2
        height = b2
        A2 = (60._rk*PI/180._rk)
        theta = A2*b1
        
!        !Case 3 
!        height = b3
!        A3 = (80._rk*PI/180._rk)
!        theta = A3*b1

        !Center of motion nodeinates
        xc = 1._rk/3._rk
        yc = 0._rk

        !Get the reference frame nodeinates of the grid point
        x0 = node(1)
        y0 = node(2)

        !Rotate about the center of motion
        x_ale =  cos(theta)*(x0-xc)+sin(theta)*(y0-yc) + xc
        y_ale = -sin(theta)*(x0-xc)+cos(theta)*(y0-yc) + yc

        !Translate vertically
        y_ale = y_ale + height



        val(1) = node(1)*blend_val+x_ale*(ONE-blend_val)
        val(2) = node(2)*blend_val+y_ale*(ONE-blend_val)
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
        class(hpaf_case2_blended_pmmf),     intent(inout)  :: self
        real(rk),                       intent(in)  :: time
        real(rk),                  intent(in)  :: node(3)

        integer(ik)                                 :: ivar
        real(rk)                                   :: val(3)
        real(rk)                                    :: b1, b2, b3, db1dt, db2dt, db3dt, &
                                                        A2, A3, x0, y0, x_ale, y_ale, &
                                                        xc, yc, height, theta, dheightdt, dthetadt, &
                                                        blend_rigid_radius, blend_transition_length, &
                                                        blend_dist, blend_val


        blend_rigid_radius = self%get_option('blend_rigid_radius')
        blend_transition_lenght = self%get_option('blend_transition_length')


        blend_dist = sqrt(val(1)**TWO+val(2)**TWO)-blend_rigid_radius

        if (blend_dist < ZERO) then
            blend_val = ZERO
        else if (blend_val > ONE) then
            blend_val = ONE
        else
            blend_val =  transition_polynomial_order5(blend_dist)
        end if



        b1 = time**TWO*(time**TWO-4._rk*time+4._rk)
        b2 = time**TWO*(3._rk-time)/4._rk
        b3 = time**THREE*(-8._rk*time**THREE+51._rk*time**TWO-111._rk*time+84._rk)/16._rk

        db1dt = TWO*time*(time**TWO-4._rk*time+4._rk) + &
                time**TWO*(TWO*time-4._rk)
        db2dt = TWO*time*(3._rk-time)/4._rk + &
                time**TWO*(-1._rk)/4._rk
        db3dt = THREE*time**TWO*(-8._rk*time**THREE+51._rk*time**TWO-111._rk*time+84._rk)/16._rk + &
                time**THREE*(-24._rk*time**TWO+102._rk*time-111._rk)/16._rk


!        !Case 1
!        height = b2
!        theta = ZERO
!        dheightdt = db2dt
!        dthetadt = ZERO
!
        !Case 2
        height = b2
        A2 = (60._rk*PI/180._rk)
        theta = A2*b1
        dheightdt = db2dt
        dthetadt = A2*db1dt
!
!        !Case 3 
!        height = b3
!        A3 = (80._rk*PI/180._rk)
!        theta = A3*b1
!        dheightdt = db3dt
!        dthetadt = A3*db1dt

        !Center of motion nodeinates
        xc = 1._rk/3._rk
        yc = 0._rk

        !Get the reference frame nodeinates of the grid point
        x0 = node(1)
        y0 = node(2)

        !Rotate about the center of motion
        x_ale = -sin(theta)*dthetadt*(x0-xc)+cos(theta)*dthetadt*(y0-yc) 
        y_ale = -cos(theta)*dthetadt*(x0-xc)-sin(theta)*dthetadt*(y0-yc)

        !Translate vertically
        y_ale = y_ale + dheightdt



        val(1) = (ONE-blend_val)*x_ale
        val(2) = (ONE-blend_val)*y_ale
        val(3) = ZERO 
 
        
    end function compute_vel
    !**********************************************************************************


end module pmmf_hpaf_case2_blended
