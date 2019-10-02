module pmmf_fic_case2
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FOUR, FIVE, EIGHT, PI
    use type_prescribed_mesh_motion_function,  only: prescribed_mesh_motion_function_t
    implicit none
    private



    !>  Flow In Cylinder - Case2 prescribed motion. 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/02/2019 
    !!
    !-------------------------------------------------------------------------------
    type, extends(prescribed_mesh_motion_function_t), public :: fic_case2_pmmf
        private

    contains

        procedure   :: init
        procedure   :: compute_pos
        procedure   :: compute_vel

    end type fic_case2_pmmf
    !********************************************************************************



contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/02/2019 
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(fic_case2_pmmf),  intent(inout)  :: self

        ! Set function name
        call self%set_name("Flow In Cylinder - Case2")

    end subroutine init
    !*************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/02/2019
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function compute_pos(self,time,node) result(val)
        class(fic_case2_pmmf),  intent(inout)   :: self
        real(rk),               intent(in)      :: time
        real(rk),               intent(in)      :: node(3)

        integer(ik) :: ivar
        real(rk)    :: val(3)

        real(rk)    :: alpha, Atheta, x0, y0, x_ale, y_ale, &
                       xc, yc, height, theta

        Atheta = PI
        alpha = time**TWO*(THREE-time)/FOUR

        !Case 2
        height = ZERO
        theta  = Atheta*alpha

        !Center of motion nodeinates
        xc = ZERO
        yc = ZERO

        !Get the reference frame nodeinates of the grid point
        x0 = node(1)
        y0 = node(2)


        !Rotation about center
        x_ale =  cos(theta)*(x0-xc)+sin(theta)*(y0-yc) + xc
        y_ale = -sin(theta)*(x0-xc)+cos(theta)*(y0-yc) + yc

        !Translation
        y_ale = y_ale + height



        val(1) = x_ale
        val(2) = y_ale
        val(3) = node(3)
        
    end function compute_pos
    !**********************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/02/2019
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function compute_vel(self,time,node) result(val)
        class(fic_case2_pmmf),  intent(inout)   :: self
        real(rk),               intent(in)      :: time
        real(rk),               intent(in)      :: node(3)

        integer(ik) :: ivar
        real(rk)    :: val(3)
        real(rk)    :: alpha, dalphadt, &
                       Atheta, x0, y0, vx_ale, vy_ale, &
                       xc, yc, height, theta, dheightdt, dthetadt


        Atheta   = PI
        alpha    = time**TWO*(THREE-time)/FOUR
        dalphadt = TWO*time*(THREE-time)/FOUR + &
                   time**TWO*(-ONE)/FOUR


        !Case 2
        height    = ZERO
        theta     = Atheta*alpha
        dheightdt = ZERO
        dthetadt  = Atheta*dalphadt


        !Center of motion nodeinates
        xc = ZERO
        yc = ZERO


        !Get the reference frame nodeinates of the grid point
        x0 = node(1)
        y0 = node(2)


        ! Rotation about (xc,yc)
        vx_ale = -sin(theta)*dthetadt*(x0-xc)+cos(theta)*dthetadt*(y0-yc) 
        vy_ale = -cos(theta)*dthetadt*(x0-xc)-sin(theta)*dthetadt*(y0-yc)

        !Translate vertically
        vy_ale = vy_ale + dheightdt

        val(1) = vx_ale
        val(2) = vy_ale
        val(3) = ZERO 
 
        
    end function compute_vel
    !**********************************************************************************


end module pmmf_fic_case2
