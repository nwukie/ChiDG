module pmmf_fic_case3
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FOUR, FIVE, EIGHT, PI
    use type_prescribed_mesh_motion_function,  only: prescribed_mesh_motion_function_t
    implicit none
    private



    !>  Flow In Cylinder - case3 prescribed motion. 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/02/2019 
    !!
    !-------------------------------------------------------------------------------
    type, extends(prescribed_mesh_motion_function_t), public :: fic_case3_pmmf
        private

    contains

        procedure   :: init
        procedure   :: compute_pos
        procedure   :: compute_vel

    end type fic_case3_pmmf
    !********************************************************************************



contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/02/2019 
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(fic_case3_pmmf),  intent(inout)  :: self

        ! Set function name
        call self%set_name("Flow In Cylinder - Case3")

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
        class(fic_case3_pmmf),  intent(inout)   :: self
        real(rk),               intent(in)      :: time
        real(rk),               intent(in)      :: node(3)

        real(rk)    :: val(3)
        real(rk)    :: Aa, alpha, psi, x0, y0, r0, theta0, &
                       a, b, e, r, x_ale, y_ale

        Aa    = 1.5
        alpha = time**TWO*(THREE-time)/FOUR
        psi   = ONE + (Aa - ONE)*alpha

        ! Get the reference frame nodeinates of the grid point
        x0     = node(1)
        y0     = node(2)
        r0     = sqrt(x0*x0 + y0*y0)
        theta0 = atan2(y0,x0)

        a = psi * r0
        b = r0 / psi
        e = sqrt(ONE - psi**(-FOUR))
        r = b/(sqrt(ONE - (e*cos(theta0))**TWO))


        x_ale = r*cos(theta0)
        y_ale = r*sin(theta0)


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
        class(fic_case3_pmmf),  intent(inout)   :: self
        real(rk),               intent(in)      :: time
        real(rk),               intent(in)      :: node(3)

        real(rk)    :: val(3)
        real(rk)    :: Aa, alpha, psi, x0, y0, r0, theta0, &
                       a, b, e, r, dalphadt, dpsidt, h, dhdt, drdt, &
                       vr, vx_ale, vy_ale, f, dedt, dbdt, dfdt

        Aa       = 1.5
        alpha    = time**TWO*(THREE-time)/FOUR
        psi      = ONE + (Aa - ONE)*alpha

        ! Get the reference frame nodeinates of the grid point
        x0     = node(1)
        y0     = node(2)
        r0     = sqrt(x0*x0 + y0*y0)
        theta0 = atan2(y0,x0)

        ! Location
        a = psi * r0
        b = r0 / psi
        e = sqrt(ONE - psi**(-FOUR))
        r = b/(sqrt(ONE - (e*cos(theta0))**TWO))
        f = sqrt(ONE - e*e*cos(theta0)*cos(theta0))


        ! Velocity 
        dalphadt = TWO*time*(THREE-time)/FOUR - (time**TWO)/FOUR
        dpsidt   = (Aa - ONE)*dalphadt

        dedt = TWO*dpsidt/((psi**5._rk)*sqrt(ONE - psi**(-FOUR)))
        dfdt = - cos(theta0)**TWO * e * dedt / (sqrt(ONE - e*e*cos(theta0)*cos(theta0)))
        dbdt = -r0*dpsidt/(psi*psi)

        drdt = (dbdt*f - b*dfdt)/(f*f)
!        h        = sqrt( (psi**FOUR)*(sin(theta0)**TWO) + cos(theta0)**TWO)
!        dhdt     = TWO*psi*psi*psi*dpsidt/h
!        drdt     = (r0*dpsidt*h - r0*psi*dhdt)/(h*h)
        vr       = drdt


        
        ! Convert to cartesian: assumes theta is not a function of (t)
        vx_ale = vr*cos(theta0)
        vy_ale = vr*sin(theta0)

        print*, sign(1._rk,vr), sign(1._rk,r - r0), vx_ale, vy_ale


        val(1) = vx_ale
        val(2) = vy_ale
        val(3) = ZERO 
 
        
    end function compute_vel
    !**********************************************************************************


end module pmmf_fic_case3
