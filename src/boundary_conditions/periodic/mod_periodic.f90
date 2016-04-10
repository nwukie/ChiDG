module mod_periodic
#include <messenger.h>
    use mod_kinds,      only: rk
    use mod_constants,  only: ZERO, TWO, PI
    use type_point,     only: point_t
    use type_face,      only: face_t
    implicit none



contains



    !>  Compute cartesian coordinate offset for periodic-chimera faces.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/31/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine compute_periodic_offset( face, gq_node, offset_x, offset_y, offset_z )
        type(face_t),   intent(in)      :: face
        type(point_t),  intent(in)      :: gq_node
        real(rk),       intent(inout)   :: offset_x
        real(rk),       intent(inout)   :: offset_y
        real(rk),       intent(inout)   :: offset_z

        real(rk)    :: x1, y1, r1, theta1, x2, y2, theta2



        !
        ! Periodicity offset methods
        !
        if ( allocated(face%periodic_type) ) then



            if ( face%periodic_type == 'cartesian' ) then
                !print*, 'Computing cartesian periodic offset'
                offset_x = face%chimera_offset_x
                offset_y = face%chimera_offset_y
                offset_z = face%chimera_offset_z



            else if ( face%periodic_type == 'cylindrical' ) then
                !print*, 'Computing cylindrical periodic offset'
                x1 = gq_node%c1_
                y1 = gq_node%c2_
                r1     = sqrt(x1**TWO + y1**TWO)
                theta1 = atan2(y1,x1)

                theta2 = theta1 + (PI/180._rk)*face%chimera_offset_theta
                x2     = r1*cos(theta2)
                y2     = r1*sin(theta2)


                offset_x = x2-x1
                offset_y = y2-y1
                offset_z = ZERO


            else

                call chidg_signal(FATAL,"mod_periodic::compute_periodic_offset - Face has invalid value for face%periodic_type")

            end if





        else

            !
            ! No periodic type was allocated so no offset.
            !
            !print*, 'No periodic offset'
            offset_x = ZERO
            offset_y = ZERO
            offset_z = ZERO

        end if





    end subroutine compute_periodic_offset
    !*********************************************************************************






end module mod_periodic
