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
    subroutine compute_periodic_offset( face, gq_node, offset_1, offset_2, offset_3 )
        type(face_t),   intent(in)      :: face
        type(point_t),  intent(in)      :: gq_node
        real(rk),       intent(inout)   :: offset_1
        real(rk),       intent(inout)   :: offset_2
        real(rk),       intent(inout)   :: offset_3




        !
        ! Define periodicity offset
        !
        if ( face%periodic_offset ) then

            offset_1 = face%chimera_offset_1
            offset_2 = face%chimera_offset_2
            offset_3 = face%chimera_offset_3

        else

            offset_1 = ZERO
            offset_2 = ZERO
            offset_3 = ZERO

        end if





    end subroutine compute_periodic_offset
    !*********************************************************************************






end module mod_periodic
