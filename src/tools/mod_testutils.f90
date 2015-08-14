module mod_testutils
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, ONE, TWO, THREE
    use type_point,     only: point_t



    implicit none


contains

    subroutine meshgen(string,pts)
        character(*),               intent(in)      :: string
        type(point_t), allocatable, intent(inout)   :: pts(:,:,:)


        select case (trim(string))
            case ('3x3x3','333')
                call meshgen_3x3x3_linear(pts)
            case default
                call signal(FATAL,'String identifying mesh generation routine was not recognized')
        end select


    end subroutine







    subroutine meshgen_3x3x3_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik), parameter      :: npt = 64
        integer(ik)                 :: ipt_xi, ipt_eta, ipt_zeta, ipt
        integer(ik)                 :: neqns, nterms_s, nterms_c
        real(rk), dimension(npt)    :: x,y,z

        !> elements (3x3x3) - linear
        !!
        !!            *-------*-------*-------*
        !!           /       /       /       /|
        !!          *-------*-------*-------* |
        !!         /       /       /       /| *
        !!        *-------*-------*-------* |/|
        !!       /       /       /       /| * |
        !!      *-------*-------*-------* |/| *
        !!      |       |       |       | * |/|
        !!      |       |       |       |/| * |
        !!      *-------*-------*-------* |/| *
        !!      |       |       |       | * |/
        !!      |       |       |       |/| *
        !!      *-------*-------*-------* |/
        !!      |       |       |       | *
        !!      |       |       |       |/
        !!      *-------*-------*-------*
        !!
        !!
        x = [ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE]

        y = [ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, THREE, THREE, THREE, THREE, &
             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, THREE, THREE, THREE, THREE, &
             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, THREE, THREE, THREE, THREE, &
             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, THREE, THREE, THREE, THREE]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, &
             TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, &
             THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE]


        !> Allocate point storage
        allocate(pts(4,4,4))

        ipt = 1
        do ipt_zeta = 1,4
            do ipt_eta = 1,4
                do ipt_xi = 1,4
                    call pts(ipt_xi,ipt_eta,ipt_zeta)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do

    end subroutine



end module mod_testutils
