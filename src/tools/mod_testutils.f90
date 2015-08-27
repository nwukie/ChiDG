module mod_testutils
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, ONE, TWO, THREE
    use type_point,     only: point_t



    implicit none


contains


    !>  Generate a set of points for a mesh. String input calls specialized
    !!  procedure for generating the points
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]      string   Character string used to select a specialized meshgen call
    !!  @param[inout]   pts      points_t array of rank-3 that gets allocated, filled, and returned
    !--------------------------------------------------------------------
    subroutine meshgen(string,pts)
        character(*),               intent(in)      :: string
        type(point_t), allocatable, intent(inout)   :: pts(:,:,:)


        select case (trim(string))
            case ('3x3x3','333')
                call meshgen_3x3x3_linear(pts)

            case ('2x2x2','222')
                call meshgen_2x2x2_linear(pts)

            case ('2x2x1','221')
                call meshgen_2x2x1_linear(pts)

            case ('3x3x1','331')
                call meshgen_3x3x1_linear(pts)

            case default
                call signal(FATAL,'String identifying mesh generation routine was not recognized')
        end select


    end subroutine






    !> Generate a set of points defining a 2x2x2 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_2x2x2_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik), parameter      :: npt = 27
        integer(ik)                 :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr
        real(rk), dimension(npt)    :: x,y,z

        ! elements (2x2x2) - linear
        !
        !          *-------*-------*
        !         /       /       /|
        !        *-------*-------* |
        !       /       /       /| *
        !      *-------*-------* |/|
        !      |       |       | * |
        !      |       |       |/| *
        !      *-------*-------* |/
        !      |       |       | *
        !      |       |       |/
        !      *-------*-------*
        !
        !
        x = [ZERO, ONE, TWO, ZERO, ONE, TWO, ZERO, ONE, TWO, &
             ZERO, ONE, TWO, ZERO, ONE, TWO, ZERO, ONE, TWO, &
             ZERO, ONE, TWO, ZERO, ONE, TWO, ZERO, ONE, TWO]

        y = [ZERO, ZERO, ZERO, ONE, ONE, ONE, TWO, TWO, TWO, &
             ZERO, ZERO, ZERO, ONE, ONE, ONE, TWO, TWO, TWO, &
             ZERO, ZERO, ZERO, ONE, ONE, ONE, TWO, TWO, TWO]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, &
             TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO]


        ! Allocate point storage
        allocate(pts(3,3,3), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        do ipt_zeta = 1,3
            do ipt_eta = 1,3
                do ipt_xi = 1,3
                    call pts(ipt_xi,ipt_eta,ipt_zeta)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do

    end subroutine













    !> Generate a set of points defining a 2x2x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_2x2x1_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik), parameter      :: npt = 18
        integer(ik)                 :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr
        real(rk), dimension(npt)    :: x,y,z

        ! elements (2x2x1) - linear
        !
        !          *-------*
        !         /       /|
        !        *-------* |
        !       /       /| * 
        !      *-------* |/|
        !      |       | * |
        !      |       |/| * 
        !      *-------* |/
        !      |       | * 
        !      |       |/ 
        !      *-------*
        !
        !
        x = [ZERO, ONE, TWO, ZERO, ONE, TWO, ZERO, ONE, TWO, &
             ZERO, ONE, TWO, ZERO, ONE, TWO, ZERO, ONE, TWO]
             

        y = [ZERO, ZERO, ZERO, ONE, ONE, ONE, TWO, TWO, TWO, &
             ZERO, ZERO, ZERO, ONE, ONE, ONE, TWO, TWO, TWO]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE]


        ! Allocate point storage
        allocate(pts(3,3,2), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        do ipt_zeta = 1,2
            do ipt_eta = 1,3
                do ipt_xi = 1,3
                    call pts(ipt_xi,ipt_eta,ipt_zeta)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do

    end subroutine















    !> Generate a set of points defining a 3x3x3 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_3x3x3_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik), parameter      :: npt = 64
        integer(ik)                 :: ipt_xi, ipt_eta, ipt_zeta, ipt
        real(rk), dimension(npt)    :: x,y,z

        ! elements (3x3x3) - linear
        !
        !            *-------*-------*-------*
        !           /       /       /       /|
        !          *-------*-------*-------* |
        !         /       /       /       /| *
        !        *-------*-------*-------* |/|
        !       /       /       /       /| * |
        !      *-------*-------*-------* |/| *
        !      |       |       |       | * |/|
        !      |       |       |       |/| * |
        !      *-------*-------*-------* |/| *
        !      |       |       |       | * |/
        !      |       |       |       |/| *
        !      *-------*-------*-------* |/
        !      |       |       |       | *
        !      |       |       |       |/
        !      *-------*-------*-------*
        !
        !
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


        ! Allocate point storage
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








    !> Generate a set of points defining a 3x3x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_3x3x1_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik), parameter      :: npt = 32
        integer(ik)                 :: ipt_xi, ipt_eta, ipt_zeta, ipt
        real(rk), dimension(npt)    :: x,y,z

        ! elements (3x3x1) - linear
        !
        !            *-------*
        !           /       /| 
        !          *-------* |
        !         /       /| *   
        !        *-------* |/|
        !       /       /| * | 
        !      *-------* |/| *
        !      |       | * |/| 
        !      |       |/| * |
        !      *-------* |/| *
        !      |       | * |/  
        !      |       |/| *  
        !      *-------* |/
        !      |       | *
        !      |       |/ 
        !      *-------*
        !
        !
        x = [ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE]

        y = [ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, THREE, THREE, THREE, THREE, &
             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, THREE, THREE, THREE, THREE]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE]


        ! Allocate point storage
        allocate(pts(4,4,2))

        ipt = 1
        do ipt_zeta = 1,2
            do ipt_eta = 1,4
                do ipt_xi = 1,4
                    call pts(ipt_xi,ipt_eta,ipt_zeta)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do

    end subroutine






















end module mod_testutils
