module mod_testutils
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX
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
            case ('1x1x1','111')
                call meshgen_1x1x1_linear(pts)

            case ('1x1x1_unit','111u')
                call meshgen_1x1x1_unit_linear(pts)

            case ('3x3x3','333')
                call meshgen_3x3x3_linear(pts)

            case ('3x3x3_unit','333u')
                call meshgen_3x3x3_unit_linear(pts)

            case ('2x2x2','222')
                call meshgen_2x2x2_linear(pts)

            case ('2x2x1','221')
                call meshgen_2x2x1_linear(pts)

            case ('3x3x1','331')
                call meshgen_3x3x1_linear(pts)

            case ('4x1x1','411')
                call meshgen_4x1x1_linear(pts)

            case ('3x1x1','311')
                call meshgen_3x1x1_linear(pts)

            case ('2x1x1','211')
                call meshgen_2x1x1_linear(pts)

            case ('40x15x1')
                call meshgen_40x15x1_linear(pts)

            case ('15x15x1')
                call meshgen_15x15x1_linear(pts)

            case ('15x15x2')
                call meshgen_15x15x2_linear(pts)

            case ('15x15x3')
                call meshgen_15x15x3_linear(pts)

            case default
                call chidg_signal(FATAL,'String identifying mesh generation routine was not recognized')
        end select


    end subroutine











    !> Generate a set of points defining a 1x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_1x1x1_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, npts_y, npts_z
        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (1x1x1) - linear

        npts_x = 2
        npts_y = 2
        npts_z = 2

        dx = 1._rk
        dy = 1._rk
        dz = 1._rk

        ! Allocate point storage
        allocate(pts(npts_x,npts_y,npts_z), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        z = ZERO
        do ipt_zeta = 1,npts_z
            y = ZERO

            do ipt_eta = 1,npts_y
                x = ZERO

                do ipt_xi = 1,npts_x
                    call pts(ipt_xi,ipt_eta,ipt_zeta)%set(x, y, z)
                    ipt = ipt + 1
                    x = x + dx
                end do

                y = y + dy
            end do

            z = z + dz
        end do

    end subroutine














    !> Generate a set of points defining a 1x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_1x1x1_unit_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, npts_y, npts_z
        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (1x1x1) - linear

        npts_x = 2
        npts_y = 2
        npts_z = 2

        dx = 2._rk
        dy = 2._rk
        dz = 2._rk

        ! Allocate point storage
        allocate(pts(npts_x,npts_y,npts_z), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        z = -ONE
        do ipt_zeta = 1,npts_z
            y = -ONE

            do ipt_eta = 1,npts_y
                x = -ONE

                do ipt_xi = 1,npts_x
                    call pts(ipt_xi,ipt_eta,ipt_zeta)%set(x, y, z)
                    ipt = ipt + 1
                    x = x + dx
                end do

                y = y + dy
            end do

            z = z + dz
        end do

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
















    !> Generate a set of points defining a 3x3x3 unit-element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_3x3x3_unit_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik), parameter      :: npt = 64
        integer(ik)                 :: ipt_xi, ipt_eta, ipt_zeta, ipt
        real(rk), dimension(npt)    :: x,y,z

        ! elements (3x3x3) - unit linear
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
        x = [ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, &
             ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, &
             ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, &
             ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX]

        y = [ZERO, ZERO, ZERO, ZERO, TWO, TWO, TWO, TWO, FOUR, FOUR, FOUR, FOUR, SIX, SIX, SIX, SIX, &
             ZERO, ZERO, ZERO, ZERO, TWO, TWO, TWO, TWO, FOUR, FOUR, FOUR, FOUR, SIX, SIX, SIX, SIX, &
             ZERO, ZERO, ZERO, ZERO, TWO, TWO, TWO, TWO, FOUR, FOUR, FOUR, FOUR, SIX, SIX, SIX, SIX, &
             ZERO, ZERO, ZERO, ZERO, TWO, TWO, TWO, TWO, FOUR, FOUR, FOUR, FOUR, SIX, SIX, SIX, SIX]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
             TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, &
             FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, &
             SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX]


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
















    !> Generate a set of points defining a 4x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_4x1x1_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik), parameter      :: npt = 20
        integer(ik)                 :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr
        real(rk), dimension(npt)    :: x,y,z

        ! elements (4x1x1) - linear
        !
        !      *------*------*------*------*
        !      |      |      |      |      | 
        !      |      |      |      |      | 
        !      *------*------*------*------*
        !



        x = [ZERO, ONE, TWO, THREE, FOUR, &       
             ZERO, ONE, TWO, THREE, FOUR, &       
             ZERO, ONE, TWO, THREE, FOUR, &       
             ZERO, ONE, TWO, THREE, FOUR]
             

        y = [ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, &
             ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, &
             ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, &
             ONE, ONE, ONE, ONE, ONE]



        ! Allocate point storage
        allocate(pts(5,2,2), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        do ipt_zeta = 1,2
            do ipt_eta = 1,2
                do ipt_xi = 1,5
                    call pts(ipt_xi,ipt_eta,ipt_zeta)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do

    end subroutine












    !> Generate a set of points defining a 2x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_2x1x1_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik), parameter      :: npt = 12
        integer(ik)                 :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr
        real(rk), dimension(npt)    :: x,y,z

        ! elements (3x1x1) - linear
        !
        !      *------*------*
        !      |      |      | 
        !      |      |      | 
        !      *------*------*
        !



        x = [ZERO, ONE, TWO, &       
             ZERO, ONE, TWO, &       
             ZERO, ONE, TWO, &       
             ZERO, ONE, TWO]
             

        y = [ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, &
             ZERO, ZERO, ZERO, &
             ONE, ONE, ONE]

        z = [ZERO, ZERO, ZERO, &
             ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, &
             ONE, ONE, ONE]



        ! Allocate point storage
        allocate(pts(3,2,2), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        do ipt_zeta = 1,2
            do ipt_eta = 1,2
                do ipt_xi = 1,3
                    call pts(ipt_xi,ipt_eta,ipt_zeta)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do

    end subroutine
























    !> Generate a set of points defining a 3x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_3x1x1_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik), parameter      :: npt = 16
        integer(ik)                 :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr
        real(rk), dimension(npt)    :: x,y,z

        ! elements (3x1x1) - linear
        !
        !      *------*------*------*
        !      |      |      |      | 
        !      |      |      |      | 
        !      *------*------*------*
        !



        x = [ZERO, ONE, TWO, THREE, &       
             ZERO, ONE, TWO, THREE, &       
             ZERO, ONE, TWO, THREE, &       
             ZERO, ONE, TWO, THREE]
             

        y = [ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, &
             ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE]

        z = [ZERO, ZERO, ZERO, ZERO, &
             ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, &
             ONE, ONE, ONE, ONE]



        ! Allocate point storage
        allocate(pts(4,2,2), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        do ipt_zeta = 1,2
            do ipt_eta = 1,2
                do ipt_xi = 1,4
                    call pts(ipt_xi,ipt_eta,ipt_zeta)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do

    end subroutine



















    !> Generate a set of points defining a 4x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_40x15x1_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, npts_y, npts_z
        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (40x15x1) - linear

        npts_x = 41
        npts_y = 16
        npts_z = 2

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 1._rk

        ! Allocate point storage
        allocate(pts(npts_x,npts_y,npts_z), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        z = ZERO
        do ipt_zeta = 1,npts_z
            y = ZERO

            do ipt_eta = 1,npts_y
                x = ZERO

                do ipt_xi = 1,npts_x
                    call pts(ipt_xi,ipt_eta,ipt_zeta)%set(x, y, z)
                    ipt = ipt + 1
                    x = x + dx
                end do

                y = y + dy
            end do

            z = z + dz
        end do

    end subroutine













    !> Generate a set of points defining a 4x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_15x15x1_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, npts_y, npts_z
        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (15x15x1) - linear

        npts_x = 16
        npts_y = 16
        npts_z = 2

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 1.0_rk

        ! Allocate point storage
        allocate(pts(npts_x,npts_y,npts_z), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        z = ZERO
        do ipt_zeta = 1,npts_z
            y = ZERO

            do ipt_eta = 1,npts_y
                x = ZERO

                do ipt_xi = 1,npts_x
                    call pts(ipt_xi,ipt_eta,ipt_zeta)%set(x, y, z)
                    ipt = ipt + 1
                    x = x + dx
                end do

                y = y + dy
            end do

            z = z + dz
        end do

    end subroutine












    !> Generate a set of points defining a 15x15x2 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_15x15x2_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, npts_y, npts_z
        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (15x15x2) - linear

        npts_x = 16
        npts_y = 16
        npts_z = 3

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 0.5_rk

        ! Allocate point storage
        allocate(pts(npts_x,npts_y,npts_z), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        z = ZERO
        do ipt_zeta = 1,npts_z
            y = ZERO

            do ipt_eta = 1,npts_y
                x = ZERO

                do ipt_xi = 1,npts_x
                    call pts(ipt_xi,ipt_eta,ipt_zeta)%set(x, y, z)
                    ipt = ipt + 1
                    x = x + dx
                end do

                y = y + dy
            end do

            z = z + dz
        end do

    end subroutine



















    !> Generate a set of points defining a 15x15x3 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_15x15x3_linear(pts)
        type(point_t), allocatable, intent(inout)  :: pts(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, npts_y, npts_z
        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (15x15x3) - linear

        npts_x = 16
        npts_y = 16
        npts_z = 4

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 0.5_rk

        ! Allocate point storage
        allocate(pts(npts_x,npts_y,npts_z), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        z = ZERO
        do ipt_zeta = 1,npts_z
            y = ZERO

            do ipt_eta = 1,npts_y
                x = ZERO

                do ipt_xi = 1,npts_x
                    call pts(ipt_xi,ipt_eta,ipt_zeta)%set(x, y, z)
                    ipt = ipt + 1
                    x = x + dx
                end do

                y = y + dy
            end do

            z = z + dz
        end do

    end subroutine








































end module mod_testutils
