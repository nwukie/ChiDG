module mod_grid
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE, TWO, ZERO, TWO_DIM, THREE_DIM, NFACES, &
                                  XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, NO_PROC, NO_DIFF
    use mod_polynomial,     only: polynomial_val
    use type_densematrix,   only: densematrix_t
    use mod_inv
    implicit none

    ! coordinate mapping matrices
    integer(ik),         parameter      :: NORDER = 4
    type(densematrix_t), save,  target  :: ELEM_MAP_2D(NORDER)  !< array of matrices
    type(densematrix_t), save,  target  :: ELEM_MAP_3D(NORDER)  !< array of matrices

    integer(ik),         parameter      :: NFACE_CORNERS = 4
    integer(ik),         save           :: FACE_CORNERS(NFACES,NFACE_CORNERS,NORDER)    !< array of indices for element corners

    integer(ik),         save           :: ELEMENT_VERTEX_INDICES(8, NORDER)


    logical :: uninitialized = .true.

contains



    !>  Call grid initialization routines. This is called by chidg%init('env'). So, should
    !!  not need to call this explicity.
    !!
    !!      - calls computes_element_mappings
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine initialize_grid()

        if (uninitialized) then
            call compute_element_mappings(TWO_DIM)
            call compute_element_mappings(THREE_DIM)
            call compute_face_corner_indices()
            call find_element_vertex_indices(THREE_DIM)
        end if

        uninitialized = .false.

    end subroutine initialize_grid
    !****************************************************************************************************





    !>  Compute matrix to convert element discrete coordinates to modal coordinates. Initializes the
    !!  array of denseblock matrices in ELEM_MAP.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  TODO: TEST
    !!  TODO: Generalize better for spatial dimension
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine compute_element_mappings(spacedim)
        use type_point, only: point_t
        integer(ik),    intent(in)  :: spacedim

        type(point_t),  allocatable :: nodes(:)
        real(rk),       allocatable :: xi(:),eta(:),zeta(:)
        integer(ik)                 :: npts_1d(7), npts_2d(7), npts_3d(7)
        integer(ik)                 :: ierr, iorder, iterm, inode, ipt
        integer(ik)                 :: ixi,  ieta, izeta

        ! Mapping order
        ! [linear, quadratic, cubic, quartic]
        npts_1d = [2,3,4,5,6,7,8]           ! Number of points defining an edge
        npts_2d = npts_1d*npts_1d           ! Number of points defining an element in 2D
        npts_3d = npts_1d*npts_1d*npts_1d   ! Number of points defining an element in 3D


        ! Loop through and compute mapping for different element types
        do iorder = 1,NORDER

            ! Initialize mapping for reference element.
            if ( spacedim == THREE_DIM ) then
                call ELEM_MAP_3D(iorder)%init(npts_3d(iorder),1,0,0,0,0,NO_PROC,0,0,NO_DIFF)
            else if ( spacedim == TWO_DIM ) then
                call ELEM_MAP_2D(iorder)%init(npts_2d(iorder),1,0,0,0,0,NO_PROC,0,0,NO_DIFF)
            end if


            ! Allocate storage for nodes and coordinates.
            if ( spacedim == THREE_DIM ) then
                allocate(nodes(npts_3d(iorder)),  &
                         xi(npts_1d(iorder)),     &
                         eta(npts_1d(iorder)),    &
                         zeta(npts_1d(iorder)), stat=ierr)
                if (ierr /= 0) call AllocationError

            else if ( spacedim == TWO_DIM ) then
                allocate(nodes(npts_2d(iorder)),  &
                         xi(npts_1d(iorder)),     &
                         eta(npts_1d(iorder)),    &
                         zeta(npts_1d(iorder)), stat=ierr)
                if (ierr /= 0) call AllocationError

            end if


            ! Compute 1d point coordinates in each direction
            do ipt = 1,npts_1d(iorder)
                xi(ipt)   = -ONE + (real(ipt,rk)-ONE)*(TWO/(real(npts_1d(iorder),rk)-ONE))
                eta(ipt)  = -ONE + (real(ipt,rk)-ONE)*(TWO/(real(npts_1d(iorder),rk)-ONE))
                zeta(ipt) = -ONE + (real(ipt,rk)-ONE)*(TWO/(real(npts_1d(iorder),rk)-ONE))
            end do


            ! Set up reference mesh nodes in each direction
            inode = 1
            if ( spacedim == THREE_DIM ) then

                do izeta = 1,npts_1d(iorder)
                    do ieta = 1,npts_1d(iorder)
                        do ixi = 1,npts_1d(iorder)
                            call nodes(inode)%set(xi(ixi), eta(ieta), zeta(izeta))
                            inode = inode + 1
                        end do
                    end do
                end do

            else if ( spacedim == TWO_DIM ) then

                do ieta = 1,npts_1d(iorder)
                    do ixi = 1,npts_1d(iorder)
                        call nodes(inode)%set(xi(ixi), eta(ieta), ZERO)
                        inode = inode + 1
                    end do
                end do

            else
                call chidg_signal(FATAL,"mod_grid::compute_elements_mappings - Invalid spacedim")

            end if


            ! Compute the values of each mapping term at each mesh point
            if ( spacedim == THREE_DIM ) then
                do iterm = 1,npts_3d(iorder)
                    do inode = 1,npts_3d(iorder)
                        ELEM_MAP_3D(iorder)%mat(inode,iterm) = polynomial_val(3,npts_3d(iorder),iterm,[nodes(inode)%c1_, nodes(inode)%c2_, nodes(inode)%c3_])
                    end do
                end do

            else if ( spacedim == TWO_DIM ) then
                do iterm = 1,npts_2d(iorder)
                    do inode = 1,npts_2d(iorder)
                        ELEM_MAP_2D(iorder)%mat(inode,iterm) = polynomial_val(2,npts_2d(iorder),iterm,[nodes(inode)%c1_, nodes(inode)%c2_, nodes(inode)%c3_])
                    end do
                end do

            end if 

            
            ! Invert matrix so that it can multiply a vector of
            ! element points to compute the mode amplitudes of the x,y mappings
            if ( spacedim == THREE_DIM ) then
                ELEM_MAP_3D(iorder)%mat = inv(ELEM_MAP_3D(iorder)%mat)
            else if ( spacedim == TWO_DIM ) then
                ELEM_MAP_2D(iorder)%mat = inv(ELEM_MAP_2D(iorder)%mat)
            end if


            ! Dellocate variables for next iteration in loop
            deallocate(nodes, xi, eta, zeta)

        end do


    end subroutine compute_element_mappings
    !**********************************************************************************************************




    !>  Return a matrix to compute the coordinate expansion.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !----------------------------------------------------------------------------------------------------------
    function get_element_mapping(spacedim,iorder) result(matrix)
        integer(ik),    intent(in) :: spacedim
        integer(ik),    intent(in) :: iorder

        real(rk),   allocatable :: matrix(:,:)

        if ( allocated(ELEM_MAP_2D(iorder)%mat) .and. allocated(ELEM_MAP_3D(iorder)%mat) ) then


            if ( spacedim == TWO_DIM ) then
                matrix = ELEM_MAP_2D(iorder)%mat
            else if ( spacedim == THREE_DIM ) then
                matrix = ELEM_MAP_3D(iorder)%mat
            end if

        else
            call chidg_signal(FATAL,"get_element_mapping: selected element mapping is not allocated. Probably need to call chidg%init('env'). Or, maybe your mapping is incorrect.")
        end if

    end function get_element_mapping
    !**********************************************************************************************************









    !> Compute the indices of the nodes that are the face corners
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/23/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------------
    subroutine compute_face_corner_indices()
        integer(ik) :: iface, icorner, iorder, base_corners(4)
        integer(ik) :: corner_one, corner_two, corner_three, corner_four, &
                       corner_five, corner_six, corner_seven, corner_eight

        do iface = 1,NFACES
            do iorder = 1,NORDER

                ! Still Plot3D faces, but CGNS nodes
                select case ( iface )
                    case ( XI_MIN )
                        face_corners(iface,1,iorder) = 1
                        face_corners(iface,2,iorder) = 5
                        face_corners(iface,3,iorder) = 8
                        face_corners(iface,4,iorder) = 4
                    case ( XI_MAX )
                        face_corners(iface,1,iorder) = 2
                        face_corners(iface,2,iorder) = 3
                        face_corners(iface,3,iorder) = 7
                        face_corners(iface,4,iorder) = 6
                    case ( ETA_MIN )
                        face_corners(iface,1,iorder) = 1
                        face_corners(iface,2,iorder) = 2
                        face_corners(iface,3,iorder) = 6
                        face_corners(iface,4,iorder) = 5
                    case ( ETA_MAX )
                        face_corners(iface,1,iorder) = 3
                        face_corners(iface,2,iorder) = 4
                        face_corners(iface,3,iorder) = 8
                        face_corners(iface,4,iorder) = 7
                    case ( ZETA_MIN )
                        face_corners(iface,1,iorder) = 1
                        face_corners(iface,2,iorder) = 4
                        face_corners(iface,3,iorder) = 3
                        face_corners(iface,4,iorder) = 2
                    case ( ZETA_MAX )
                        face_corners(iface,1,iorder) = 5
                        face_corners(iface,2,iorder) = 6
                        face_corners(iface,3,iorder) = 7
                        face_corners(iface,4,iorder) = 8
                    case default
                        call chidg_signal(FATAL,"compute_face_corner_indices")
                end select

            end do !iorder
        end do !iface

    end subroutine compute_face_corner_indices
    !**********************************************************************************************************




    !>
    !!
    !!  @author Eric Wolf
    !!
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine find_element_vertex_indices(spacedim) 
        integer(ik),    intent(in)  :: spacedim

        integer(ik) :: ierr, iorder

        ! Loop through and compute mapping for different element types
        do iorder = 1,NORDER

            ! Initialize element vertex indices
            select case (spacedim)
                case (THREE_DIM)

                    element_vertex_indices(1,iorder) = 1
                    element_vertex_indices(2,iorder) = 2
                    element_vertex_indices(3,iorder) = 3
                    element_vertex_indices(4,iorder) = 4
                    element_vertex_indices(5,iorder) = 5
                    element_vertex_indices(6,iorder) = 6
                    element_vertex_indices(7,iorder) = 7
                    element_vertex_indices(8,iorder) = 8

                case default
                    call chidg_signal(FATAL,"mod_grid::find_element_vertex_indices - Invalid spacedim")

            end select

        end do !iorder


    end subroutine find_element_vertex_indices
    !**********************************************************************************************************








end module mod_grid
