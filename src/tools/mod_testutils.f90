module mod_testutils
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX
    use mod_chidg_mpi,              only: IRANK
    use mod_plot3d_utilities,       only: get_block_points_plot3d, &
                                          get_block_elements_plot3d, &
                                          get_block_boundary_faces_plot3d
    use mod_hdf_utilities,          only: initialize_file_hdf, add_domain_hdf, &
                                          open_domain_hdf, close_domain_hdf, &
                                          set_bc_patch_hdf, add_bc_state_hdf, &
                                          set_contains_grid_hdf, close_file_hdf
    use mod_bc,                     only: create_bc
    use type_point,                 only: point_t
    use type_bc_state,              only: bc_state_t
    use type_domain_connectivity,   only: domain_connectivity_t
    use hdf5
    implicit none


contains



    !>  Create an actual ChiDG-formatted grid file that could be
    !!  read in by a test. Also with initialized boundary conditions.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !---------------------------------------------------------------------------
    subroutine create_mesh_file(selector, filename)
        character(*),   intent(in)  :: selector
        character(*),   intent(in)  :: filename

        integer(ik) :: ierr

        ! Open HDF5.
        call h5open_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"create_mesh_file: h5open_f")



        ! Generate grid file base on selector case.
        select case (trim(selector))
            case("D2_E8_M1 : Overlapping : Matching")
                call create_mesh_file_D2E8M1_overlapping_matching(filename)

            case default
                call chidg_signal(FATAL,"create_mesh_file: There was no valid case that matched the incoming string")

        end select


        ! Close HDF5.
        call h5close_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"create_mesh_file: h5close_f")


    end subroutine create_mesh_file
    !***************************************************************************







    !>  Generate a set of points for a mesh. String input calls specialized
    !!  procedure for generating the points
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/24/2016
    !!
    !!  @param[in]      string          Character string used to select a specialized 
    !!                                  meshgen call
    !!  @param[inout]   nodes           Array of node coordinates for the grid
    !!  @param[inout]   connectivity    Connectivity data for the grid
    !--------------------------------------------------------------------
    subroutine create_mesh(string,nodes,connectivity)
        character(*),                   intent(in)      :: string
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity


        select case (trim(string))
            case ('1x1x1','111')
                call meshgen_1x1x1_linear(nodes,connectivity)

            case ('1x1x1_unit','111u')
                call meshgen_1x1x1_unit_linear(nodes,connectivity)

            case ('3x3x3','333')
                call meshgen_3x3x3_linear(nodes,connectivity)

            case ('3x3x3_unit','333u')
                call meshgen_3x3x3_unit_linear(nodes,connectivity)

            case ('2x2x2','222')
                call meshgen_2x2x2_linear(nodes,connectivity)

            case ('2x2x1','221')
                call meshgen_2x2x1_linear(nodes,connectivity)

            case ('3x3x1','331')
                call meshgen_3x3x1_linear(nodes,connectivity)

            case ('4x1x1','411')
                call meshgen_4x1x1_linear(nodes,connectivity)

            case ('3x1x1','311')
                call meshgen_3x1x1_linear(nodes,connectivity)

            case ('2x1x1','211')
                call meshgen_2x1x1_linear(nodes,connectivity)

            case ('40x15x1')
                call meshgen_40x15x1_linear(nodes,connectivity)

            case ('15x15x1')
                call meshgen_15x15x1_linear(nodes,connectivity)

            case ('15x15x2')
                call meshgen_15x15x2_linear(nodes,connectivity)

            case ('15x15x3')
                call meshgen_15x15x3_linear(nodes,connectivity)


            case default
                call chidg_signal(FATAL,'String identifying mesh generation routine was not recognized')
        end select


        !
        ! Generate nodes
        !
!        mapping = 1
!        idomain = 1
!        nodes    = get_block_points_plot3d(xcoords,ycoords,zcoords)
!        elements = get_block_elements_plot3d(xcoords,ycoords,zcoords,mapping,idomain)


    end subroutine create_mesh
    !****************************************************************************











    !> Generate a set of points defining a 1x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/24/2016
    !!
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_1x1x1_linear(nodes,connectivity)
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, &
                       npts_y, npts_z, npts, nelem, npts_element
        integer(ik) :: ielem_x, ielem_y, ielem_z, ipt_x, ipt_y, ipt_z, &
                       ipt_elem, ielem, xstart, ystart, zstart
        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (1x1x1) - linear

        npts_x = 2
        npts_y = 2
        npts_z = 2

        dx = 1._rk
        dy = 1._rk
        dz = 1._rk

        ! Allocate point storage
        npts = npts_x*npts_y*npts_z
        nelem = (npts_x-1)*(npts_y-1)*(npts_z-1)
        npts_element = 8    !linear
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError

        call connectivity%init(nelem,npts)

        ipt = 1
        z = ZERO
        do ipt_zeta = 1,npts_z
            y = ZERO
            do ipt_eta = 1,npts_y
                x = ZERO
                do ipt_xi = 1,npts_x


                    call nodes(ipt)%set(x, y, z)
                    ipt = ipt + 1


                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do

        ielem = 1
        do ielem_z = 1,npts_z-1
            do ielem_y = 1,npts_y-1
                do ielem_x = 1,npts_x-1
                    call connectivity%data(ielem)%init(1)
                    connectivity%data(ielem)%data(1) = 1        ! idomain
                    connectivity%data(ielem)%data(2) = ielem    ! ielement
                    connectivity%data(ielem)%data(3) = 1        ! mapping
                    call connectivity%data(ielem)%set_element_partition(IRANK)

                    xstart = 1 + (ielem_x-1)
                    ystart = 1 + (ielem_y-1)
                    zstart = 1 + (ielem_z-1)
                    
                    
                    ipt_elem = 1
                    do ipt_z = zstart,(zstart+1)
                        do ipt_y = ystart,(ystart+1)
                            do ipt_x = xstart,(xstart+1)
                                ipt = ipt_x  +  (ipt_y-1)*npts_x  +  (ipt_z-1)*(npts_x*npts_y)     
                                !connectivity(ielem,2+ipt_elem) = ipt
                                call connectivity%data(ielem)%set_element_node(ipt_elem, ipt)
                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do

                    ielem = ielem + 1
                end do
            end do
        end do





    end subroutine meshgen_1x1x1_linear
    !**************************************************************************************














    !> Generate a set of points defining a 1x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_1x1x1_unit_linear(nodes,connectivity)
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, &
                       npts_y, npts_z, npts, nelem, npts_element
        integer(ik) :: ielem_x, ielem_y, ielem_z, ipt_x, ipt_y, ipt_z, &
                       ipt_elem, ielem, xstart, ystart, zstart
        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (1x1x1) - linear

        npts_x = 2
        npts_y = 2
        npts_z = 2

        dx = 2._rk
        dy = 2._rk
        dz = 2._rk

        ! Allocate point storage
        npts = npts_x*npts_y*npts_z
        nelem = (npts_x-1)*(npts_y-1)*(npts_z-1)
        npts_element = 8    !linear

        !allocate(nodes(npts), connectivity(nelem,2+npts_element), stat=ierr)
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError

        call connectivity%init(nelem,npts)




        ipt = 1
        z = -ONE
        do ipt_zeta = 1,npts_z
            y = -ONE
            do ipt_eta = 1,npts_y
                x = -ONE
                do ipt_xi = 1,npts_x

                    call nodes(ipt)%set(x, y, z)
                    ipt = ipt + 1

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do



        ielem = 1
        do ielem_z = 1,npts_z-1
            do ielem_y = 1,npts_y-1
                do ielem_x = 1,npts_x-1
                    call connectivity%data(ielem)%init(1)
                    connectivity%data(ielem)%data(1) = 1        ! idomain
                    connectivity%data(ielem)%data(2) = ielem    ! ielement
                    connectivity%data(ielem)%data(3) = 1        ! mapping
                    call connectivity%data(ielem)%set_element_partition(IRANK)

                    xstart = 1 + (ielem_x-1)
                    ystart = 1 + (ielem_y-1)
                    zstart = 1 + (ielem_z-1)
                    
                    
                    ipt_elem = 1
                    do ipt_z = zstart,(zstart+1)
                        do ipt_y = ystart,(ystart+1)
                            do ipt_x = xstart,(xstart+1)
                                ipt = ipt_x  +  (ipt_y-1)*npts_x  +  (ipt_z-1)*(npts_x*npts_y)     
                                !connectivity(ielem,2+ipt_elem) = ipt
                                call connectivity%data(ielem)%set_element_node(ipt_elem, ipt)
                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do

                    ielem = ielem + 1
                end do
            end do
        end do







    end subroutine meshgen_1x1x1_unit_linear
    !**************************************************************************************












    !> Generate a set of points defining a 2x2x2 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_2x2x2_linear(nodes,connectivity)
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, &
                       npts_y, npts_z, npts, nelem, npts_element
        integer(ik) :: ielem_x, ielem_y, ielem_z, ipt_x, ipt_y, ipt_z, &
                       ipt_elem, ielem, xstart, ystart, zstart

        integer(ik), parameter      :: npt = 27
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
        npts_x = 3
        npts_y = 3
        npts_z = 3
        npts = npts_x*npts_y*npts_z
        nelem = (npts_x-1)*(npts_y-1)*(npts_z-1)
        npts_element = 8    !linear

        !allocate(nodes(npts), connectivity(nelem,2+npts_element), stat=ierr)
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError

        call connectivity%init(nelem,npts)





        ipt = 1
        do ipt_zeta = 1,3
            do ipt_eta = 1,3
                do ipt_xi = 1,3
                    call nodes(ipt)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do




        ielem = 1
        do ielem_z = 1,npts_z-1
            do ielem_y = 1,npts_y-1
                do ielem_x = 1,npts_x-1
                    call connectivity%data(ielem)%init(1)
                    connectivity%data(ielem)%data(1) = 1        ! idomain
                    connectivity%data(ielem)%data(2) = ielem    ! ielement
                    connectivity%data(ielem)%data(3) = 1        ! mapping
                    call connectivity%data(ielem)%set_element_partition(IRANK)

                    xstart = 1 + (ielem_x-1)
                    ystart = 1 + (ielem_y-1)
                    zstart = 1 + (ielem_z-1)
                    
                    
                    ipt_elem = 1
                    do ipt_z = zstart,(zstart+1)
                        do ipt_y = ystart,(ystart+1)
                            do ipt_x = xstart,(xstart+1)
                                ipt = ipt_x  +  (ipt_y-1)*npts_x  +  (ipt_z-1)*(npts_x*npts_y)     
                                !connectivity(ielem,2+ipt_elem) = ipt
                                call connectivity%data(ielem)%set_element_node(ipt_elem, ipt)
                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do

                    ielem = ielem + 1
                end do
            end do
        end do





    end subroutine meshgen_2x2x2_linear
    !**************************************************************************************













    !> Generate a set of points defining a 2x2x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_2x2x1_linear(nodes,connectivity)
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, &
                       npts_y, npts_z, npts, nelem, npts_element
        integer(ik) :: ielem_x, ielem_y, ielem_z, ipt_x, ipt_y, ipt_z, &
                       ipt_elem, ielem, xstart, ystart, zstart

        integer(ik), parameter      :: npt = 18
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
        npts_x = 3
        npts_y = 3
        npts_z = 2
        npts = npts_x*npts_y*npts_z
        nelem = (npts_x-1)*(npts_y-1)*(npts_z-1)
        npts_element = 8    !linear

        !allocate(nodes(npts), connectivity(nelem,2+npts_element), stat=ierr)
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError

        call connectivity%init(nelem,npts)




        ipt = 1
        do ipt_zeta = 1,2
            do ipt_eta = 1,3
                do ipt_xi = 1,3
                    call nodes(ipt)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do




        ielem = 1
        do ielem_z = 1,npts_z-1
            do ielem_y = 1,npts_y-1
                do ielem_x = 1,npts_x-1
                    call connectivity%data(ielem)%init(1)
                    connectivity%data(ielem)%data(1) = 1        ! idomain
                    connectivity%data(ielem)%data(2) = ielem    ! ielement
                    connectivity%data(ielem)%data(3) = 1        ! mapping
                    call connectivity%data(ielem)%set_element_partition(IRANK)

                    xstart = 1 + (ielem_x-1)
                    ystart = 1 + (ielem_y-1)
                    zstart = 1 + (ielem_z-1)
                    
                    
                    ipt_elem = 1
                    do ipt_z = zstart,(zstart+1)
                        do ipt_y = ystart,(ystart+1)
                            do ipt_x = xstart,(xstart+1)
                                ipt = ipt_x  +  (ipt_y-1)*npts_x  +  (ipt_z-1)*(npts_x*npts_y)     
                                !connectivity(ielem,2+ipt_elem) = ipt
                                call connectivity%data(ielem)%set_element_node(ipt_elem, ipt)
                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do

                    ielem = ielem + 1
                end do
            end do
        end do





    end subroutine meshgen_2x2x1_linear
    !***************************************************************************************















    !> Generate a set of points defining a 3x3x3 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_3x3x3_linear(nodes,connectivity)
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, &
                       npts_y, npts_z, npts, nelem, npts_element
        integer(ik) :: ielem_x, ielem_y, ielem_z, ipt_x, ipt_y, ipt_z, &
                       ipt_elem, ielem, xstart, ystart, zstart

        integer(ik), parameter      :: npt = 64
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
        x = [ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE]

        y = [ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
             THREE, THREE, THREE, THREE, &
             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
             THREE, THREE, THREE, THREE, &
             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
             THREE, THREE, THREE, THREE, &
             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
             THREE, THREE, THREE, THREE]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
             ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, & 
             ONE, ONE, ONE, &
             TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, &
             TWO, TWO, TWO, &
             THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, &
             THREE, THREE, THREE, THREE, THREE, THREE]


        ! Allocate point storage
        npts_x = 4
        npts_y = 4
        npts_z = 4
        npts = npts_x*npts_y*npts_z
        nelem = (npts_x-1)*(npts_y-1)*(npts_z-1)
        npts_element = 8    !linear


        !allocate(nodes(npts), connectivity(nelem,2+npts_element), stat=ierr)
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError

        call connectivity%init(nelem,npts)


        ipt = 1
        do ipt_zeta = 1,4
            do ipt_eta = 1,4
                do ipt_xi = 1,4
                    call nodes(ipt)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do


        ielem = 1
        do ielem_z = 1,npts_z-1
            do ielem_y = 1,npts_y-1
                do ielem_x = 1,npts_x-1
                    call connectivity%data(ielem)%init(1)
                    connectivity%data(ielem)%data(1) = 1        ! idomain
                    connectivity%data(ielem)%data(2) = ielem    ! ielement
                    connectivity%data(ielem)%data(3) = 1        ! mapping
                    call connectivity%data(ielem)%set_element_partition(IRANK)

                    xstart = 1 + (ielem_x-1)
                    ystart = 1 + (ielem_y-1)
                    zstart = 1 + (ielem_z-1)
                    
                    
                    ipt_elem = 1
                    do ipt_z = zstart,(zstart+1)
                        do ipt_y = ystart,(ystart+1)
                            do ipt_x = xstart,(xstart+1)
                                ipt = ipt_x  +  (ipt_y-1)*npts_x  +  (ipt_z-1)*(npts_x*npts_y)     
                                !connectivity(ielem,2+ipt_elem) = ipt
                                call connectivity%data(ielem)%set_element_node(ipt_elem, ipt)
                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do

                    ielem = ielem + 1
                end do
            end do
        end do



    end subroutine meshgen_3x3x3_linear
    !***************************************************************************************
















    !> Generate a set of points defining a 3x3x3 unit-element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_3x3x3_unit_linear(nodes,connectivity)
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, npts_y, &
                       npts_z, npts, nelem, npts_element
        integer(ik) :: ielem_x, ielem_y, ielem_z, ipt_x, ipt_y, ipt_z, ipt_elem, &
                       ielem, xstart, ystart, zstart

        integer(ik), parameter      :: npt = 64
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
        npts_x = 4
        npts_y = 4
        npts_z = 4
        npts = npts_x*npts_y*npts_z
        nelem = (npts_x-1)*(npts_y-1)*(npts_z-1)
        npts_element = 8    !linear

        !allocate(nodes(npts), connectivity(nelem,2+npts_element), stat=ierr)
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError

        call connectivity%init(nelem,npts)




        ipt = 1
        do ipt_zeta = 1,4
            do ipt_eta = 1,4
                do ipt_xi = 1,4
                    call nodes(ipt)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do




        ielem = 1
        do ielem_z = 1,npts_z-1
            do ielem_y = 1,npts_y-1
                do ielem_x = 1,npts_x-1
                    call connectivity%data(ielem)%init(1)
                    connectivity%data(ielem)%data(1) = 1        ! idomain
                    connectivity%data(ielem)%data(2) = ielem    ! ielement
                    connectivity%data(ielem)%data(3) = 1        ! mapping
                    call connectivity%data(ielem)%set_element_partition(IRANK)

                    xstart = 1 + (ielem_x-1)
                    ystart = 1 + (ielem_y-1)
                    zstart = 1 + (ielem_z-1)
                    
                    
                    ipt_elem = 1
                    do ipt_z = zstart,(zstart+1)
                        do ipt_y = ystart,(ystart+1)
                            do ipt_x = xstart,(xstart+1)
                                ipt = ipt_x  +  (ipt_y-1)*npts_x  +  (ipt_z-1)*(npts_x*npts_y)     
                                !connectivity(ielem,2+ipt_elem) = ipt
                                call connectivity%data(ielem)%set_element_node(ipt_elem, ipt)
                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do

                    ielem = ielem + 1
                end do
            end do
        end do





    end subroutine meshgen_3x3x3_unit_linear
    !***************************************************************************************












    !> Generate a set of points defining a 3x3x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_3x3x1_linear(nodes,connectivity)
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, &
                       npts_y, npts_z, npts, nelem, npts_element
        integer(ik) :: ielem_x, ielem_y, ielem_z, ipt_x, ipt_y, ipt_z, &
                       ipt_elem, ielem, xstart, ystart, zstart

        integer(ik), parameter      :: npt = 32
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
        npts_x = 4
        npts_y = 4
        npts_z = 2
        npts = npts_x*npts_y*npts_z
        nelem = (npts_x-1)*(npts_y-1)*(npts_z-1)
        npts_element = 8    !linear

        !allocate(nodes(npts), connectivity(nelem,2+npts_element), stat=ierr)
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError

        call connectivity%init(nelem,npts)



        ipt = 1
        do ipt_zeta = 1,2
            do ipt_eta = 1,4
                do ipt_xi = 1,4
                    call nodes(ipt)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do


        ielem = 1
        do ielem_z = 1,npts_z-1
            do ielem_y = 1,npts_y-1
                do ielem_x = 1,npts_x-1
                    call connectivity%data(ielem)%init(1)
                    connectivity%data(ielem)%data(1) = 1        ! idomain
                    connectivity%data(ielem)%data(2) = ielem    ! ielement
                    connectivity%data(ielem)%data(3) = 1        ! mapping
                    call connectivity%data(ielem)%set_element_partition(IRANK)

                    xstart = 1 + (ielem_x-1)
                    ystart = 1 + (ielem_y-1)
                    zstart = 1 + (ielem_z-1)
                    
                    
                    ipt_elem = 1
                    do ipt_z = zstart,(zstart+1)
                        do ipt_y = ystart,(ystart+1)
                            do ipt_x = xstart,(xstart+1)
                                ipt = ipt_x  +  (ipt_y-1)*npts_x  +  (ipt_z-1)*(npts_x*npts_y)     
                                !connectivity(ielem,2+ipt_elem) = ipt
                                call connectivity%data(ielem)%set_element_node(ipt_elem, ipt)
                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do

                    ielem = ielem + 1
                end do
            end do
        end do



    end subroutine meshgen_3x3x1_linear
    !***************************************************************************************
















    !> Generate a set of points defining a 4x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_4x1x1_linear(nodes,connectivity)
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, npts_y, npts_z, npts, nelem, npts_element
        integer(ik) :: ielem_x, ielem_y, ielem_z, ipt_x, ipt_y, ipt_z, ipt_elem, ielem, xstart, ystart, zstart

        integer(ik), parameter      :: npt = 20
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
        npts_x = 5
        npts_y = 2
        npts_z = 2
        npts = npts_x*npts_y*npts_z
        nelem = (npts_x-1)*(npts_y-1)*(npts_z-1)
        npts_element = 8    !linear

        !allocate(nodes(npts), connectivity(nelem,2+npts_element), stat=ierr)
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError
        call connectivity%init(nelem,npts)



        ipt = 1
        do ipt_zeta = 1,2
            do ipt_eta = 1,2
                do ipt_xi = 1,5
                    call nodes(ipt)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do



        ielem = 1
        do ielem_z = 1,npts_z-1
            do ielem_y = 1,npts_y-1
                do ielem_x = 1,npts_x-1
                    call connectivity%data(ielem)%init(1)
                    connectivity%data(ielem)%data(1) = 1        ! idomain
                    connectivity%data(ielem)%data(2) = ielem    ! ielement
                    connectivity%data(ielem)%data(3) = 1        ! mapping
                    call connectivity%data(ielem)%set_element_partition(IRANK)

                    xstart = 1 + (ielem_x-1)
                    ystart = 1 + (ielem_y-1)
                    zstart = 1 + (ielem_z-1)
                    
                    
                    ipt_elem = 1
                    do ipt_z = zstart,(zstart+1)
                        do ipt_y = ystart,(ystart+1)
                            do ipt_x = xstart,(xstart+1)
                                ipt = ipt_x  +  (ipt_y-1)*npts_x  +  (ipt_z-1)*(npts_x*npts_y)     
                                !connectivity(ielem,2+ipt_elem) = ipt
                                call connectivity%data(ielem)%set_element_node(ipt_elem, ipt)
                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do

                    ielem = ielem + 1
                end do
            end do
        end do



    end subroutine meshgen_4x1x1_linear
    !***************************************************************************************












    !> Generate a set of points defining a 2x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_2x1x1_linear(nodes,connectivity)
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, &
                       npts_y, npts_z, npts, nelem, npts_element
        integer(ik) :: ielem_x, ielem_y, ielem_z, ipt_x, ipt_y, ipt_z, &
                       ipt_elem, ielem, xstart, ystart, zstart

        integer(ik), parameter      :: npt = 12
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
        npts_x = 3
        npts_y = 2
        npts_z = 2
        npts = npts_x*npts_y*npts_z
        nelem = (npts_x-1)*(npts_y-1)*(npts_z-1)
        npts_element = 8    !linear

        !allocate(nodes(npts), connectivity(nelem,2+npts_element), stat=ierr)
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError
        call connectivity%init(nelem,npts)




        ipt = 1
        do ipt_zeta = 1,2
            do ipt_eta = 1,2
                do ipt_xi = 1,3
                    call nodes(ipt)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do


        ielem = 1
        do ielem_z = 1,npts_z-1
            do ielem_y = 1,npts_y-1
                do ielem_x = 1,npts_x-1
                    call connectivity%data(ielem)%init(1)
                    connectivity%data(ielem)%data(1) = 1        ! idomain
                    connectivity%data(ielem)%data(2) = ielem    ! ielement
                    connectivity%data(ielem)%data(3) = 1        ! mapping
                    call connectivity%data(ielem)%set_element_partition(IRANK)

                    xstart = 1 + (ielem_x-1)
                    ystart = 1 + (ielem_y-1)
                    zstart = 1 + (ielem_z-1)
                    
                    
                    ipt_elem = 1
                    do ipt_z = zstart,(zstart+1)
                        do ipt_y = ystart,(ystart+1)
                            do ipt_x = xstart,(xstart+1)
                                ipt = ipt_x  +  (ipt_y-1)*npts_x  +  (ipt_z-1)*(npts_x*npts_y)     
                                !connectivity(ielem,2+ipt_elem) = ipt
                                call connectivity%data(ielem)%set_element_node(ipt_elem, ipt)
                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do

                    ielem = ielem + 1
                end do
            end do
        end do



    end subroutine meshgen_2x1x1_linear
    !**************************************************************************************
























    !> Generate a set of points defining a 3x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_3x1x1_linear(nodes,connectivity)
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, npts_y, &
                       npts_z, npts, nelem, npts_element
        integer(ik) :: ielem_x, ielem_y, ielem_z, ipt_x, ipt_y, ipt_z, ipt_elem, &
                       ielem, xstart, ystart, zstart

        integer(ik), parameter      :: npt = 16
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
        npts_x = 4
        npts_y = 2
        npts_z = 2
        npts = npts_x*npts_y*npts_z
        nelem = (npts_x-1)*(npts_y-1)*(npts_z-1)
        npts_element = 8    !linear
        !allocate(nodes(npts), connectivity(nelem,2+npts_element), stat=ierr)
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError
        call connectivity%init(nelem,npts)

        ipt = 1
        do ipt_zeta = 1,2
            do ipt_eta = 1,2
                do ipt_xi = 1,4
                    call nodes(ipt)%set(x(ipt), y(ipt), z(ipt))
                    ipt = ipt + 1
                end do
            end do
        end do


        ielem = 1
        do ielem_z = 1,npts_z-1
            do ielem_y = 1,npts_y-1
                do ielem_x = 1,npts_x-1
                    call connectivity%data(ielem)%init(1)
                    connectivity%data(ielem)%data(1) = 1        ! idomain
                    connectivity%data(ielem)%data(2) = ielem    ! ielement
                    connectivity%data(ielem)%data(3) = 1        ! mapping
                    call connectivity%data(ielem)%set_element_partition(IRANK)

                    xstart = 1 + (ielem_x-1)
                    ystart = 1 + (ielem_y-1)
                    zstart = 1 + (ielem_z-1)
                    
                    
                    ipt_elem = 1
                    do ipt_z = zstart,(zstart+1)
                        do ipt_y = ystart,(ystart+1)
                            do ipt_x = xstart,(xstart+1)
                                ipt = ipt_x  +  (ipt_y-1)*npts_x  +  (ipt_z-1)*(npts_x*npts_y)     
                                !connectivity(ielem,2+ipt_elem) = ipt
                                call connectivity%data(ielem)%set_element_node(ipt_elem, ipt)
                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do

                    ielem = ielem + 1
                end do
            end do
        end do




    end subroutine meshgen_3x1x1_linear
    !**************************************************************************************



















    !> Generate a set of points defining a 4x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_40x15x1_linear(nodes,connectivity)
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, npts_y, &
                       npts_z, npts, nelem, npts_element
        integer(ik) :: ielem_x, ielem_y, ielem_z, ipt_x, ipt_y, ipt_z, ipt_elem, &
                       ielem, xstart, ystart, zstart

        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (40x15x1) - linear

        npts_x = 41
        npts_y = 16
        npts_z = 2

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 1._rk

        ! Allocate point storage
        npts  = npts_x*npts_y*npts_z
        nelem = (npts_x-1)*(npts_y-1)*(npts_z-1)
        npts_element = 8    !linear
        !allocate(nodes(npts), connectivity(nelem,2+npts_element), stat=ierr)
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError
        call connectivity%init(nelem,npts)

        ipt = 1
        z = ZERO
        do ipt_zeta = 1,npts_z
            y = ZERO
            do ipt_eta = 1,npts_y
                x = ZERO
                do ipt_xi = 1,npts_x

                    call nodes(ipt)%set(x, y, z)
                    ipt = ipt + 1

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do



        ielem = 1
        do ielem_z = 1,npts_z-1
            do ielem_y = 1,npts_y-1
                do ielem_x = 1,npts_x-1
                    call connectivity%data(ielem)%init(1)
                    connectivity%data(ielem)%data(1) = 1        ! idomain
                    connectivity%data(ielem)%data(2) = ielem    ! ielement
                    connectivity%data(ielem)%data(3) = 1        ! mapping
                    call connectivity%data(ielem)%set_element_partition(IRANK)

                    xstart = 1 + (ielem_x-1)
                    ystart = 1 + (ielem_y-1)
                    zstart = 1 + (ielem_z-1)
                    
                    
                    ipt_elem = 1
                    do ipt_z = zstart,(zstart+1)
                        do ipt_y = ystart,(ystart+1)
                            do ipt_x = xstart,(xstart+1)
                                ipt = ipt_x  +  (ipt_y-1)*npts_x  +  (ipt_z-1)*(npts_x*npts_y)     
                                !connectivity(ielem,2+ipt_elem) = ipt
                                call connectivity%data(ielem)%set_element_node(ipt_elem, ipt)
                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do

                    ielem = ielem + 1
                end do
            end do
        end do



    end subroutine meshgen_40x15x1_linear
    !**************************************************************************************













    !> Generate a set of points defining a 4x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_15x15x1_linear(nodes,connectivity)
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, &
                       npts_y, npts_z, npts, nelem, npts_element
        integer(ik) :: ielem_x, ielem_y, ielem_z, ipt_x, ipt_y, ipt_z, &
                       ipt_elem, ielem, xstart, ystart, zstart

        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (15x15x1) - linear

        npts_x = 16
        npts_y = 16
        npts_z = 2

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 1.0_rk

        ! Allocate point storage
        npts  = npts_x*npts_y*npts_z
        nelem = (npts_x-1)*(npts_y-1)*(npts_z-1)
        npts_element = 8    !linear
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError
        call connectivity%init(nelem,npts)

        ipt = 1
        z = ZERO
        do ipt_zeta = 1,npts_z
            y = ZERO
            do ipt_eta = 1,npts_y
                x = ZERO
                do ipt_xi = 1,npts_x

                    call nodes(ipt)%set(x, y, z)
                    ipt = ipt + 1

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do




        ielem = 1
        do ielem_z = 1,npts_z-1
            do ielem_y = 1,npts_y-1
                do ielem_x = 1,npts_x-1
                    call connectivity%data(ielem)%init(1)
                    connectivity%data(ielem)%data(1) = 1        ! idomain
                    connectivity%data(ielem)%data(2) = ielem    ! ielement
                    connectivity%data(ielem)%data(3) = 1        ! mapping
                    call connectivity%data(ielem)%set_element_partition(IRANK)

                    xstart = 1 + (ielem_x-1)
                    ystart = 1 + (ielem_y-1)
                    zstart = 1 + (ielem_z-1)
                    
                    
                    ipt_elem = 1
                    do ipt_z = zstart,(zstart+1)
                        do ipt_y = ystart,(ystart+1)
                            do ipt_x = xstart,(xstart+1)
                                ipt = ipt_x  +  (ipt_y-1)*npts_x  +  (ipt_z-1)*(npts_x*npts_y)     
                                !connectivity(ielem,2+ipt_elem) = ipt
                                call connectivity%data(ielem)%set_element_node(ipt_elem, ipt)
                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do

                    ielem = ielem + 1
                end do
            end do
        end do



    end subroutine meshgen_15x15x1_linear
    !***************************************************************************************












    !> Generate a set of points defining a 15x15x2 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_15x15x2_linear(nodes,connectivity)
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, npts_y, &
                       npts_z, npts, nelem, npts_element
        integer(ik) :: ielem_x, ielem_y, ielem_z, ipt_x, ipt_y, ipt_z, ipt_elem, &
                       ielem, xstart, ystart, zstart

        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (15x15x2) - linear

        npts_x = 16
        npts_y = 16
        npts_z = 3

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 0.5_rk

        ! Allocate point storage
        npts  = npts_x*npts_y*npts_z
        nelem = (npts_x-1)*(npts_y-1)*(npts_z-1)
        npts_element = 8    !linear
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError
        call connectivity%init(nelem,npts)

        ipt = 1
        z = ZERO
        do ipt_zeta = 1,npts_z
            y = ZERO
            do ipt_eta = 1,npts_y
                x = ZERO
                do ipt_xi = 1,npts_x

                    call nodes(ipt)%set(x, y, z)
                    ipt = ipt + 1

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do





        ielem = 1
        do ielem_z = 1,npts_z-1
            do ielem_y = 1,npts_y-1
                do ielem_x = 1,npts_x-1
                    call connectivity%data(ielem)%init(1)
                    connectivity%data(ielem)%data(1) = 1        ! idomain
                    connectivity%data(ielem)%data(2) = ielem    ! ielement
                    connectivity%data(ielem)%data(3) = 1        ! mapping
                    call connectivity%data(ielem)%set_element_partition(IRANK)

                    xstart = 1 + (ielem_x-1)
                    ystart = 1 + (ielem_y-1)
                    zstart = 1 + (ielem_z-1)
                    
                    
                    ipt_elem = 1
                    do ipt_z = zstart,(zstart+1)
                        do ipt_y = ystart,(ystart+1)
                            do ipt_x = xstart,(xstart+1)
                                ipt = ipt_x  +  (ipt_y-1)*npts_x  +  (ipt_z-1)*(npts_x*npts_y)     
                                !connectivity(ielem,2+ipt_elem) = ipt
                                call connectivity%data(ielem)%set_element_node(ipt_elem, ipt)
                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do

                    ielem = ielem + 1
                end do
            end do
        end do







    end subroutine meshgen_15x15x2_linear
    !***************************************************************************************



















    !> Generate a set of points defining a 15x15x3 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_15x15x3_linear(nodes,connectivity)
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, npts_x, npts_y, &
                       npts_z, npts, nelem, npts_element
        integer(ik) :: ielem_x, ielem_y, ielem_z, ipt_x, ipt_y, ipt_z, ipt_elem, &
                       ielem, xstart, ystart, zstart

        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (15x15x3) - linear

        npts_x = 16
        npts_y = 16
        npts_z = 4

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 0.5_rk

        ! Allocate point storage
        npts  = npts_x*npts_y*npts_z
        nelem = (npts_x-1)*(npts_y-1)*(npts_z-1)
        npts_element = 8    !linear
        !allocate(nodes(npts), connectivity(nelem,2+npts_element), stat=ierr)
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError
        call connectivity%init(nelem,npts)

        ipt = 1
        z = ZERO
        do ipt_zeta = 1,npts_z
            y = ZERO
            do ipt_eta = 1,npts_y
                x = ZERO
                do ipt_xi = 1,npts_x

                    call nodes(ipt)%set(x, y, z)
                    ipt = ipt + 1

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do



        ielem = 1
        do ielem_z = 1,npts_z-1
            do ielem_y = 1,npts_y-1
                do ielem_x = 1,npts_x-1
                    call connectivity%data(ielem)%init(1)
                    connectivity%data(ielem)%data(1) = 1        ! idomain
                    connectivity%data(ielem)%data(2) = ielem    ! ielement
                    connectivity%data(ielem)%data(3) = 1        ! mapping
                    call connectivity%data(ielem)%set_element_partition(IRANK)

                    xstart = 1 + (ielem_x-1)
                    ystart = 1 + (ielem_y-1)
                    zstart = 1 + (ielem_z-1)
                    
                    
                    ipt_elem = 1
                    do ipt_z = zstart,(zstart+1)
                        do ipt_y = ystart,(ystart+1)
                            do ipt_x = xstart,(xstart+1)
                                ipt = ipt_x  +  (ipt_y-1)*npts_x  +  (ipt_z-1)*(npts_x*npts_y)     
                                call connectivity%data(ielem)%set_element_node(ipt_elem, ipt)
                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do

                    ielem = ielem + 1
                end do
            end do
        end do




    end subroutine meshgen_15x15x3_linear
    !**************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine create_mesh_file_D2E8M1_overlapping_matching(filename)
        character(*),   intent(in)  :: filename

        class(bc_state_t),  allocatable                 :: bc_state
        character(8)                                    :: faces(5)
        integer(HID_T)                                  :: file_id, dom1_id, dom2_id, bcface_id
        integer(ik)                                     :: spacedim, mapping, bcface, ierr
        type(point_t),  allocatable                     :: nodes1(:), nodes2(:)
        integer(ik),    allocatable                     :: elements1(:,:), elements2(:,:) 
        integer(ik),    allocatable                     :: faces1(:,:), faces2(:,:)
        real(rk),       allocatable, dimension(:,:,:)   :: xcoords1, xcoords2, ycoords, zcoords
        real(rk)                                        :: xmax


        ! Create/initialize file
        file_id = initialize_file_hdf(filename)
        

        ! Generate coordinates for first block
!        call meshgen_2x2x2_linear(xcoords,ycoords,zcoords)


        !
        ! Create second block by copying and translating first block.
        !
        xmax = maxval(xcoords1)
        xcoords2 = xcoords1 + 0.95*xmax


        !
        ! Get nodes/elements
        !
        mapping = 1
        nodes1    = get_block_points_plot3d(xcoords1,ycoords,zcoords)
        nodes2    = get_block_points_plot3d(xcoords2,ycoords,zcoords)
        elements1 = get_block_elements_plot3d(xcoords1,ycoords,zcoords,mapping,idomain=1)
        elements2 = get_block_elements_plot3d(xcoords2,ycoords,zcoords,mapping,idomain=2)


        !
        ! Add domains
        !
        spacedim = 3
        call add_domain_hdf(file_id,"01",nodes1,elements1,"Linear Advection",spacedim)
        call add_domain_hdf(file_id,"02",nodes2,elements2,"Linear Advection",spacedim)


        !
        ! Set boundary conditions patch connectivities
        !
        dom1_id = open_domain_hdf(file_id,"01")
        dom2_id = open_domain_hdf(file_id,"02")

        do bcface = 1,6
            ! Get face node indices for boundary 'bcface'
            faces1 = get_block_boundary_faces_plot3d(xcoords1,ycoords,zcoords,mapping,bcface)
            faces2 = get_block_boundary_faces_plot3d(xcoords2,ycoords,zcoords,mapping,bcface)

            ! Set bc patch face indices
            call set_bc_patch_hdf(dom1_id,faces1,bcface)
            call set_bc_patch_hdf(dom2_id,faces2,bcface)
        end do !bcface


        !
        ! Create bc_state, "Scalar Extrapolate"
        !
        call create_bc("Scalar Extrapolate", bc_state)


        !
        ! Set boundary condition states for Domain 1: Leave XI_MAX empty
        !
        faces = ["XI_MIN  ", "ETA_MIN ", "ETA_MAX ", "ZETA_MIN", "ZETA_MAX"]
        do bcface = 1,size(faces)
            call h5gopen_f(dom1_id,"BoundaryConditions/"//trim(adjustl(faces(bcface))),bcface_id,ierr)
            call add_bc_state_hdf(bcface_id,bc_state)
            call h5gclose_f(bcface_id,ierr)
        end do

        !
        ! Set boundary condition states for Domain 1: Leave XI_MIN empty
        !
        faces = ["XI_MAX  ", "ETA_MIN ", "ETA_MAX ", "ZETA_MIN", "ZETA_MAX"]
        do bcface = 1,size(faces)
            call h5gopen_f(dom2_id,"BoundaryConditions/"//trim(adjustl(faces(bcface))),bcface_id,ierr)
            call add_bc_state_hdf(bcface_id,bc_state)
            call h5gclose_f(bcface_id,ierr)
        end do


        ! Set 'Contains Grid'
        call set_contains_grid_hdf(file_id,"True")

        ! Close file
        call close_domain_hdf(dom1_id)
        call close_domain_hdf(dom2_id)
        call close_file_hdf(file_id)

    end subroutine create_mesh_file_D2E8M1_overlapping_matching
    !*************************************************************************************










end module mod_testutils
