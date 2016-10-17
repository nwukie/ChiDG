module mod_plot3d_utilities
#include <messenger.h>
    use mod_kinds,      only: ik, rk
    use mod_constants,  only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    implicit none





contains

    !----------------------------------------------------------------------------------------
    !!
    !!  API for extracting information from Plot3D block-structured grids
    !!
    !!  Procedures:
    !!  -----------
    !!  get_block_nodes_plot3d                    - Return nodes
    !!  get_block_elements_plot3d                 - Return element connectivities
    !!  get_block_boundary_faces_plot3d           - Return face connectivities for boundary
    !!  check_block_mapping_conformation_plot3d   - Check mesh conforms to agglomeration
    !!      
    !!
    !****************************************************************************************




    !>  Return a linear(1D) array of nodes for the grid.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/16/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_block_nodes_plot3d(xcoords,ycoords,zcoords) result(nodes)
        real(rk),   intent(in)  :: xcoords(:,:,:)
        real(rk),   intent(in)  :: ycoords(:,:,:)
        real(rk),   intent(in)  :: zcoords(:,:,:)

        type(point_t),  allocatable :: nodes(:)
        integer(ik)                 :: nnodes, i,j,k, ierr



        nnodes = size(xcoords,1)*size(xcoords,2)*size(xcoords,3)

        
        allocate(nodes(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError


        inode = 1
        do k = 1,size(xcoords,3)
            do j = 1,size(xcoords,2)
                do i = 1,size(xcoords,1)

                    nodes(inode)%c1_ = xcoords(i,j,k)    
                    nodes(inode)%c2_ = ycoords(i,j,k)    
                    nodes(inode)%c3_ = zcoords(i,j,k)    
    
                    inode = inode + 1

                end do ! i
            end do ! j
        end do ! k


    end function get_block_nodes_plot3d
    !****************************************************************************************






    !>  Given coordinate arrays for a block-structured grid, return an array of element
    !!  indices for the block.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/15/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_block_elements_plot3d(xcoords,ycoords,zcoords,mapping,idomain) result(elements)
        real(rk),       intent(in)  :: xcoords(:,:,:)
        real(rk),       intent(in)  :: ycoords(:,:,:)
        real(rk),       intent(in)  :: zcoords(:,:,:)
        integer(ik),    intent(in)  :: mapping
        integer(ik),    intent(in)  :: idomain

        integer(ik) :: npt_i, npt_j, npt_k, nelem_i, nelem_j, nelem_k, nelem, &
                       info_size, npts_1d, npts_element, ielem, ielem_i, ielem_j, ielem_k, &
                       istart_i, istart_j, istart_k, ipt_i, ipt_j, ipt_k, ipt, ipt_elem, ierr

        integer(ik), allocatable    :: elements(:,:)

        !
        ! Dimensions for reading plot3d grid
        !
        npt_i = size(xcoords,1)
        npt_j = size(xcoords,2)
        npt_k = size(xcoords,3)


        !
        ! Determine number of points for geometry features
        !
        npts_1d      = mapping+1
        npts_element = npts_1d * npts_1d * npts_1d


        !
        ! Compute number of elements in current block
        !
        nelem_i = (npt_i-1)/mapping
        nelem_j = (npt_j-1)/mapping
        nelem_k = (npt_k-1)/mapping
        nelem   = nelem_i * nelem_j * nelem_k



        !
        ! Generate element connectivities
        !
        info_size = 3  ! idomain, ielem, elem_type, ipt_1, ipt_2, ipt_3, ...


        allocate(elements(nelem, info_size+npts_element), stat=ierr)
        if (ierr /= 0) call AllocationError

        ielem = 1
        do ielem_k = 1,nelem_k
            do ielem_j = 1,nelem_j
                do ielem_i = 1,nelem_i

                    ! Set element info
                    elements(ielem,1) = idomain
                    elements(ielem,2) = ielem
                    elements(ielem,3) = mapping

                    ! Get starting point
                    istart_i = 1 + ((ielem_i-1)*mapping) 
                    istart_j = 1 + ((ielem_j-1)*mapping) 
                    istart_k = 1 + ((ielem_k-1)*mapping) 

                    !
                    ! For the current element, compute node indices
                    !
                    ipt=1       ! Global point index
                    ipt_elem=1  ! Local-element point index
                    do ipt_k = istart_k,(istart_k + mapping)
                        do ipt_j = istart_j,(istart_j + mapping)
                            do ipt_i = istart_i,(istart_i + mapping)

                                ipt = ipt_i  +  (ipt_j-1)*npt_i  +  (ipt_k-1)*(npt_i*npt_j)
                                elements(ielem,info_size+ipt_elem) = ipt

                                ipt_elem = ipt_elem + 1
                            end do
                        end do
                    end do


                    ielem = ielem + 1
                end do
            end do
        end do


    end function get_block_elements_plot3d
    !****************************************************************************************





    !>  Given coordinate arrays for a block-structured grid and a face, return
    !!  an array of face node indices for the specified boundary.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/15/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_block_boundary_faces_plot3d(xcoords,ycoords,zcoords,mapping,bcface) result(faces)
        real(rk),       intent(in)  :: xcoords(:,:,:)
        real(rk),       intent(in)  :: ycoords(:,:,:)
        real(rk),       intent(in)  :: zcoords(:,:,:)
        integer(ik),    intent(in)  :: mapping
        integer(ik),    intent(in)  :: bcface

        integer(ik), allocatable    :: faces(:,:)

        integer(ik) :: ipt_i, ipt_j, ipt_k, ipt, ipt_face, npt_i, npt_j, npt_k, npts, &
                       iface_i, iface_j, iface_k, iface, nelem_i, nelem_j, nelem_k, nelem, &
                       nfaces_xi, nfaces_eta, nfaces_zeta, npts_1d, npts_face, npts_element, &
                       pointstart_i, pointstart_j, pointstart_k 


        !
        ! Check block conforms to agglomeration routine for higher-order elements
        !
        call check_block_mapping_conformation_plot3d(xcoords,ycoords,zcoords,mapping)


        !
        ! Point dimensions for block
        !
        npt_i = size(xcoords,1)
        npt_j = size(xcoords,2)
        npt_k = size(xcoords,3)
        npts  = npt_i * npt_j * npt_k


        !
        ! Compute number of points for geometry features
        !
        npts_1d      = mapping+1
        npts_face    = npts_1d * npts_1d
        npts_element = npts_1d * npts_1d * npts_1d


        !
        ! Compute number of elements in each direction and total nelem
        !
        nelem_i = (npt_i-1)/mapping
        nelem_j = (npt_j-1)/mapping
        nelem_k = (npt_k-1)/mapping
        nelem   = nelem_i * nelem_j * nelem_k


        !
        ! Get number of faces in each direction
        !
        nfaces_xi   = nelem_j * nelem_k
        nfaces_eta  = nelem_i * nelem_k
        nfaces_zeta = nelem_i * nelem_j



        select case (bcface)
            case (XI_MIN)
                allocate(faces(nfaces_xi,npts_face))
                ipt_i   = 1
                iface = 1
                do iface_k = 1,nelem_k
                    do iface_j = 1,nelem_j
                            pointstart_j = 1 + (iface_j-1)*mapping
                            pointstart_k = 1 + (iface_k-1)*mapping
                            ipt=1       ! Global point index
                            ipt_face=1  ! Local-element point index
                            do ipt_k = pointstart_k,(pointstart_k + mapping)
                                do ipt_j = pointstart_j,(pointstart_j + mapping)
                                    ipt = ipt_i  +  (ipt_j-1)*npt_i  +  (ipt_k-1)*(npt_i*npt_j)
                                    faces(iface,ipt_face) = ipt
                                    ipt_face = ipt_face + 1
                                end do
                            end do
                            iface = iface + 1
                    end do
                end do



            case (XI_MAX)
                allocate(faces(nfaces_xi,npts_face))
                ipt_i   = npt_i
                iface = 1
                do iface_k = 1,nelem_k
                    do iface_j = 1,nelem_j
                            pointstart_j = 1 + (iface_j-1)*mapping
                            pointstart_k = 1 + (iface_k-1)*mapping
                            ipt=1       ! Global point index
                            ipt_face=1  ! Local-element point index
                            do ipt_k = pointstart_k,(pointstart_k + mapping)
                                do ipt_j = pointstart_j,(pointstart_j + mapping)
                                    ipt = ipt_i  +  (ipt_j-1)*npt_i  +  (ipt_k-1)*(npt_i*npt_j)
                                    faces(iface,ipt_face) = ipt
                                    ipt_face = ipt_face + 1
                                end do
                            end do
                            iface = iface + 1
                    end do
                end do



            case (ETA_MIN)
                allocate(faces(nfaces_eta,npts_face))
                ipt_j   = 1
                iface = 1
                do iface_k = 1,nelem_k
                    do iface_i = 1,nelem_i
                            pointstart_i = 1 + (iface_i-1)*mapping
                            pointstart_k = 1 + (iface_k-1)*mapping
                            ipt=1       ! Global point index
                            ipt_face=1  ! Local-element point index
                            do ipt_k = pointstart_k,(pointstart_k + mapping)
                                do ipt_i = pointstart_i,(pointstart_i + mapping)
                                    ipt = ipt_i  +  (ipt_j-1)*npt_i  +  (ipt_k-1)*(npt_i*npt_j)
                                    faces(iface,ipt_face) = ipt
                                    ipt_face = ipt_face + 1
                                end do
                            end do
                            iface = iface + 1
                    end do
                end do



            case (ETA_MAX)
                allocate(faces(nfaces_eta,npts_face))
                ipt_j   = npt_j
                iface = 1
                do iface_k = 1,nelem_k
                    do iface_i = 1,nelem_i
                            pointstart_i = 1 + (iface_i-1)*mapping
                            pointstart_k = 1 + (iface_k-1)*mapping
                            ipt=1       ! Global point index
                            ipt_face=1  ! Local-element point index
                            do ipt_k = pointstart_k,(pointstart_k + mapping)
                                do ipt_i = pointstart_i,(pointstart_i + mapping)
                                    ipt = ipt_i  +  (ipt_j-1)*npt_i  +  (ipt_k-1)*(npt_i*npt_j)
                                    faces(iface,ipt_face) = ipt
                                    ipt_face = ipt_face + 1
                                end do
                            end do
                            iface = iface + 1
                    end do
                end do



            case (ZETA_MIN)
                allocate(faces(nfaces_zeta,npts_face))
                ipt_k   = 1
                iface = 1
                do iface_j = 1,nelem_j
                    do iface_i = 1,nelem_i
                            pointstart_i = 1 + (iface_i-1)*mapping
                            pointstart_j = 1 + (iface_j-1)*mapping
                            ipt=1       ! Global point index
                            ipt_face=1  ! Local-element point index
                            do ipt_j = pointstart_j,(pointstart_j + mapping)
                                do ipt_i = pointstart_i,(pointstart_i + mapping)
                                    ipt = ipt_i  +  (ipt_j-1)*npt_i  +  (ipt_k-1)*(npt_i*npt_j)
                                    faces(iface,ipt_face) = ipt
                                    ipt_face = ipt_face + 1
                                end do
                            end do
                            iface = iface + 1
                    end do
                end do



            case (ZETA_MAX)
                allocate(faces(nfaces_zeta,npts_face))
                ipt_k   = npt_k
                iface = 1
                do iface_j = 1,nelem_j
                    do iface_i = 1,nelem_i
                            pointstart_i = 1 + (iface_i-1)*mapping
                            pointstart_j = 1 + (iface_j-1)*mapping
                            ipt=1       ! Global point index
                            ipt_face=1  ! Local-element point index
                            do ipt_j = pointstart_j,(pointstart_j + mapping)
                                do ipt_i = pointstart_i,(pointstart_i + mapping)
                                    ipt = ipt_i  +  (ipt_j-1)*npt_i  +  (ipt_k-1)*(npt_i*npt_j)
                                    faces(iface,ipt_face) = ipt
                                    ipt_face = ipt_face + 1
                                end do
                            end do
                            iface = iface + 1
                    end do
                end do


            case default
                call chidg_signal(FATAL, "get_block_boundary_faces_plot3d: Invalid block face to get faces from")

        end select




    end function get_block_boundary_faces_plot3d
    !****************************************************************************************









    !>  Given coordinate arrays for a block-structured grid, check that the point
    !!  counts in each direction conform to the rule for agglomerating elements
    !!  in order to create higher-order elements.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/15/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine check_block_mapping_conformation_plot3d(xcoords,ycoords,zcoords,mapping)
        real(rk),       intent(in)  :: xcoords(:,:,:)
        real(rk),       intent(in)  :: ycoords(:,:,:)
        real(rk),       intent(in)  :: zcoords(:,:,:)
        integer(ik),    intent(in)  :: mapping

        integer(ik) :: ipt, npt_i, npt_j, npt_k, npts_1d, nelem_i, nelem_j, nelem_k


        ! Point dimensions for block
        npt_i   = size(xcoords,1)
        npt_j   = size(xcoords,2)
        npt_k   = size(xcoords,3)
        npts_1d = mapping+1


        !
        ! Test that block conforms to element mapping via agglomeration
        !
        !
        ! Count number of elements in each direction and check block conforms to
        ! the agglomeration rule for higher-order elements
        !
        nelem_i = 0
        ipt = 1
        do while (ipt < npt_i)
            nelem_i = nelem_i + 1
            ipt = ipt + (npts_1d-1)
        end do
        if (ipt > npt_i) call chidg_signal(FATAL,"Block mesh does not conform to agglomeration routine in 'i'")

        nelem_j = 0
        ipt = 1
        do while (ipt < npt_j)
            nelem_j = nelem_j + 1
            ipt = ipt + (npts_1d-1)
        end do
        if (ipt > npt_j) call chidg_signal(FATAL,"Block mesh does not conform to agglomeration routine in 'j'")

        nelem_k = 0
        ipt = 1
        do while (ipt < npt_k)
            nelem_k = nelem_k + 1
            ipt = ipt + (npts_1d-1)
        end do
        if (ipt > npt_k) call chidg_signal(FATAL,"Block mesh does not conform to agglomeration routine in 'k'")



    end subroutine check_block_mapping_conformation_plot3d
    !****************************************************************************************






end module mod_plot3d_utilities
