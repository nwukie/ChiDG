!>  Utility to convert a plot3d grid to an hdf5 file
!!  which is read by the ChiDG code
!!
!!  @author Nathan A. Wukie
!!  @date   5/18/2016
!!
!!
!!
!!
!--------------------------------------------------------------------------------------------
module mod_chidg_convert_p3d_hdf5
#include <messenger.h>
    use mod_kinds,          only: rk,ik, rdouble
    use mod_constants,      only: IO_DESTINATION, XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use mod_hdf_utilities,  only: initialize_chidg_file_hdf, set_storage_version_major_hdf, &
                                  set_storage_version_minor_hdf, set_ndomains_hdf, set_domain_index_hdf, &
                                  set_domain_mapping_hdf, set_domain_dimensionality_hdf, set_domain_equation_set_hdf, &
                                  set_contains_grid_hdf, STORAGE_FORMAT_MAJOR, STORAGE_FORMAT_MINOR
    use hdf5
    use h5lt
    implicit none





contains


    !>  Specific routine for converting plot3d grid to hdf5
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/18/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine chidg_convert_p3d_hdf5(filename)
        character(*),   intent(in)  :: filename


        ! File, group vars
        character(1024)             :: file_prefix, hdf_file, blockgroup, blockname
        logical                     :: file_exists

        ! HDF5 vars
        integer(HID_T)              :: file_id,   Block_id,  Grid_id,  BC_id
        integer(HID_T)              :: xspace_id, yspace_id, zspace_id
        integer(HID_T)              :: xset_id,   yset_id,   zset_id
        integer(HID_T)              :: element_space_id, element_set_id
        integer(HID_T)              :: ximin_id,  ximax_id,  etamin_id,  etamax_id,  zetamin_id,  zetamax_id
        integer(HID_T)              :: bc_ximin_space_id, bc_ximax_space_id, bc_etamin_space_id, bc_etamax_space_id, &
                                       bc_zetamin_space_id, bc_zetamax_space_id
        integer(HID_T)              :: bc_ximin_set_id, bc_ximax_set_id, bc_etamin_set_id, bc_etamax_set_id, &
                                       bc_zetamin_set_id, bc_zetamax_set_id
        integer(HSIZE_T)            :: dims_rank_one(1), dims_rank_two(2), adim

        ! Plot3d vars
        integer(ik)                 :: i, j, k, ext_loc, fileunit, ipt_elem, ibc, iface, iface_i, iface_j, iface_k, ipt_face
        integer(ik)                 :: npts, npt_i, npt_j, npt_k, npts_1d
        integer(ik)                 :: ielem, ielem_i, ielem_j, ielem_k, istart_i, istart_j, istart_k
        integer(ik)                 :: nelem_i, nelem_j, nelem_k
        integer(ik)                 :: nfaces_xi, nfaces_eta, nfaces_zeta
        integer(ik)                 :: ierr, igrid, nblks, mapping, spacedim, ipt, ipt_i, ipt_j, ipt_k
        integer(ik)                 :: info_size, npts_element, npts_face
        integer(ik)                 :: nelem, pointstart_i, pointend_i, pointstart_j, pointend_j, pointstart_k, pointend_k
        integer(ik),    allocatable :: blkdims(:,:)
        real(rdouble),  allocatable :: xcoords(:,:,:), ycoords(:,:,:), zcoords(:,:,:)
        real(rdouble),  allocatable :: xcoords_linear(:), ycoords_linear(:), zcoords_linear(:)
        integer,        allocatable :: elements(:,:), faces(:,:)

        ! equation set string
        character(len=1024)         :: eqnset_string

        ! boundary condition types
        integer(ik)                 :: ximin_bc, ximax_bc, etamin_bc, etamax_bc, zetamin_bc, zetamax_bc


        !
        ! Print ChiDG Header
        !
        !call print_header()
    
        !
        ! Send output to screen, not file.
        !
        IO_DESTINATION = 'screen'



        !
        ! Check if input file exists
        !
        inquire(file=filename, exist=file_exists)
        if (file_exists) then
            call write_line("Found "//trim(filename))
            ext_loc = index(filename,'.')           ! get location of extension
            file_prefix = filename(1:(ext_loc-1))   ! save file prefix w/o extension
            hdf_file = trim(file_prefix)//'.h5'     ! set hdf5 filename with extension
        else
            call chidg_signal_one(FATAL,"File not found.",filename)
        end if


        !
        ! Initialize HDF5
        !
        call h5open_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"Error: h5open_f")


        !
        ! Create base ChiDG file and get file identifier
        !
        file_id = initialize_chidg_file_hdf(file_prefix)


        !
        ! Add file major.minor version numbers as attributes
        !
        call set_storage_version_major_hdf(file_id,STORAGE_FORMAT_MAJOR)
        call set_storage_version_minor_hdf(file_id,STORAGE_FORMAT_MINOR)




        !
        ! Read plot3d grid
        !
        open(newunit=fileunit, file=trim(filename), form='unformatted')
        read(fileunit) nblks
        call write_line(nblks," grid blocks", delimiter=" ")



        !
        ! Add number of grid domains as attribute
        !
        call set_ndomains_hdf(file_id,nblks)



        !
        ! Make space for storing dimensions of each block domain from Plot3D file.
        !
        allocate(blkdims(3,nblks),stat=ierr)
        if (ierr /= 0) call AllocationError

        !
        ! Read index dimensions from Plot3D file for each block
        !
        read(fileunit) (blkdims(1,igrid), blkdims(2,igrid), blkdims(3,igrid), igrid=1,nblks)



        !
        ! Loop through grid domain and for each domain, create an HDF5 group (D_$BLOCKNAME)
        !
        do igrid = 1,nblks

            ! Read spacedim from user
            spacedim = 0
            do while ( (spacedim < 2) .or. (spacedim > 3) )
                call write_line("Enter number of spatial dimensions for block: ", igrid, delimiter=" ")
                call write_line("Key -- ( 2 = 2D, 3 = 3D )")
                read*, spacedim
            end do


            ! Read mapping from user
            call write_line("Enter mapping for block: ", igrid, delimiter=" ")
            call write_line("Key -- (1 = linear, 2 = quadratic, 3 = cubic, 4 = quartic, 5 = quintic, 6 = sextic, 7 = septic )")
            read*, mapping
            npts_1d      = mapping+1
            npts_face    = npts_1d * npts_1d
            npts_element = npts_1d * npts_1d * npts_1d



            !
            ! Dimensions for reading plot3d grid
            !
            npt_i = blkdims(1,igrid)
            npt_j = blkdims(2,igrid)
            npt_k = blkdims(3,igrid)
            npts  = blkdims(1,igrid) * blkdims(2,igrid) * blkdims(3,igrid)



            !
            ! Test that mesh conforms to element mapping via agglomeration
            !
            !
            ! Count number of elements in each direction and check mesh conforms to
            ! the agglomeration rule for higher-order elements
            !
            nelem_i = 0
            ipt = 1
            do while (ipt < npt_i)
                nelem_i = nelem_i + 1
                ipt = ipt + (npts_1d-1)
            end do
            if (ipt > npt_i) call chidg_signal(FATAL,"Mesh does not conform to agglomeration routine in 'i'")

            nelem_j = 0
            ipt = 1
            do while (ipt < npt_j)
                nelem_j = nelem_j + 1
                ipt = ipt + (npts_1d-1)
            end do
            if (ipt > npt_j) call chidg_signal(FATAL,"Mesh does not conform to agglomeration routine in 'j'")

            nelem_k = 0
            ipt = 1
            do while (ipt < npt_k)
                nelem_k = nelem_k + 1
                ipt = ipt + (npts_1d-1)
            end do
            if (ipt > npt_k) call chidg_signal(FATAL,"Mesh does not conform to agglomeration routine in 'k'")



            !
            ! Dimensions for writing HDF5 grid
            !
            dims_rank_one = npts



            !
            ! Compute number of elements in current block
            !
            nelem_i = (npt_i-1)/mapping
            nelem_j = (npt_j-1)/mapping
            nelem_k = (npt_k-1)/mapping
            nelem   = nelem_i * nelem_j * nelem_k


            allocate(xcoords(npt_i,npt_j,npt_k),ycoords(npt_i,npt_j,npt_k),zcoords(npt_i,npt_j,npt_k),   &
                     xcoords_linear(npts),   ycoords_linear(npts),   zcoords_linear(npts), stat=ierr)
            if (ierr /= 0) stop "memory allocation error: plot3d_to_hdf5"

            read(fileunit) ((( xcoords(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                           ((( ycoords(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                           ((( zcoords(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k)

            !
            ! Create a domain-group for the current block domain
            !
            write(blockname, '(I2.2)') igrid
            blockgroup = "D_"//trim(blockname)
            call h5gcreate_f(file_id, trim(blockgroup), Block_id, ierr)
            if (ierr /= 0) stop "Error: h5gcreate_f"


            !
            ! Write domain attributes
            !
            call set_domain_index_hdf(Block_id,igrid)
            call set_domain_mapping_hdf(Block_id,mapping)
            call set_domain_dimensionality_hdf(Block_id, spacedim)


            !
            ! Create a grid-group within the current block domain
            !
            call h5gcreate_f(Block_id, "Grid", Grid_id, ierr)
            if (ierr /= 0) stop "Error: h5gcreate_f"


            !
            ! Create dataspaces for grid coordinates
            !
            call h5screate_simple_f(1, dims_rank_one, xspace_id, ierr)
            if (ierr /= 0) stop "Error: h5screate_simple_f"
            call h5screate_simple_f(1, dims_rank_one, yspace_id, ierr)
            if (ierr /= 0) stop "Error: h5screate_simple_f"
            call h5screate_simple_f(1, dims_rank_one, zspace_id, ierr)
            if (ierr /= 0) stop "Error: h5screate_simple_f"


            !
            ! Create datasets for grid coordinates
            !
            call h5dcreate_f(Grid_id, 'CoordinateX', H5T_NATIVE_DOUBLE, xspace_id, xset_id, ierr)
            if (ierr /= 0) stop "Error: h5dcreate_f"
            call h5dcreate_f(Grid_id, 'CoordinateY', H5T_NATIVE_DOUBLE, yspace_id, yset_id, ierr)
            if (ierr /= 0) stop "Error: h5dcreate_f"
            call h5dcreate_f(Grid_id, 'CoordinateZ', H5T_NATIVE_DOUBLE, zspace_id, zset_id, ierr)
            if (ierr /= 0) stop "Error: h5dcreate_f"




            !
            ! Re-order coordinates to be linear arrays
            !
            ipt = 1
            do ipt_k = 1,npt_k
                do ipt_j = 1,npt_j
                    do ipt_i = 1,npt_i
                        xcoords_linear(ipt) = xcoords(ipt_i,ipt_j,ipt_k)
                        ycoords_linear(ipt) = ycoords(ipt_i,ipt_j,ipt_k)
                        zcoords_linear(ipt) = zcoords(ipt_i,ipt_j,ipt_k)
                        ipt = ipt + 1
                    end do
                end do
            end do



            !
            ! Write coordinates to datasets
            !
            call h5dwrite_f(xset_id, H5T_NATIVE_DOUBLE, xcoords_linear, dims_rank_one, ierr)
            if (ierr /= 0) stop "Error: h5dwrite_f"
            call h5dwrite_f(yset_id, H5T_NATIVE_DOUBLE, ycoords_linear, dims_rank_one, ierr)
            if (ierr /= 0) stop "Error: h5dwrite_f"
            call h5dwrite_f(zset_id, H5T_NATIVE_DOUBLE, zcoords_linear, dims_rank_one, ierr)
            if (ierr /= 0) stop "Error: h5dwrite_f"




            !
            ! Generate element connectivities
            !
            info_size = 3               ! idomain, ielem, elem_type, ipt_1, ipt_2, ipt_3, ...
            dims_rank_two(1) = nelem
            dims_rank_two(2) = info_size + npts_element



            allocate(elements(dims_rank_two(1), dims_rank_two(2)), stat=ierr)
            if (ierr /= 0) call AllocationError

            ielem = 1
            do ielem_k = 1,nelem_k
                do ielem_j = 1,nelem_j
                    do ielem_i = 1,nelem_i

                        ! Set element info
                        elements(ielem,1) = igrid
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




            !
            ! Create dataspace for element connectivity
            !
            call h5screate_simple_f(2, dims_rank_two, element_space_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5screate_simple_f")




            !
            ! Create dataset for element connectivity
            !
            call h5dcreate_f(Grid_id, 'Elements', H5T_NATIVE_INTEGER, element_space_id, element_set_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5dcreate_f")



            !
            ! Write element connectivities
            !
            call h5dwrite_f(element_set_id, H5T_NATIVE_INTEGER, elements, dims_rank_two, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5dwrite_f")



            !
            ! Read equationset from user
            !
            call write_line("Setting equation set for domain ", igrid, delimiter=" ")
            call write_line("Enter equation set: ")
            !read*, eqnset_string
            read(*,"(A1024)") eqnset_string


            !
            ! Write equationset attribute
            !
            call set_domain_equation_set_hdf(Block_id,trim(eqnset_string))



            !
            ! Create a boundary condition-group within the current block domain
            !
            call h5gcreate_f(Block_id, "BoundaryConditions", BC_id, ierr)
            if (ierr /= 0) stop "Error: h5gcreate_f"


            !
            ! Create empty groups for boundary conditions
            !
            call h5gcreate_f(BC_id,"XI_MIN",ximin_id, ierr)
            if (ierr /= 0) stop "Error: h5gcreate_f"
            call h5gcreate_f(BC_id,"XI_MAX",ximax_id, ierr)
            if (ierr /= 0) stop "Error: h5gcreate_f"
            call h5gcreate_f(BC_id,"ETA_MIN",etamin_id, ierr)
            if (ierr /= 0) stop "Error: h5gcreate_f"
            call h5gcreate_f(BC_id,"ETA_MAX",etamax_id, ierr)
            if (ierr /= 0) stop "Error: h5gcreate_f"
            call h5gcreate_f(BC_id,"ZETA_MIN",zetamin_id, ierr)
            if (ierr /= 0) stop "Error: h5gcreate_f"
            call h5gcreate_f(BC_id,"ZETA_MAX",zetamax_id, ierr)
            if (ierr /= 0) stop "Error: h5gcreate_f"



            !
            ! Create dataspaces for boundary condition connectivity
            !
            nfaces_xi   = nelem_j * nelem_k
            nfaces_eta  = nelem_i * nelem_k
            nfaces_zeta = nelem_i * nelem_j

            dims_rank_two(1) = nfaces_xi
            dims_rank_two(2) = npts_face
            call h5screate_simple_f(2, dims_rank_two, bc_ximin_space_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5screate_simple_f")
            call h5screate_simple_f(2, dims_rank_two, bc_ximax_space_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5screate_simple_f")

            dims_rank_two(1) = nfaces_eta
            dims_rank_two(2) = npts_face
            call h5screate_simple_f(2, dims_rank_two, bc_etamin_space_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5screate_simple_f")
            call h5screate_simple_f(2, dims_rank_two, bc_etamax_space_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5screate_simple_f")

            dims_rank_two(1) = nfaces_zeta
            dims_rank_two(2) = npts_face
            call h5screate_simple_f(2, dims_rank_two, bc_zetamin_space_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5screate_simple_f")
            call h5screate_simple_f(2, dims_rank_two, bc_zetamax_space_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5screate_simple_f")



            !
            ! Create datasets for boundary condition connectivity
            !
            call h5dcreate_f(ximin_id,   'Faces', H5T_NATIVE_INTEGER, bc_ximin_space_id,   bc_ximin_set_id,   ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5dcreate_f")
            call h5dcreate_f(ximax_id,   'Faces', H5T_NATIVE_INTEGER, bc_ximax_space_id,   bc_ximax_set_id,   ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5dcreate_f")
            call h5dcreate_f(etamin_id,  'Faces', H5T_NATIVE_INTEGER, bc_etamin_space_id,  bc_etamin_set_id,  ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5dcreate_f")
            call h5dcreate_f(etamax_id,  'Faces', H5T_NATIVE_INTEGER, bc_etamax_space_id,  bc_etamax_set_id,  ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5dcreate_f")
            call h5dcreate_f(zetamin_id, 'Faces', H5T_NATIVE_INTEGER, bc_zetamin_space_id, bc_zetamin_set_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5dcreate_f")
            call h5dcreate_f(zetamax_id, 'Faces', H5T_NATIVE_INTEGER, bc_zetamax_space_id, bc_zetamax_set_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"h5dcreate_f")





            !
            ! Generate and write boundary condition connectivities
            !
            do ibc = 1,6

                select case (ibc)
                    case (XI_MIN)
                        ipt_i   = 1
                        allocate(faces(nfaces_xi,npts_face))

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


                        dims_rank_two(1) = nfaces_xi
                        dims_rank_two(2) = npts_face
                        call h5dwrite_f(bc_ximin_set_id, H5T_NATIVE_INTEGER, faces, dims_rank_two, ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,"h5dwrite_f")
                        deallocate(faces)


                    case (XI_MAX)
                        ipt_i   = npt_i
                        allocate(faces(nfaces_xi,npts_face))

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


                        dims_rank_two(1) = nfaces_xi
                        dims_rank_two(2) = npts_face
                        call h5dwrite_f(bc_ximax_set_id, H5T_NATIVE_INTEGER, faces, dims_rank_two, ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,"h5dwrite_f")
                        deallocate(faces)



                    case (ETA_MIN)
                        ipt_j   = 1
                        allocate(faces(nfaces_eta,npts_face))

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


                        dims_rank_two(1) = nfaces_eta
                        dims_rank_two(2) = npts_face
                        call h5dwrite_f(bc_etamin_set_id, H5T_NATIVE_INTEGER, faces, dims_rank_two, ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,"h5dwrite_f")
                        deallocate(faces)


                    case (ETA_MAX)
                        ipt_j   = npt_j
                        allocate(faces(nfaces_eta,npts_face))

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


                        dims_rank_two(1) = nfaces_eta
                        dims_rank_two(2) = npts_face
                        call h5dwrite_f(bc_etamax_set_id, H5T_NATIVE_INTEGER, faces, dims_rank_two, ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,"h5dwrite_f")
                        deallocate(faces)




                    case (ZETA_MIN)
                        ipt_k   = 1
                        allocate(faces(nfaces_zeta,npts_face))

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


                        dims_rank_two(1) = nfaces_zeta
                        dims_rank_two(2) = npts_face
                        call h5dwrite_f(bc_zetamin_set_id, H5T_NATIVE_INTEGER, faces, dims_rank_two, ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,"h5dwrite_f")
                        deallocate(faces)



                    case (ZETA_MAX)
                        ipt_k   = npt_k
                        allocate(faces(nfaces_zeta,npts_face))

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


                        dims_rank_two(1) = nfaces_zeta
                        dims_rank_two(2) = npts_face
                        call h5dwrite_f(bc_zetamax_set_id, H5T_NATIVE_INTEGER, faces, dims_rank_two, ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,"h5dwrite_f")
                        deallocate(faces)

                end select





            end do !ibc









            !
            ! Close boundary condition faces connectivity datasets
            !
            call h5dclose_f(bc_ximin_set_id,   ierr)
            call h5dclose_f(bc_ximax_set_id,   ierr)
            call h5dclose_f(bc_etamin_set_id,  ierr)
            call h5dclose_f(bc_etamax_set_id,  ierr)
            call h5dclose_f(bc_zetamin_set_id, ierr)
            call h5dclose_f(bc_zetamax_set_id, ierr)

            !
            ! Close boundary condition faces connectivity dataspaces
            !
            call h5sclose_f(bc_ximin_space_id,   ierr)
            call h5sclose_f(bc_ximax_space_id,   ierr)
            call h5sclose_f(bc_etamin_space_id,  ierr)
            call h5sclose_f(bc_etamax_space_id,  ierr)
            call h5sclose_f(bc_zetamin_space_id, ierr)
            call h5sclose_f(bc_zetamax_space_id, ierr)






            !
            ! Close boundary condition groups
            !
            call h5gclose_f(ximin_id,   ierr) 
            call h5gclose_f(ximax_id,   ierr) 
            call h5gclose_f(etamin_id,  ierr) 
            call h5gclose_f(etamax_id,  ierr) 
            call h5gclose_f(zetamin_id, ierr) 
            call h5gclose_f(zetamax_id, ierr) 


            !
            ! Close datasets
            !
            call h5dclose_f(xset_id,ierr)
            if (ierr /= 0) stop "Error: h5dclose_f"
            call h5dclose_f(yset_id,ierr)
            if (ierr /= 0) stop "Error: h5dclose_f"
            call h5dclose_f(zset_id,ierr)
            if (ierr /= 0) stop "Error: h5dclose_f"
            call h5dclose_f(element_set_id,ierr)
            if (ierr /= 0) stop "Error: h5dclose_f"


            !
            ! Close dataspaces
            !
            call h5sclose_f(xspace_id,ierr)
            if (ierr /= 0) stop "Error: h5sclose_f"
            call h5sclose_f(yspace_id,ierr)
            if (ierr /= 0) stop "Error: h5sclose_f"
            call h5sclose_f(zspace_id,ierr)
            if (ierr /= 0) stop "Error: h5sclose_f"
            call h5sclose_f(element_space_id,ierr)
            if (ierr /= 0) stop "Error: h5sclose_f"


            !
            ! Close active groups
            !
            call h5gclose_f(Grid_id, ierr)
            if (ierr /= 0) stop "Error: h5gclose_f"
            call h5gclose_f(BC_id, ierr)
            if (ierr /= 0) stop "Error: h5gclose_f"
            call h5gclose_f(Block_id, ierr)
            if (ierr /= 0) stop "Error: h5gclose_f"


            deallocate(zcoords_linear,ycoords_linear,xcoords_linear,zcoords,ycoords,xcoords,elements)
        end do



        !
        ! Add grid/solution attributes. Indicating the file contains a grid, and no solution.
        !
        call set_contains_grid_hdf(file_id,"True")


        !
        ! Close files and interfaces
        !
        close(fileunit)                 ! Close plot3d file
        call h5fclose_f(file_id,ierr)   ! Close hdf5 file
        call h5close_f(ierr)            ! Close hdf5 interface


        
        !
        ! Exit message
        !
        call write_line("Saved ", trim(file_prefix)//'.h5', delimiter=" ")


    end subroutine chidg_convert_p3d_hdf5
    !*****************************************************************************************











end module mod_chidg_convert_p3d_hdf5
