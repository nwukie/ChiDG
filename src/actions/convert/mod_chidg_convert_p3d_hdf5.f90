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
    use mod_hdf_utilities,  only: initialize_file_hdf, set_ndomains_hdf, set_domain_index_hdf, &
                                  set_domain_mapping_hdf, set_domain_dimensionality_hdf, set_domain_equation_set_hdf, &
                                  set_contains_grid_hdf, set_domain_coordinates_hdf, set_domain_elements_hdf, &
                                  set_bc_patch_hdf
    use mod_plot3d_utilities,   only: get_block_elements_plot3d, get_block_boundary_faces_plot3d, &
                                      check_block_mapping_conformation_plot3d
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
        integer(HID_T)              :: file_id, block_id,  bc_id

        ! Plot3d vars
        integer(ik)                 :: i, j, k, ext_loc, fileunit, bcface
        integer(ik)                 :: npts, npt_i, npt_j, npt_k
        integer(ik)                 :: ierr, igrid, nblks, mapping, spacedim
        integer(ik),    allocatable :: blkdims(:,:)
        real(rdouble),  allocatable :: xcoords(:,:,:), ycoords(:,:,:), zcoords(:,:,:)
        integer,        allocatable :: elements(:,:), faces(:,:)

        ! equation set string
        character(len=1024)         :: eqnset_string


    
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
        file_id = initialize_file_hdf(file_prefix)


        !
        ! Open plot3d grid, read number of domains
        !
        open(newunit=fileunit, file=trim(filename), form='unformatted')
        read(fileunit) nblks
        call write_line(nblks," grid blocks", delimiter=" ")


        !
        ! Set number of domains
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



            !
            ! Dimensions for reading plot3d grid
            !
            npt_i = blkdims(1,igrid)
            npt_j = blkdims(2,igrid)
            npt_k = blkdims(3,igrid)
            npts  = blkdims(1,igrid) * blkdims(2,igrid) * blkdims(3,igrid)


            !
            ! Read block coordinates
            !
            if (allocated(xcoords)) deallocate(xcoords,ycoords,zcoords)
            allocate(xcoords(npt_i,npt_j,npt_k),ycoords(npt_i,npt_j,npt_k),zcoords(npt_i,npt_j,npt_k), stat=ierr)
            if (ierr /= 0) stop "memory allocation error: plot3d_to_hdf5"

            read(fileunit) ((( xcoords(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                           ((( ycoords(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                           ((( zcoords(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k)


            !
            ! Check mesh conforms to agglomeration routine for higher-order elements
            !
            call check_block_mapping_conformation_plot3d(xcoords,ycoords,zcoords,mapping)


            !
            ! Create a domain-group for the current block domain
            !
            write(blockname, '(I2.2)') igrid
            blockgroup = "D_"//trim(blockname)
            call h5gcreate_f(file_id, trim(blockgroup), block_id, ierr)
            if (ierr /= 0) stop "Error: h5gcreate_f"


            !
            ! Write domain attributes
            !
            call set_domain_index_hdf(block_id,igrid)
            call set_domain_mapping_hdf(block_id,mapping)
            call set_domain_dimensionality_hdf(block_id, spacedim)


            !
            ! Write coordinates
            !
            call set_domain_coordinates_hdf(block_id,xcoords,ycoords,zcoords)


            !
            ! Generate and set element connectivities
            !
            elements = get_block_elements_plot3d(xcoords,ycoords,zcoords,mapping,igrid)
            call set_domain_elements_hdf(block_id,elements)


            !
            ! Read equation set from user
            !
            call write_line("Setting equation set for domain ", igrid, delimiter=" ")
            call write_line("Enter equation set: ")
            read(*,"(A1024)") eqnset_string


            !
            ! Write equation set attribute
            !
            call set_domain_equation_set_hdf(block_id,trim(eqnset_string))



            !
            ! Create a boundary condition-group within the current block domain
            !
            call h5gcreate_f(block_id, "BoundaryConditions", bc_id, ierr)
            if (ierr /= 0) stop "Error: h5gcreate_f"
            call h5gclose_f(bc_id, ierr)
            if (ierr /= 0) stop "Error: h5gclose_f"



            !
            ! Generate and write boundary condition connectivities
            !
            do bcface = 1,6

                ! Get face node indices for boundary 'bcface'
                faces = get_block_boundary_faces_plot3d(xcoords,ycoords,zcoords,mapping,bcface)

                ! Set bc patch face indices
                call set_bc_patch_hdf(block_id,faces,bcface)

            end do !bcface


            !
            ! Close block
            !
            call h5gclose_f(block_id, ierr)
            if (ierr /= 0) stop "Error: h5gclose_f"

        end do !igrid



        !
        ! Set 'Contains Grid'
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
