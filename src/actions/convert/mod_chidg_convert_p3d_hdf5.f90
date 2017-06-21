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
    use mod_kinds,              only: rk,ik, rdouble
    use mod_constants,          only: IO_DESTINATION, TWO
    use mod_equations,          only: equation_builder_factory
    use mod_hdf_utilities,      only: initialize_file_hdf, open_file_hdf,                   &
                                      set_domain_equation_set_hdf, set_contains_grid_hdf,   &
                                      set_domain_coordinates_hdf, set_patch_hdf,            &
                                      create_patch_hdf, close_patch_hdf,                    &
                                      add_domain_hdf, open_domain_hdf, close_domain_hdf,    &
                                      close_file_hdf, close_hdf, open_hdf
    use mod_plot3d_utilities,   only: get_block_elements_plot3d, get_block_boundary_faces_plot3d, &
                                      check_block_mapping_conformation_plot3d, get_block_points_plot3d
    use type_chidg,             only: chidg_t
    use hdf5,                   only: HID_T
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

        ! ChiDG environment for messages and loggin
        type(chidg_t)               :: chidg

        ! File, group vars
        character(1024)             :: file_prefix, hdf_file, blockgroup, blockname
        character(len=8)            :: bc_face_strings(6)
        character(:),   allocatable :: bc_face_string
        logical                     :: file_exists

        ! HDF5 vars
        integer(HID_T)              :: file_id, dom_id, patch_id

        ! Plot3d vars
        integer(ik)                 :: i, j, k, ext_loc, fileunit, bcface
        integer(ik)                 :: npts, npt_i, npt_j, npt_k
        integer(ik)                 :: ierr, igrid, nblks, mapping, system
        integer(ik),    allocatable :: blkdims(:,:)
        real(rdouble),  allocatable :: coordsx(:,:,:), coordsy(:,:,:), coordsz(:,:,:)
        real(rdouble),  allocatable :: coords1(:,:,:), coords2(:,:,:), coords3(:,:,:)
        integer,        allocatable :: elements(:,:), faces(:,:)
        real(rk),       allocatable :: nodes(:,:)
        character(:),   allocatable :: coord_system

        ! equation set string
        character(len=1024)         :: eqnset_string


        !
        ! Initialize logs and infrastructure
        !
        call chidg%start_up('mpi')
        call chidg%start_up('core')
    
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
        call open_hdf()


        !
        ! Create base ChiDG file and get file identifier
        !
        call initialize_file_hdf(hdf_file)
        file_id = open_file_hdf(hdf_file)


        !
        ! Open plot3d grid, read number of domains
        !
        open(newunit=fileunit, file=trim(filename), form='unformatted')
        read(fileunit) nblks
        call write_line(nblks," grid blocks", delimiter=" ")


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

            ! Read mapping from user
            call write_line("Enter mapping for block: ", igrid, delimiter=" ")
            call write_line("Key -- (1 = linear, 2 = quadratic, 3 = cubic, 4 = quartic, 5 = quintic, 6 = sextic, 7 = septic )")
            read*, mapping


            ! Read coordinate system
            call write_line("Enter coordinate system to use: ")
            call write_line("Key: (1 = Cartesian, 2 = Cylindrical)")
            read*, system


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
            if (allocated(coordsx)) deallocate(coordsx,coordsy,coordsz)
            allocate(coordsx(npt_i,npt_j,npt_k),coordsy(npt_i,npt_j,npt_k),coordsz(npt_i,npt_j,npt_k), stat=ierr)
            if (ierr /= 0) stop "memory allocation error: plot3d_to_hdf5"

            read(fileunit) ((( coordsx(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                           ((( coordsy(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                           ((( coordsz(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k)



            !
            ! Transform coordinates if necessary
            !
            if (system == 1) then
                coords1 = coordsx
                coords2 = coordsy
                coords3 = coordsz
                coord_system = 'Cartesian'
            else if (system == 2) then
                coords1 = sqrt(coordsx**TWO + coordsy**TWO)
                coords2 = atan2(coordsy,coordsx)
                coords3 = coordsz
                coord_system = 'Cylindrical'
            else
                call chidg_signal(FATAL,"chidg convert: Invalid coordinate system.")
            end if





            !
            ! Check mesh conforms to agglomeration routine for higher-order elements
            !
            call check_block_mapping_conformation_plot3d(coords1,coords2,coords3,mapping)


            !
            ! Create a domain-group for the current block domain
            !
            write(blockname, '(I2.2)') igrid


            !
            ! Get nodes,elements from block
            !
            nodes    = get_block_points_plot3d(coords1,coords2,coords3)
            elements = get_block_elements_plot3d(coords1,coords2,coords3,mapping,igrid)



            !
            ! Read equation set from user
            !
            call write_line("Setting equation set for domain ", igrid, delimiter=" ")
            call write_line("Enter equation set(? to list): ")

            do
                read(*,'(A1024)', iostat=ierr) eqnset_string
                if ( (ierr/=0)  ) print*, "Invalid input. Try again :)"

                if (trim(eqnset_string) == '?') then
                    call equation_builder_factory%list()
                else
                    if ( equation_builder_factory%has(trim(eqnset_string)) ) exit
                    print*, "We didn't find '"//trim(eqnset_string)//"' registered in ChiDG :/. Try another equation and remember you can enter '?' to list the registered equation sets."
                end if

            end do


            !
            ! Add new domain to file
            !
            call add_domain_hdf(file_id,trim(blockname),nodes,elements,coord_system,trim(eqnset_string))


            !
            ! Create a boundary condition-group within the current block domain
            !
            dom_id = open_domain_hdf(file_id,trim(blockname))


            !
            ! Generate and write boundary condition connectivities
            !
            bc_face_strings = ["XI_MIN  ","XI_MAX  ","ETA_MIN ","ETA_MAX ","ZETA_MIN","ZETA_MAX"]

            do bcface = 1,6

                ! Get face node indices for boundary 'bcface'
                faces = get_block_boundary_faces_plot3d(coords1,coords2,coords3,mapping,bcface)


                ! Set bc patch face indices
                bc_face_string  = trim(bc_face_strings(bcface))
                patch_id = create_patch_hdf(dom_id,bc_face_string)
                call set_patch_hdf(patch_id,faces)
                call close_patch_hdf(patch_id)

            end do !bcface


            call close_domain_hdf(dom_id)


        end do !igrid



        !
        ! Set 'Contains Grid'
        !
        call set_contains_grid_hdf(file_id,"True")


        !
        ! Close files and interfaces
        !
        close(fileunit)
        call close_file_hdf(file_id)
        call close_hdf()

        
        !
        ! Exit message
        !
        call write_line("Saved ", trim(file_prefix)//'.h5', delimiter=" ")


    end subroutine chidg_convert_p3d_hdf5
    !*****************************************************************************************











end module mod_chidg_convert_p3d_hdf5
