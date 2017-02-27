module mod_hdf_utilities
#include <messenger.h>
    use mod_kinds,              only: rk, ik, rdouble
    use mod_constants,          only: NFACES, TWO_DIM, THREE_DIM
    use mod_file_utilities,     only: delete_file
    use mod_bc,                 only: check_bc_state_registered, create_bc
    use mod_string,             only: string_t
    use mod_function,           only: create_function
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM, GLOBAL_MASTER
    use type_function,          only: function_t
    use type_svector,           only: svector_t
    use type_bc_state,          only: bc_state_t
    use type_point,             only: point_t
    use type_file_properties,   only: file_properties_t
    use type_chidg_data,        only: chidg_data_t
    use hdf5
    use h5lt
    implicit none

    

    !
    ! HDF5 storage format
    !
    integer, parameter :: STORAGE_FORMAT_MAJOR = 1
    integer, parameter :: STORAGE_FORMAT_MINOR = 3


    ! Attribute sizes
    integer(HSIZE_T), parameter :: SIZE_ONE = 1

    
    ! HDF library status
    logical    :: HDF_is_open     = .false.
    integer    :: HDF_nfiles_open = 0


contains

    !----------------------------------------------------------------------------------------
    !!
    !!  ChiDG HDF File Format API
    !!
    !!  HDF:
    !!  ---------------------------
    !!  open_hdf
    !!  close_hdf
    !!
    !!  File:
    !!  ---------------------------
    !!  initialize_file_hdf
    !!  initialize_file_structure_hdf
    !!  open_file_hdf
    !!  close_file_hdf
    !!  check_file_storage_version_hdf
    !!
    !!  set_storage_version_major_hdf
    !!  set_storage_version_minor_hdf
    !!  get_storage_version_major_hdf
    !!  get_storage_version_minor_hdf
    !!
    !!  get_properties_hdf
    !!
    !!  set_contains_grid_hdf
    !!  get_contains_grid_hdf
    !!
    !!  set_contains_solution_hdf
    !!  get_contains_solution_hdf
    !!
    !!
    !!  Domain-level routines:
    !!  ---------------------------
    !!  add_domain_hdf
    !!  create_domain_hdf
    !!  open_domain_hdf
    !!  close_domain_hdf
    !!
    !!  set_ndomains_hdf
    !!  get_ndomains_hdf
    !!
    !!  get_domain_name_hdf
    !!  get_domain_names_hdf
    !!
    !!  set_domain_index_hdf
    !!  get_domain_index_hdf
    !!  get_domain_indices_hdf
    !!
    !!  set_domain_coordinates_hdf
    !!  get_domain_coordinates_hdf
    !!
    !!  set_domain_coordinate_system_hdf
    !!  get_domain_coordinate_system_hdf
    !!
    !!  set_domain_connectivity_hdf
    !!  get_domain_connectivity_hdf
    !!
    !!  get_domain_nelements_hdf
    !!  get_domain_nnodes_hdf
    !!
    !!  set_coordinate_order_hdf
    !!  get_coordinate_order_hdf
    !!  get_coordinate_orders_hdf
    !!
    !!  set_domain_field_order_hdf
    !!  get_field_order_hdf
    !!  get_field_orders_hdf
    !!
    !!  set_domain_dimensionality_hdf
    !!  get_domain_dimensionality_hdf
    !!  get_domain_dimensionalities_hdf
    !!
    !!  set_domain_mapping_hdf
    !!  get_domain_mapping_hdf
    !!
    !!  set_domain_equation_set_hdf
    !!  get_domain_equation_set_hdf
    !!  get_domain_equation_sets_hdf
    !!
    !!
    !!  Boundary Conditions:
    !!  ---------------------------
    !!  create_bc_state_group_hdf
    !!  get_nbc_state_groups_hdf
    !!  get_bc_state_group_names_hdf
    !!
    !!  set_bc_patch_hdf
    !!  add_bc_state_hdf
    !!  get_bc_states_hdf
    !!  add_bc_properties_hdf
    !!  get_bc_properties_hdf
    !!  get_nbc_states_hdf
    !!  get_bc_state_names_hdf
    !!
    !!  remove_bc_state_group_hdf
    !!  remove_bc_state_hdf
    !!  remove_bc_property_hdf
    !!
    !!  get_bcnames_hdf
    !!
    !!  Utilities:
    !!  ----------------------------
    !!  check_bc_state_exists_hdf
    !!  check_bc_property_exists_hdf
    !!  delete_group_attributes_hdf
    !!  check_attribute_exists_hdf
    !!  check_link_exists_hdf
    !!  check_file_exists_hdf
    !!  check_file_has_extension_hdf
    !!      
    !!
    !****************************************************************************************




    !>  Handle opening the HDF library.
    !!
    !!  This way, it will only open the library if it needs opened, incase it is called
    !!  multiple times.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/18/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine open_hdf()

        integer(ik) :: ierr

        if (.not. HDF_is_open) then
            call h5open_f(ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"open_hdf: h5open_f did not execute successfully.")
        end if

        HDF_is_open = .true.

    end subroutine open_hdf
    !***************************************************************************************





    !>  Handle closing the HDF library.
    !!
    !!  This way, it will only close the library if it hasn't already been closed, incase 
    !!  it is called multiple times.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/18/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine close_hdf()

        integer(ik) :: ierr

        if (HDF_is_open .and. (HDF_nfiles_open==0)) then
            call h5close_f(ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"close_hdf: h5close_f did not execute successfully.")
        end if

        HDF_is_open = .false.

    end subroutine close_hdf
    !***************************************************************************************






    !>  Create a ChiDG-format file, with initialized format structure.
    !!  Return an HDF file identifier
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   9/25/2016
    !!
    !!  @param  filename    String for the file to be created, without extension
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine initialize_file_hdf(filename)
        character(*),   intent(in)  :: filename

        character(:),   allocatable :: filename_init
        integer(HID_T)              :: fid
        integer(ik)                 :: ierr, loc
        logical                     :: file_exists


        !
        ! Append extension if we need to
        !
        loc = index(filename,".h5")
        if (loc == 0) then
            filename_init = trim(filename)//".h5"
        else
            filename_init = trim(filename)
        end if
        

        !
        ! Check if input file already exists
        !
        file_exists = check_file_exists_hdf(filename_init)
        if (file_exists) then
            call write_line("Found "//trim(filename_init)//" that already exists. Deleting it to create new file...")
            call delete_file(trim(filename_init))
        end if


        !
        ! Create file
        !
        call open_hdf()
        call h5fcreate_f(trim(filename_init), H5F_ACC_TRUNC_F, fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"initialize_file_hdf: Error h5fcreate_f.")
        call write_line("File created: "//trim(filename_init))





        !
        ! Set storage formet
        !
        call set_storage_version_major_hdf(fid,STORAGE_FORMAT_MAJOR)
        call set_storage_version_minor_hdf(fid,STORAGE_FORMAT_MINOR)


        !
        ! Set contains status for grid/solution 
        !
        call set_contains_grid_hdf(fid,"False")
        call set_contains_solution_hdf(fid,"False")


        !
        ! Set "ndomains"
        !
        call set_ndomains_hdf(fid,0)


        call h5fclose_f(fid,ierr)

    end subroutine initialize_file_hdf
    !****************************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/23/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine initialize_file_structure_hdf(fid,data)
        integer(HID_T),     intent(in)  :: fid
        type(chidg_data_t), intent(in)  :: data

        integer(ik)                 :: idom
        integer(HID_T)              :: domain_id
        character(:),   allocatable :: domain_name



        do idom = 1,data%ndomains()

            ! Create domain group
            domain_name = data%info(idom)%name
            call create_domain_hdf(fid,domain_name)


            ! Set additional attributes
            domain_id = open_domain_hdf(fid,trim(domain_name))
            call set_domain_dimensionality_hdf(domain_id, data%get_dimensionality())
            call set_domain_equation_set_hdf(domain_id,data%eqnset(idom)%get_name())
            call close_domain_hdf(domain_id)

        end do !idom


    end subroutine initialize_file_structure_hdf
    !***************************************************************************************










    !>  Open a ChiDG-formatted HDF file and return an HDF file identifier
    !!
    !!      - Check file existence
    !!      - Open HDF interface
    !!      - Open file
    !!      - Check version
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   10/13/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function open_file_hdf(filename) result(fid)
        character(*),   intent(in)      :: filename

        character(:),   allocatable :: filename_open, user_msg
        integer(HID_T)  :: fid
        integer         :: ierr, loc
        logical         :: file_exists

        !
        ! Append extension if we need to
        !
        loc = index(filename,".h5")
        if (loc == 0) then
            filename_open = trim(filename)//".h5"
        else
            filename_open = trim(filename)
        end if

        !  Check file exists
        inquire(file=filename_open, exist=file_exists)
        if (.not. file_exists) then
            call chidg_signal(FATAL,"open_file_hdf: Could not find file: "//trim(filename_open))
        end if



        !
        !  Open input file using default properties.
        !
        call open_hdf()
        call h5fopen_f(filename_open, H5F_ACC_RDWR_F, fid, ierr)
        user_msg = "open_file_hdf: h5fopen_f, There was an error opening the file: "//trim(filename_open)
        if (ierr /= 0) call chidg_signal(FATAL,user_msg)


        !
        ! Check file format major.minor version
        !
        call check_file_storage_version_hdf(fid)


        HDF_nfiles_open = HDF_nfiles_open + 1

    end function open_file_hdf
    !*****************************************************************************************





    


    !>  Close ChiDG-formatted HDF file
    !!
    !!      - Close file
    !!      - Close HDF interface
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/15/2016
    !!
    !----------------------------------------------------------------------------------------
    subroutine close_file_hdf(fid)
        integer(HID_T), intent(in)  :: fid

        integer :: ierr

        !  Close file and Fortran interface
        call h5fclose_f(fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"close_file_hdf: error closing file.")
    
        HDF_nfiles_open = HDF_nfiles_open - 1

        call close_hdf()

    end subroutine close_file_hdf
    !****************************************************************************************










    !>  Check the format version of a ChiDG-formatted HDF file and compare it against
    !!  the version used by the current library. 
    !!
    !!  Errors if incompatible.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/25/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine check_file_storage_version_hdf(fid)
        integer(HID_T), intent(in)  :: fid

        integer(ik)                 :: file_major_version, file_minor_version, user_option
        character(:), allocatable   :: msg
        logical                     :: read_user_input

        !
        ! Get file major.minor version
        !
        file_major_version = get_storage_version_major_hdf(fid)
        file_minor_version = get_storage_version_minor_hdf(fid)


        !
        ! Test against current version
        !
        if ( (file_major_version /= STORAGE_FORMAT_MAJOR) .or. &
             (file_minor_version /= STORAGE_FORMAT_MINOR) ) then

            msg = "The storage format of the file being worked with &
                   does not match the storage format in the ChiDG library being used. This &
                   probably means the file was generated with another version of the ChiDG &
                   library and may not be compatible with the library currently being used. &
                   You could try a few things here. 1: regenerate the with with the ChiDG &
                   library being used. 2: Use a different version of the ChiDG library &
                   that uses a storage format for the file being used. 3: Full-speed ahead! &
                   Proceed anyways and try your luck!"//NEW_LINE('A')//"     &
                   Options: Exit(1), Continue(2)."

            call chidg_signal(MSG,msg)
            !call write_line(msg)

            read_user_input = .true.
            do while(read_user_input)
                
                read(*,*) user_option

                if (user_option == 1) then
                    call chidg_abort()
                    stop
                else if (user_option == 2) then
                    read_user_input = .false.
                else
                    call write_line("Valid inputs are: 1,2")
                    read_user_input = .true.
                end if

            end do


        end if
        

    end subroutine check_file_storage_version_hdf
    !*****************************************************************************************







    !>  Set the major version index for a ChiDG-formatted HDF file
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/24/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_storage_version_major_hdf(fid,major_version)
        integer(HID_T), intent(in)  :: fid
        integer(ik),    intent(in)  :: major_version

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(fid, "/", 'Major Version', [major_version], SIZE_ONE, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_storage_version_major_hdf: error setting the major version number")

    end subroutine set_storage_version_major_hdf
    !****************************************************************************************




    !>  Set the minor version index for a ChiDG-formatted HDF file
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/24/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_storage_version_minor_hdf(fid,minor_version)
        integer(HID_T), intent(in)  :: fid
        integer(ik),    intent(in)  :: minor_version

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(fid, "/", 'Minor Version', [minor_version], SIZE_ONE, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_storage_version_minor_hdf: error setting the minor version number")

    end subroutine set_storage_version_minor_hdf
    !****************************************************************************************



    !>  Return the major version index of a ChiDG-formatted HDF file
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/24/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_storage_version_major_hdf(fid) result(major_version)
        integer(HID_T), intent(in)  :: fid

        integer(ik) :: major_version, ierr, attr_status
        integer, dimension(1) :: buf


        call check_attribute_exists_hdf(fid,"Major Version","Soft Fail",attr_status)

        if (attr_status == 0) then
            call h5ltget_attribute_int_f(fid, "/", 'Major Version', buf, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_storage_version_major_hdf: error getting the major version number")
            major_version = int(buf(1),kind=ik)
        else
            major_version = -1
        end if

    end function get_storage_version_major_hdf
    !****************************************************************************************



    !>  Return the minor version index of a ChiDG-formatted HDF file
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/24/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_storage_version_minor_hdf(fid) result(minor_version)
        integer(HID_T), intent(in)  :: fid

        integer(ik) :: minor_version, ierr, attr_status
        integer, dimension(1) :: buf


        call check_attribute_exists_hdf(fid,"Minor Version","Soft Fail",attr_status)

        if (attr_status == 0) then
            call h5ltget_attribute_int_f(fid, "/", 'Minor Version', buf, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_storage_version_minor_hdf: error getting the minor version number")
            minor_version = int(buf(1),kind=ik)
        else
            minor_version = -1
        end if

    end function get_storage_version_minor_hdf
    !****************************************************************************************










    !>  Return a properties instance containing information about an hdf file
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!
    !!  @param[in]  filename    Character string containing a filename for a file that gets interrogated
    !!  @result     prop        file_properties_t instance that gets returned with file information
    !!
    !----------------------------------------------------------------------------------------
    function get_properties_hdf(filename) result(prop)
        character(*),   intent(in)  :: filename

        integer(HID_T)              :: fid
        integer(ik)                 :: ierr, idom
        integer(ik)                 :: nterms_1d
        logical                     :: fileexists = .false.

        type(file_properties_t)     :: prop


        ! 
        !  Check file exists
        !
        inquire(file=filename, exist=fileexists)
        if (.not. fileexists) then
            call chidg_signal(FATAL,'get_properties_hdf: Could not find grid file')
        end if


        !
        !  Initialize Fortran interface.
        !
        call h5open_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'get_properties_hdf - h5open_f: HDF5 Fortran interface had an error during initialization')



        !
        !  Open input file using default properties.
        !
        call h5fopen_f(filename, H5F_ACC_RDWR_F, fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'get_properties_hdf - h5fopen_f: There was an error opening the grid file.')



        !
        ! Get number of domains
        !
        prop%ndomains = get_ndomains_hdf(fid)
        call prop%set_ndomains(prop%ndomains)


        !
        ! Check for Grid and Solution contents
        !
        prop%contains_grid      = get_contains_grid_hdf(fid)
        prop%contains_solution  = get_contains_solution_hdf(fid)


        !
        ! Get domain names
        !
        prop%domain_names = get_domain_names_hdf(fid)



        !
        ! Get order of coordinate and solution expansions
        !
        if ( prop%contains_grid ) then
            prop%order_c = get_coordinate_orders_hdf(fid,prop%domain_names)
        end if

        if ( prop%contains_solution ) then
            !prop%order_s = get_solution_orders_hdf(fid,prop%domain_names)
            prop%order_s = get_solution_orders_hdf(fid)
        end if


        !
        ! Get number of spatial dimensions
        !
        prop%spacedim = get_domain_dimensionalities_hdf(fid,prop%domain_names)




        !
        ! Compute number of terms in the polynomial expansions for each domain
        !
        do idom = 1,prop%ndomains
            

            nterms_1d = (prop%order_c(idom) + 1)
            if ( prop%spacedim(idom) == THREE_DIM ) then
                prop%nterms_c(idom) = nterms_1d * nterms_1d * nterms_1d
            else if ( prop%spacedim(idom) == TWO_DIM ) then
                prop%nterms_c(idom) = nterms_1d * nterms_1d
            end if


 
            nterms_1d = (prop%order_s(idom) + 1)
            if ( prop%spacedim(idom) == THREE_DIM ) then
                prop%nterms_s(idom) = nterms_1d * nterms_1d * nterms_1d
            else if ( prop%spacedim(idom) == TWO_DIM ) then
                prop%nterms_s(idom) = nterms_1d * nterms_1d
            end if


        end do ! idom


        !
        ! Get equation set for each domain
        !
        prop%eqnset = get_domain_equation_sets_hdf(fid, prop%domain_names)




        !
        ! Close file
        !
        call h5fclose_f(fid,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_properties_hdf: h5fclose.")        
        call h5close_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_properties_hdf: h5close.")        

    end function get_properties_hdf
    !****************************************************************************************










    !>  Add a domain to the file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine add_domain_hdf(fid,domain_name,nodes,elements,coord_system,equation_set,spacedim)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: domain_name
        type(point_t),  intent(in)  :: nodes(:)
        integer(ik),    intent(in)  :: elements(:,:)
        character(*),   intent(in)  :: coord_system
        character(*),   intent(in)  :: equation_set
        integer(ik),    intent(in)  :: spacedim


        integer(HID_T)  :: dom_id, grid_id, bc_id, var_id
        integer(ik)     :: mapping, ierr


        !
        ! Create new domain
        !
        call create_domain_hdf(fid,domain_name)


        !
        ! Open Domain
        !
        dom_id = open_domain_hdf(fid,domain_name)


        !
        ! Write domain attributes
        !
        mapping  = elements(1,3)
        call set_domain_mapping_hdf(dom_id,mapping)
        call set_domain_dimensionality_hdf(dom_id, spacedim)


        ! Set nodes
        call set_domain_coordinates_hdf(dom_id,nodes)

        ! Set elements
        call set_domain_connectivity_hdf(dom_id,elements)

        ! Set coordinate system
        call set_domain_coordinate_system_hdf(dom_id,coord_system)

        ! Write equation set attribute
        call set_domain_equation_set_hdf(dom_id,trim(equation_set))


        !
        ! Close groups
        !
        call close_domain_hdf(dom_id)


    end subroutine add_domain_hdf
    !****************************************************************************************










    !>  Create a new Domain group.
    !!
    !!  Activities:
    !!      - Create a new domain group
    !!      - Increment number of domains in the file
    !!      - Set the domain index
    !!      - Close the domain group
    !!
    !!  Convention:
    !!      - Group Prefix: D_
    !!      - Example: D_MyDomain
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/23/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine create_domain_hdf(fid,domain_name)
        integer(HID_T),     intent(in)  :: fid
        character(*),       intent(in)  :: domain_name

        integer(ik)     :: ndomains, ierr
        integer(HID_T)  :: domain_id, grid_id, bc_id, var_id
        logical         :: domain_exists

        !
        ! Check if the domain group already exists
        !
        domain_exists = check_link_exists_hdf(fid,"D_"//trim(domain_name))

        
        if (.not. domain_exists) then

            ! Create domain group
            call h5gcreate_f(fid, "D_"//trim(domain_name), domain_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"create_domain_hdf: h5gcreate_f")

            call h5gcreate_f(domain_id, "Grid", grid_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"create_domain_hdf: h5gcreate_f")
            call h5gclose_f(grid_id,ierr)

            call h5gcreate_f(domain_id, "BoundaryConditions", bc_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"create_domain_hdf: h5gcreate_f")
            call h5gclose_f(bc_id,ierr)

            call h5gcreate_f(domain_id, "Variables", var_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"create_domain_hdf: h5gcreate_f")
            call h5gclose_f(var_id,ierr)



            ! Get current number of domains, increment, and reset ndomains.
            ndomains = get_ndomains_hdf(fid)
            ndomains = ndomains + 1
            call set_ndomains_hdf(fid,ndomains)


            ! Set domain name
            call set_domain_name_hdf(domain_id,domain_name)


            ! Close domain
            call close_domain_hdf(domain_id)

        end if

    end subroutine create_domain_hdf
    !***************************************************************************************

    






    !>  Open a domain group and return HDF group identifier.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function open_domain_hdf(fid,domainname) result(dom_id)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: domainname

        integer(HID_T)  :: dom_id
        integer(ik)     :: ierr
        logical         :: exists


        ! Check exists
        exists = check_link_exists_hdf(fid,"D_"//trim(domainname))
        if (.not. exists) call chidg_signal_one(FATAL,"open_domain_hdf: Couldn't find domain in file.","D_"//trim(domainname))


        ! If so, open.
        call h5gopen_f(fid,"D_"//trim(domainname), dom_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"open_domain_hdf: Error in h5gopen_f")


    end function open_domain_hdf
    !****************************************************************************************






    !>  Close a domain group from an HDF identifier.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine close_domain_hdf(dom_id)
        integer(HID_T), intent(in)  :: dom_id

        integer(ik) :: ierr

        call h5gclose_f(dom_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL, "close_domain_hdf: Error in h5gclose_f.")

    end subroutine close_domain_hdf
    !****************************************************************************************








    !>  Given a file identifier, set the number of domains in an hdf5 file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_ndomains_hdf(fid,ndomains)
        integer(HID_T), intent(in)  :: fid
        integer(ik),    intent(in)  :: ndomains

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(fid, "/", "ndomains", [ndomains], SIZE_ONE, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_ndomains_hdf: Error h5ltget_attribute_int_f")

    end subroutine set_ndomains_hdf
    !****************************************************************************************








    !>  Given a file identifier, return the number of domains in an hdf5 file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier
    !!
    !----------------------------------------------------------------------------------------
    function get_ndomains_hdf(fid) result(ndomains)
        integer(HID_T), intent(in)  :: fid

        integer                 :: ierr
        integer(ik)             :: ndomains
        integer, dimension(1)   :: buf

        call h5ltget_attribute_int_f(fid, "/", "ndomains", buf, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_ndomains_hdf: h5ltget_attribute_int_f had a problem getting the number of domains")
        ndomains = int(buf(1), kind=ik)

    end function get_ndomains_hdf
    !****************************************************************************************












    !>  Given a domain identifier, set the domain name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !!  @param[in]  domain_id       HDF domain identifier
    !!  @param[in]  domain_name     String holding the name to be set.
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_name_hdf(domain_id,domain_name)
        integer(HID_T), intent(in)  :: domain_id
        character(*),   intent(in)  :: domain_name

        integer(ik) :: ierr

        call h5ltset_attribute_string_f(domain_id, ".", "Domain Name", trim(domain_name), ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_name_hdf: Error h5ltset_attribute_string_f")

    end subroutine set_domain_name_hdf
    !****************************************************************************************







    !>  Given a domain identifier, return the domain name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !!  @param[in]  domain_id       HDF domain identifier
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_name_hdf(domain_id) result(domain_name)
        integer(HID_T), intent(in)  :: domain_id

        character(1024)             :: temp_string
        character(:),   allocatable :: domain_name
        integer(ik)                 :: ierr

        call h5ltget_attribute_string_f(domain_id, ".", "Domain Name", temp_string, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_name_hdf: Error h5ltget_attribute_string_f")

        domain_name = trim(temp_string)

    end function get_domain_name_hdf
    !****************************************************************************************

















    !>  Return a list of domain names from an HDF5 file identifier.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_names_hdf(fid) result(names)
        integer(HID_T),     intent(in)  :: fid

        character(:),       allocatable     :: user_msg
        character(len=1024), allocatable    :: names(:)
        character(len=1024)                 :: gname
        integer(HSIZE_T)                    :: igrp
        integer                             :: ndomains, nmembers, type
        integer                             :: idom, ierr

        !
        ! Get number of domains
        !
        ndomains = get_ndomains_hdf(fid)
        allocate(names(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        !  Get number of groups in the file root
        !
        call h5gn_members_f(fid, "/", nmembers, ierr)


        !
        !  Loop through groups and read domain names
        !
        idom = 1
        do igrp = 0,nmembers-1

            ! Get group name
            !call h5gget_obj_info_idx_f(fid,"/", igrp, gname, type, ierr)
            !call h5lget_name_by_idx_f(fid,".",H5_INDEX_CRT_ORDER_F,H5_ITER_INC_F,igrp,gname,ierr)
            call h5lget_name_by_idx_f(fid,".",H5_INDEX_NAME_F,H5_ITER_INC_F,igrp,gname,ierr)
            user_msg = "get_domain_names_hdf: Error iterating through links to detect domain groups."
            if (ierr /= 0) call chidg_signal(FATAL,user_msg)

            ! Test if group is a 'Domain'
            if (gname(1:2) == 'D_') then

                ! Store name
                names(idom) = trim(gname(3:))
                idom = idom + 1

            end if

        end do


    end function get_domain_names_hdf
    !****************************************************************************************













    !>  Set status of file attribute "/Contains Grid".
    !!
    !!  'True'/'False'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid             HDF5 file identifier.
    !!  @result     grid_status     Logical indicating if ChiDG grid exists in file, fid.
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_contains_grid_hdf(fid,status_string)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: status_string

        integer(ik) :: ierr

        call h5ltset_attribute_string_f(fid, "/", "Contains Grid", status_string, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"check_contains_grid_hdf - h5ltget_attribute_int_f")

    end subroutine set_contains_grid_hdf
    !****************************************************************************************










    !>  Return status of file attribute "/Contains Grid".
    !!
    !!  Attribute: /Contains Grid
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid             HDF5 file identifier.
    !!  @result     grid_status     Logical indicating if ChiDG grid exists in file, fid.
    !!
    !---------------------------------------------------------------------------------------
    function get_contains_grid_hdf(fid) result(grid_status)
        integer(HID_T), intent(in)  :: fid

        logical                         :: grid_status, attr_exists
        character(len=10)               :: contains_grid_attr
        character(len=:), allocatable   :: msg
        integer                         :: ierr


        call check_attribute_exists_hdf(fid,"Contains Grid")


        !
        ! Get attribute for 'Contains Grid
        !
        call h5ltget_attribute_string_f(fid, "/", "Contains Grid", contains_grid_attr, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"check_contains_grid_hdf - h5ltget_attribute_int_f")


        !
        ! Test grid attribute
        !
        if (trim(contains_grid_attr) == "True" .or. trim(contains_grid_attr) == "true") then
            grid_status = .true.
        else
            grid_status = .false.
        end if



    end function get_contains_grid_hdf
    !****************************************************************************************








    !>  Set status of file attribute "/Contains Solution".
    !!
    !!  'Yes'/'No'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid             HDF5 file identifier.
    !!  @result     grid_status     Logical indicating if ChiDG grid exists in file, fid.
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_contains_solution_hdf(fid,status_string)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: status_string

        integer(ik) :: ierr

        call h5ltset_attribute_string_f(fid, "/", "Contains Solution", status_string, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"check_contains_solution_hdf - h5ltget_attribute_int_f")

    end subroutine set_contains_solution_hdf
    !****************************************************************************************





    !>  Test if the file contains a Solution.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid                 HDF5 file identifier.
    !!  @result     solution_status     Logical indicating if ChiDG solution exists in file, fid.
    !!
    !----------------------------------------------------------------------------------------
    function get_contains_solution_hdf(fid) result(solution_status)
        integer(HID_T), intent(in)  :: fid

        logical             :: solution_status, attr_exists
        character(len=10)   :: contains_solution_attr
        integer             :: ierr

        
        call check_attribute_exists_hdf(fid,"Contains Solution")



        !
        ! Get attribute for "Contains Solution"
        !
        call h5ltget_attribute_string_f(fid, "/", 'Contains Solution', contains_solution_attr, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"check_contains_solution - h5ltget_attribute_string_f")


        !
        ! Test solution attribute
        !
        if (trim(contains_solution_attr) == "True" .or. trim(contains_solution_attr) == "true") then
            solution_status = .true.
        else
            solution_status = .false.
        end if



    end function get_contains_solution_hdf
    !****************************************************************************************

















!    !>  Return a domain name given a domain index.
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   2/3/2016
!    !!
!    !!  @param[in]  fid         HDF file identifier
!    !!  @param[in]  idom_hdf    A specified domain index to be queried. This is an attribute of each domain in 
!    !!                          the HDF file, per the ChiDG convention.
!    !!
!    !----------------------------------------------------------------------------------------
!    function get_domain_name_hdf(fid,idom_hdf) result(dname)
!        integer(HID_T),     intent(in)  :: fid
!        integer(ik),        intent(in)  :: idom_hdf
!
!        character(len=1024)                 :: dname
!        character(len=1024), allocatable    :: dnames(:)
!        integer(ik),         allocatable    :: dindices(:)
!        integer                             :: iind, ndomains
!
!        !
!        ! Get number of domains, domain names, and domain indices
!        !
!        ndomains = get_ndomains_hdf(fid)
!        dnames   = get_domain_names_hdf(fid)
!        dindices = get_domain_indices_hdf(fid)
!
!
!
!        do iind = 1,ndomains
!
!            if ( dindices(iind) == idom_hdf ) then
!               dname = dnames(iind) 
!            end if
!
!        end do
!
!
!
!    end function get_domain_name_hdf
!    !****************************************************************************************








!
!    !>  Set "Domain Index" for a domain group.
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   2/3/2016
!    !!
!    !!  @param[in]  fid     HDF file identifier
!    !!
!    !----------------------------------------------------------------------------------------
!    subroutine set_domain_index_hdf(dom_id,domain_index)
!        integer(HID_T),     intent(in)  :: dom_id
!        integer(ik),        intent(in)  :: domain_index
!
!        integer(ik)         :: ierr
!
!        call h5ltset_attribute_int_f(dom_id,".","Domain Index",[domain_index],SIZE_ONE,ierr)
!        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_index_hdf: Error h5ltset_attribute_int_f")
!
!
!    end subroutine set_domain_index_hdf
!    !****************************************************************************************
!




!    !>  Return a list of domain indices from an HDF5 file identifier. This is because, the 
!    !!  current method of detecting domains by name can change the order they are detected 
!    !!  in. So, each domain is given an idomain attribute that is independent of the order of 
!    !!  discovery from the file.
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   2/3/2016
!    !!
!    !!  @param[in]  fid     HDF file identifier
!    !!
!    !----------------------------------------------------------------------------------------
!    function get_domain_index_hdf(dom_id) result(domain_index)
!        integer(HID_T),     intent(in)  :: dom_id
!
!        integer(ik) :: domain_index, ierr
!        integer, dimension(1) :: buf
!
!        call h5ltget_attribute_int_f(dom_id,".","Domain Index",buf,ierr)
!        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_index_hdf: Error h5ltget_attribute_int_f")
!
!        domain_index = int(buf(1), kind=ik)
!
!    end function get_domain_index_hdf
!    !****************************************************************************************







!    !> Return a list of domain indices from an HDF5 file identifier. This is because, 
!    !! the current method of detecting domains by name can change the order they are 
!    !! detected in. So, each domain is given an idomain attribute that is independent 
!    !! of the order of discovery from the file.
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   2/3/2016
!    !!
!    !!  @param[in]  fid     HDF file identifier
!    !!
!    !---------------------------------------------------------------------------------------
!    function get_domain_indices_hdf(fid) result(indices)
!        integer(HID_T),     intent(in)  :: fid
!
!        integer(HID_T)                          :: did
!        integer(ik),            allocatable     :: indices(:)
!        character(len=1024),    allocatable     :: names(:)
!        integer(ik)                             :: idom, ndomains, ierr
!        integer, dimension(1)                   :: buf
!        integer(HSIZE_T)                        :: adim
!        logical                                 :: attribute_exists
!
!
!        !
!        ! Get number of domains
!        !
!        ndomains = get_ndomains_hdf(fid)
!        names    = get_domain_names_hdf(fid)
!
!
!
!        !
!        ! Allocate indices
!        !
!        allocate(indices(ndomains), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!
!        !
!        !  Loop through groups and read domain names
!        !
!        idom = 1
!        do idom = 1,ndomains
!            !
!            ! Open domain group
!            !
!            call h5gopen_f(fid,"D_"//trim(adjustl(names(idom))), did, ierr)
!            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_indices_hdf: error opening domain group")
!
!            !
!            ! Get idomain attribute from fid/domain/idomain
!            !
!            call h5aexists_f(did, 'Domain Index', attribute_exists, ierr)
!
!            
!            !
!            ! If it doesn't exist, set to the current value of idom
!            !
!            adim = 1
!            if ( .not. attribute_exists ) then
!                call h5ltset_attribute_int_f(fid, "D_"//trim(adjustl(names(idom))), 'Domain Index', [idom], adim, ierr)
!                if (ierr /= 0) call chidg_signal(FATAL,"get_domain_indices_hdf: error writing an initial domain index")
!            end if
!
!
!            !
!            ! Get value that was just set to be sure. 
!            !
!            call h5ltget_attribute_int_f(fid, "D_"//trim(adjustl(names(idom))), 'Domain Index', buf, ierr)
!            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_indices_hdf: error retrieving domain indices")
!
!            !
!            ! Set value detected to indices array that will be passed back from the function
!            !
!            indices(idom) = buf(1)
!
!
!            !
!            ! Close domain
!            !
!            call h5gclose_f(did,ierr)
!            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_indices_hdf: h5gclose")
!
!        end do ! idom
!
!
!    end function get_domain_indices_hdf
!    !****************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  dom_id        HDF file identifier of a block-domain group
    !!  @param[in]  domain_mapping  Integer specifying the block-domain mapping 
    !!                              1-linear, 2-quadratic, etc.
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_mapping_hdf(dom_id,domain_mapping)
        integer(HID_T),     intent(in)  :: dom_id
        integer(ik),        intent(in)  :: domain_mapping

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(dom_id,".","Domain Mapping",[domain_mapping],SIZE_ONE,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_mapping_hdf: Error h5ltset_attribute_int_f")

    end subroutine set_domain_mapping_hdf
    !****************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_mapping_hdf(dom_id) result(domain_mapping)
        integer(HID_T),     intent(in)  :: dom_id

        integer(ik) :: domain_mapping, ierr
        integer, dimension(1) :: buf

        call h5ltget_attribute_int_f(dom_id,".","Domain Mapping",buf,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_mapping_hdf: Error h5ltget_attribute_int_f")

        domain_mapping = int(buf(1), kind=ik)

    end function get_domain_mapping_hdf
    !****************************************************************************************









    !>  For a domain, set the domain coordinates
    !!
    !!  /D_domainname/Grid/Coordinate1
    !!  /D_domainname/Grid/Coordinate2
    !!  /D_domainname/Grid/Coordinate3
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   10/15/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_coordinates_hdf(dom_id,nodes)
        integer(HID_T), intent(in)  :: dom_id
        type(point_t),  intent(in)  :: nodes(:)

        integer(HID_T)      :: grid_id, xspace_id, yspace_id, zspace_id, xset_id, yset_id, zset_id
        integer(HSIZE_T)    :: dims_rank_one(1)
        integer(ik)         :: ierr, ipt, npts
        logical             :: exists
        real(rk), allocatable, dimension(:) :: xcoords, ycoords, zcoords

        !
        ! Create a grid-group within the current block domain
        !
        exists = check_link_exists_hdf(dom_id,"Grid")
        if (exists) then
            call h5gopen_f(dom_id,"Grid", grid_id,ierr)
        else
            call h5gcreate_f(dom_id, "Grid", grid_id, ierr)
        end if
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_hdf: h5gcreate_f")




        !
        ! Re-order coordinates to be linear arrays
        !
        npts = size(nodes)
        dims_rank_one = npts
        allocate(xcoords(npts), ycoords(npts), zcoords(npts), stat=ierr)

        ipt = 1
        do ipt = 1,npts
            xcoords(ipt) = nodes(ipt)%c1_
            ycoords(ipt) = nodes(ipt)%c2_
            zcoords(ipt) = nodes(ipt)%c3_
        end do

        !
        ! Create dataspaces for grid coordinates
        !
        call h5screate_simple_f(1, dims_rank_one, xspace_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_hdf: h5screate_simple_f")
        call h5screate_simple_f(1, dims_rank_one, yspace_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_hdf: h5screate_simple_f")
        call h5screate_simple_f(1, dims_rank_one, zspace_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_hdf: h5screate_simple_f")


        !
        ! Create datasets for grid coordinates
        !
        call h5dcreate_f(grid_id, "Coordinate1", H5T_NATIVE_DOUBLE, xspace_id, xset_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_hdf: h5dcreate_f")
        call h5dcreate_f(grid_id, "Coordinate2", H5T_NATIVE_DOUBLE, yspace_id, yset_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_hdf: h5dcreate_f")
        call h5dcreate_f(grid_id, "Coordinate3", H5T_NATIVE_DOUBLE, zspace_id, zset_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_hdf: h5dcreate_f")

        !
        ! Write coordinates to datasets
        !
        call h5dwrite_f(xset_id, H5T_NATIVE_DOUBLE, xcoords, dims_rank_one, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5dwrite_f")
        call h5dwrite_f(yset_id, H5T_NATIVE_DOUBLE, ycoords, dims_rank_one, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5dwrite_f")
        call h5dwrite_f(zset_id, H5T_NATIVE_DOUBLE, zcoords, dims_rank_one, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5dwrite_f")


        !
        ! Close datasets
        !
        call h5dclose_f(xset_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5dclose_f")
        call h5dclose_f(yset_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5dclose_f")
        call h5dclose_f(zset_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5dclose_f")

        !
        ! Close dataspaces
        !
        call h5sclose_f(xspace_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5sclose_f")
        call h5sclose_f(yspace_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5sclose_f")
        call h5sclose_f(zspace_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5sclose_f")

        !
        ! Close Grid group
        !
        call h5gclose_f(grid_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5gclose_f")


    end subroutine set_domain_coordinates_hdf
    !***************************************************************************************










    !>  For a domain, get the domain coordinates
    !!
    !!  /D_domainname/Grid/Coordinate1
    !!  /D_domainname/Grid/Coordinate2
    !!  /D_domainname/Grid/Coordinate3
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   02/14/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_coordinates_hdf(dom_id) result(nodes)
        integer(HID_T), intent(in)  :: dom_id

        type(point_t), allocatable  :: nodes(:)

        integer(HID_T)              :: grid_id, did_1, did_2, did_3, sid
        integer(HSIZE_T)            :: rank_one_dims(1), maxdims(3)
        integer(ik)                 :: ierr, ipt, npts
        logical                     :: exists

        real(rdouble), dimension(:), allocatable, target    :: pts1, pts2, pts3
        type(c_ptr)                                         :: cp_pts, cp_conn


        !
        ! Create a grid-group within the current block domain
        !
        exists = check_link_exists_hdf(dom_id,"Grid")
        if (exists) then
            call h5gopen_f(dom_id,"Grid", grid_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinates_hdf: h5gopen_f")
        else
            call chidg_signal(FATAL,"get_domain_coordinates_hdf: The current domain does not contain a 'Grid'.")
        end if


        !
        !  Open the Coordinate datasets
        !
        call h5dopen_f(grid_id, "Coordinate1", did_1, ierr, H5P_DEFAULT_F)
        call h5dopen_f(grid_id, "Coordinate2", did_2, ierr, H5P_DEFAULT_F)
        call h5dopen_f(grid_id, "Coordinate3", did_3, ierr, H5P_DEFAULT_F)


        !
        !  Get the dataspace id and dimensions
        !
        call h5dget_space_f(did_1, sid, ierr)
        call h5sget_simple_extent_dims_f(sid, rank_one_dims, maxdims, ierr)
        npts = rank_one_dims(1)


        !
        !  Read points 1,2,3
        !
        allocate(pts1(npts), pts2(npts), pts3(npts), stat=ierr)
        if (ierr /= 0) call AllocationError


        cp_pts = c_loc(pts1(1))
        call h5dread_f(did_1, H5T_NATIVE_DOUBLE, cp_pts, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinates_hdf: h5dread_f")

        cp_pts = c_loc(pts2(1))
        call h5dread_f(did_2, H5T_NATIVE_DOUBLE, cp_pts, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinates_hdf: h5dread_f")

        cp_pts = c_loc(pts3(1))
        call h5dread_f(did_3, H5T_NATIVE_DOUBLE, cp_pts, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinates_hdf: h5dread_f")



        !
        !  Accumulate pts into a single points_t matrix to initialize domain
        !
        allocate(nodes(npts), stat=ierr)
        if (ierr /= 0) call AllocationError
            

        do ipt = 1,rank_one_dims(1)
            call nodes(ipt)%set(real(pts1(ipt),rk),real(pts2(ipt),rk),real(pts3(ipt),rk))
        end do


        ! Close the Coordinate datasets
        call h5dclose_f(did_1,ierr)
        call h5dclose_f(did_2,ierr)
        call h5dclose_f(did_3,ierr)

        ! Close the dataspace id
        call h5sclose_f(sid,ierr)

        ! Close Grid group
        call h5gclose_f(grid_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinates_her: h5gclose_f")


    end function get_domain_coordinates_hdf
    !***************************************************************************************










    !>  For a domain, set the domain coordinate system.
    !!
    !!  Attribute: /D_domainname/Grid/Coordinate System
    !!  Values: 'Cartesian' or 'Cylindrical'
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   02/14/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_coordinate_system_hdf(dom_id,coord_system)
        integer(HID_T), intent(in)  :: dom_id
        character(*),   intent(in)  :: coord_system

        character(:),   allocatable :: user_msg
        integer(HID_T)              :: grid_id
        integer(ik)                 :: ierr
        logical                     :: exists, cartesian, cylindrical



        !
        ! Check valid input
        !
        cartesian   = (coord_system == 'Cartesian')
        cylindrical = (coord_system == 'Cylindrical')
        user_msg = "set_domain_coordinate_system_hdf: Invalid coordinate system. Valid systems are 'Cartesian' or 'Cylindrical'."
        if (.not. (cartesian .or. cylindrical)) call chidg_signal_one(FATAL,user_msg,coord_system)




        !
        ! Create a grid-group within the current block domain
        !
        exists = check_link_exists_hdf(dom_id,"Grid")
        if (exists) then
            call h5gopen_f(dom_id,"Grid", grid_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinate_system_hdf: h5gopen_f")
        else
            call chidg_signal(FATAL,"set_domain_coordinate_system_hdf: The current domain does not contain a 'Grid'.")
        end if


        !
        ! Set 'Coordinate System' attribute
        !
        call h5ltset_attribute_string_f(grid_id,".","Coordinate System",coord_system,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinate_system_hdf: Error setting 'Coordinate System' attribute.")



        !
        ! Close Grid group
        !
        call h5gclose_f(grid_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinate_system_hdf: h5gclose_f")


    end subroutine set_domain_coordinate_system_hdf
    !***************************************************************************************










    !>  For a domain, get the domain coordinate system.
    !!
    !!  Attribute: /D_domainname/Grid/Coordinate System
    !!  Values: 'Cartesian' or 'Cylindrical'
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   02/14/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_coordinate_system_hdf(dom_id) result(coord_system_trim)
        integer(HID_T), intent(in)  :: dom_id

        character(1024)             :: coord_system
        character(:),   allocatable :: coord_system_trim

        integer(HID_T)      :: grid_id
        integer(ik)         :: ierr
        logical             :: exists


        !
        ! Create a grid-group within the current block domain
        !
        exists = check_link_exists_hdf(dom_id,"Grid")
        if (exists) then
            call h5gopen_f(dom_id,"Grid", grid_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_system_hdf: h5gopen_f")
        else
            call chidg_signal(FATAL,"get_domain_coordinate_system_hdf: The current domain does not contain a 'Grid'.")
        end if


        !
        ! Get 'Coordinate System' attribute
        !
        call h5ltget_attribute_string_f(grid_id,".","Coordinate System",coord_system,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_system_hdf: Error getting 'Coordinate System' attribute.")

        
        ! Trim blanks
        coord_system_trim = trim(coord_system)


        !
        ! Close Grid group
        !
        call h5gclose_f(grid_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_system_hdf: h5gclose_f")


    end function get_domain_coordinate_system_hdf
    !***************************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/25/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_coordinate_order_hdf(dom_id,order)
        integer(HID_T), intent(in)  :: dom_id
        integer(ik),    intent(in)  :: order

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(dom_id,".","Coordinate Order", [order],SIZE_ONE,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_coordinate_order_hdf: Error setting 'Coordinate Order' attribute")

    end subroutine set_coordinate_order_hdf
    !****************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/25/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine get_coordinate_order_hdf(dom_id,order)
        integer(HID_T), intent(in)  :: dom_id
        integer(ik),    intent(in)  :: order

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(dom_id,".","Coordinate Order", [order],SIZE_ONE,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_coordinate_order_hdf: Error getting 'Coordinate Order' attribute")

    end subroutine get_coordinate_order_hdf
    !****************************************************************************************








    !>  Returns an array of integers that specifies the order of the coordinate expansion 
    !!  for every domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid         HDF file identifier.
    !!  @param[in]  dnames(:)   List of domain names to be interrogated. 
    !!
    !----------------------------------------------------------------------------------------
    function get_coordinate_orders_hdf(fid, dnames) result(orders)
        integer(HID_T),         intent(in)  :: fid
        character(len=1024),    intent(in)  :: dnames(:)

        integer(ik), allocatable    :: orders(:)
        integer                     :: ierr, idom, mapping
        integer, dimension(1)       :: buf


        !
        ! Allocate storage for orders
        !
        allocate(orders(size(dnames)), stat=ierr)
        if (ierr /= 0) call AllocationError

        !
        !  Loop through groups and read domains
        !
        do idom = 1,size(dnames)

            !
            !  Get coordinate mapping
            !
            call h5ltget_attribute_int_f(fid, "D_"//trim(dnames(idom)), "Domain Mapping", buf, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_coordinate_orders_hdf: h5ltget_attribute_int_f")


            !
            ! Compute number of terms in coordinate expansion
            !
            mapping = buf(1)

            orders(idom) = int(mapping, kind=ik)

        end do

    end function get_coordinate_orders_hdf
    !****************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/25/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_solution_order_hdf(dom_id,order)
        integer(HID_T), intent(in)  :: dom_id
        integer(ik),    intent(in)  :: order

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(dom_id,".","Solution Order", [order],SIZE_ONE,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_solution_order_hdf: Error setting 'Solution Order' attribute")

    end subroutine set_solution_order_hdf
    !****************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/25/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    function get_solution_order_hdf(dom_id) result(order)
        integer(HID_T), intent(in)  :: dom_id

        integer(ik) :: order, ierr
        integer, dimension(1) :: buf

        call h5ltget_attribute_int_f(dom_id,".","Solution Order",buf,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_solution_order_hdf: Error getting 'Solution Order' attribute")

        order = int(buf(1), kind=ik)

    end function get_solution_order_hdf
    !****************************************************************************************










    !>  Returns an array of integer that specifies the order of the solution expansion for 
    !!  every domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier.
    !!  @param[in]  dnames  List of domain names to be interrogated.
    !!
    !----------------------------------------------------------------------------------------
    function get_solution_orders_hdf(fid) result(orders)
        integer(HID_T),         intent(in)  :: fid

        integer(HID_T)                      :: did
        integer(ik), allocatable            :: orders(:)
        integer                             :: ierr, idom
        character(len=1024), allocatable    :: domain_names(:)


        !
        ! Get domains
        !
        domain_names = get_domain_names_hdf(fid)


        !
        ! Allocate storage for orders
        !
        allocate(orders(size(domain_names)), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        !  Loop through groups and read domains
        !
        do idom = 1,size(domain_names)

            !
            ! Open domain group, get solution order, close domain group
            !
            did = open_domain_hdf(fid,trim(domain_names(idom)))

            orders(idom) = get_solution_order_hdf(did)

            call close_domain_hdf(did)

        end do

    end function get_solution_orders_hdf
    !****************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine set_domain_dimensionality_hdf(dom_id,dimensionality)
        integer(HID_T), intent(in)  :: dom_id
        integer(ik),    intent(in)  :: dimensionality

        integer(ik)         :: ierr

        !  Get coordinate mapping
        call h5ltset_attribute_int_f(dom_id,".","Domain Dimensionality",[dimensionality],SIZE_ONE,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dimensionality_hdf: h5ltset_attribute_int_f")

    end subroutine set_domain_dimensionality_hdf
    !****************************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_dimensionality_hdf(dom_id) result(dimensionality)
        integer(HID_T), intent(in)  :: dom_id

        integer(ik) :: dimensionality, ierr
        integer, dimension(1)   :: buf

        !  Get coordinate mapping
        call h5ltget_attribute_int_f(dom_id,".","Domain Dimensionality",buf,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_dimensionality_hdf: h5ltget_attribute_int_f")

        dimensionality = int(buf(1),kind=ik)

    end function get_domain_dimensionality_hdf
    !****************************************************************************************










    !>  Returns an array of integers that specifies the number of spatial dimensions to use 
    !!  for every domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!  @param[in]  fid         HDF file identifier.
    !!  @param[in]  dnames(:)   List of domain names to be interrogated. 
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_dimensionalities_hdf(fid, dnames) result(dimensionalities)
        integer(HID_T),         intent(in)  :: fid
        character(len=1024),    intent(in)  :: dnames(:)

        integer(ik), allocatable    :: dimensionalities(:)
        integer(ik)                 :: ierr, idom
        integer, dimension(1)       :: dimensionality


        !
        ! Allocate storage for orders
        !
        allocate(dimensionalities(size(dnames)), stat=ierr)
        if (ierr /= 0) call AllocationError

        !
        !  Loop through groups and read domains
        !
        do idom = 1,size(dnames)

            call h5ltget_attribute_int_f(fid, "D_"//trim(dnames(idom)), "Domain Dimensionality", dimensionality, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_dimensionalities_hdf: Error h5ltget_attribute_int_f")

            dimensionalities(idom) = dimensionality(1)

        end do

    end function get_domain_dimensionalities_hdf
    !****************************************************************************************







    !>  Set the 'Equation Set' string for a domain. Specified by HDF Domain identifier.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_equation_set_hdf(dom_id,equation_set)
        integer(HID_T), intent(in)  :: dom_id
        character(*),   intent(in)  :: equation_set

        integer(ik) :: ierr

        !  Get coordinate mapping
        call h5ltset_attribute_string_f(dom_id,".","Equation Set",equation_set,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_equation_set_hdf: h5ltset_attribute_string_f")

    end subroutine set_domain_equation_set_hdf
    !****************************************************************************************





    !>  Get the 'Equation Set' string for a domain. Specified by HDF Domain identifier.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_equation_set_hdf(dom_id) result(equation_set)
        integer(HID_T), intent(in)  :: dom_id

        character(1024) :: equation_set
        integer(ik)     :: ierr

        !  Get coordinate mapping
        call h5ltget_attribute_string_f(dom_id,".","Equation Set",equation_set,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_equation_set_hdf: h5ltset_attribute_string_f")

    end function get_domain_equation_set_hdf
    !****************************************************************************************






    !>  Return an array of strings that are the 'Equation Set' attribute of every domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier.
    !!  @param[in]  dnames  List of domain names to be interrogated.
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_equation_sets_hdf(fid, dnames) result(eqnsets)
        integer(HID_T),         intent(in)  :: fid
        character(len=1024),    intent(in)  :: dnames(:)

        integer(HID_T)                      :: did
        logical                             :: eqnset_exists
        character(len=1024), allocatable    :: eqnsets(:)
        integer                             :: ierr, idom


        !
        ! Allocate storage for eqnsets
        !
        allocate(eqnsets(size(dnames)), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        !  Loop through groups and read domains
        !
        do idom = 1,size(dnames)
            !
            ! Open domain
            !
            call h5gopen_f(fid, "D_"//trim(dnames(idom)), did, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_equation_sets_hdf: error opening domain group.")


            !
            ! Check eqnset attribute exists.
            !
            call h5aexists_f(did, "Equation Set", eqnset_exists, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_equation_sets_hdf: error checking if 'Equation Set' exists.")


            !
            ! If eqnset doesn't exists, create attribute and set to 'empty'
            !
            if ( .not. eqnset_exists ) then
                call h5ltset_attribute_string_f(did, ".", "Equation Set", 'empty', ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_domain_equation_sets_hdf: error setting empty 'Equation Set' attribute.")
            end if


            !
            !  Get eqnset string from hdf attribute.
            !
            call h5ltget_attribute_string_f(did, ".", "Equation Set", eqnsets(idom), ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_equation_sets_hdf - h5ltget_attribute_int_f.")


            call h5gclose_f(did,ierr)

        end do

    end function get_domain_equation_sets_hdf
    !****************************************************************************************











    !>  Set the element connectivities for a block.
    !!
    !!  Accepts array of element connectivities as
    !!      elements(nelem, size_connectivity)
    !!
    !!  /D_domainname/Grid/Elements
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/15/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_connectivity_hdf(dom_id,elements)
        integer(HID_T), intent(in)  :: dom_id
        integer(ik),    intent(in)  :: elements(:,:)

        integer(ik)         :: ierr
        integer(HID_T)      :: element_set_id, element_space_id, grid_id
        integer(HSIZE_T)    :: dims_rank_two(2)
        logical             :: exists

        !
        ! Size element connectivities
        !
        dims_rank_two(1) = size(elements,1)
        dims_rank_two(2) = size(elements,2)


        !
        ! Create a grid-group within the current block domain
        !
        exists = check_link_exists_hdf(dom_id,"Grid")
        if (exists) then
            call h5gopen_f(dom_id, "Grid", grid_id, ierr)
        else
            call h5gcreate_f(dom_id,"Grid", grid_id, ierr)
        end if
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_hdf: h5gcreate_f/h5gopen_f")


        !
        ! Create dataset for element connectivity: element_set_id
        !
        exists = check_link_exists_hdf(grid_id,"Elements")
        if (exists) then
            call h5dopen_f(grid_id,"Elements", element_set_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_hdf: h5dopen_f")

            !
            ! Write element connectivities
            !
            call h5dwrite_f(element_set_id, H5T_NATIVE_INTEGER, elements, dims_rank_two, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_hdf: h5dwrite_f")

        else
            call h5screate_simple_f(2, dims_rank_two, element_space_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_hdf: h5screate_simple_f")

            call h5dcreate_f(grid_id, "Elements", H5T_NATIVE_INTEGER, element_space_id, element_set_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_hdf: h5dcreate_f")



            !
            ! Write element connectivities
            !
            call h5dwrite_f(element_set_id, H5T_NATIVE_INTEGER, elements, dims_rank_two, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_hdf: h5dwrite_f")

            call h5sclose_f(element_space_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_hdf: h5sclose_f")

        end if


        !
        ! Close dataset, dataspace, Grid group
        !
        call h5dclose_f(element_set_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_hdf: h5dclose_f")
        call h5gclose_f(grid_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_hdf: h5gclose_f")



    end subroutine set_domain_connectivity_hdf
    !****************************************************************************************












    !>  Get the element connectivities for a domain.
    !!
    !!  Returns array of element connectivities as
    !!      elements(nelem, size_connectivity)
    !!
    !!  /D_domainname/Grid/Elements
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/15/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_connectivity_hdf(dom_id) result(connectivity_ik)
        integer(HID_T), intent(in)  :: dom_id


        integer(ik)         :: ierr
        integer(HID_T)      :: grid_id, did_e, sid
        integer(HSIZE_T)    :: rank_two_dims(2), maxdims(3)
        logical             :: exists

        character(:),   allocatable         :: user_msg
        integer(ik),    allocatable         :: connectivity_ik(:,:)
        integer,        allocatable, target :: connectivity(:,:)
        type(c_ptr)                         :: cp_conn


        !
        ! Open 'Grid' group for the current domain
        !
        exists = check_link_exists_hdf(dom_id,'Grid')
        if (exists) then
            call h5gopen_f(dom_id, 'Grid', grid_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_connectivity_hdf: h5gopen_f")
        else
            call chidg_signal(FATAL,"get_domain_connectivity_hdf: Current domain does not have a 'Grid' group.")
        end if


        !
        ! Open Elements connectivity dataset
        !
        call h5dopen_f(grid_id, 'Elements', did_e, ierr, H5P_DEFAULT_F)
        user_msg = "get_domain_connectivity_hdf: h5dopen_f did not open 'Elements' dataset propertly."
        if (ierr /= 0) call chidg_signal(FATAL,user_msg)


        !
        !  Get the dataspace id and dimensions. (nelem, connectivity size)
        !
        call h5dget_space_f(did_e, sid, ierr)
        user_msg = "get_domain_connectivity_hdf: h5dget_space_f did not return 'Elements' dataspace propertly."
        if (ierr /= 0) call chidg_signal(FATAL,user_msg)

        call h5sget_simple_extent_dims_f(sid, rank_two_dims, maxdims, ierr)
        user_msg = "get_domain_connectivity_hdf: h5sget_simple_extent_dims_f did not return extent propertly."
        if (ierr == -1) call chidg_signal(FATAL,user_msg)


        !
        ! Allocate/Read connectivity
        ! 
        allocate(connectivity(rank_two_dims(1),rank_two_dims(2)), stat=ierr)
        if (ierr /= 0) call AllocationError

        cp_conn = c_loc(connectivity(1,1))
        call h5dread_f(did_e, H5T_NATIVE_INTEGER, cp_conn, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_connectivity_hdf5: h5dread_f")


        !
        ! Convert to ChiDG integer type
        !
        !   Explicit allocation to handle GCC bug:
        !       GCC/GFortran Bugzilla Bug 52162 
        !
        allocate(connectivity_ik(size(connectivity,1),size(connectivity,2)),stat=ierr)
        if (ierr /= 0) call AllocationError
        connectivity_ik = int(connectivity, ik)


        !
        ! Close identifiers
        !
        call h5dclose_f(did_e,  ierr)
        call h5sclose_f(sid,    ierr)
        call h5gclose_f(grid_id,ierr)


    end function get_domain_connectivity_hdf
    !****************************************************************************************











    !>  Given a domain identifier, return the number of elements in the domain.
    !!
    !!  TODO: Switch this to read an attribute. That way we can use this is 
    !!        solution files without grids.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !---------------------------------------------------------------------------------------
    function get_domain_nelements_hdf(domain_id) result(nelements)
        integer(HID_T), intent(in)  :: domain_id

        integer(ik)                 :: nelements, ierr
        integer(HID_T)              :: elements_id, space_id
        integer(HSIZE_T)            :: dims(2), maxdims(2)
        character(:), allocatable   :: user_msg
        logical                     :: grid_exists


        ! Check file has a grid.
        grid_exists = check_link_exists_hdf(domain_id,"Grid/Elements")
        user_msg = "get_domain_nelements_hdf: Trying to determine the number of elements in a &
                    domain without a 'Grid/Elements' group. This file probably doesn't contain &
                    a grid, so you will want to figure out how to remedy this."
        if (.not. grid_exists) call chidg_signal(FATAL,user_msg)


        ! Open 'Elements' data set.
        call h5dopen_f(domain_id,"Grid/Elements", elements_id, ierr)
        user_msg = "get_domain_nelements_hdf: There was an error opening the 'Grid/Elements' &
                    data set."
        if (ierr /= 0) call chidg_signal(FATAL,user_msg)


        ! Get the dataspace id.
        call h5dget_space_f(elements_id, space_id, ierr)
        user_msg = "get_domain_nelements_hdf: There was an error opening the 'Grid/Elements' &
                    data space."
        if (ierr /= 0) call chidg_signal(FATAL,user_msg)


        ! Get the data space dimensions.
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, ierr)
        user_msg = "get_domain_nelements_hdf: There was an error returning the dimensions of &
                    the 'Grid/Elements' data space."
        if (ierr == -1) call chidg_signal(FATAL,user_msg)


        ! Close groups.
        call h5sclose_f(space_id,ierr)
        call h5dclose_f(elements_id,ierr)


        ! Return number of elements. Size of first dimension.
        nelements = int(dims(1), kind=ik)


    end function get_domain_nelements_hdf
    !***************************************************************************************









    !>  Given a domain identifier, return the number of nodes in the domain.
    !!
    !!  TODO: Switch this to read an attribute. That way we can use this in 
    !!        solution files without grids.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !---------------------------------------------------------------------------------------
    function get_domain_nnodes_hdf(domain_id) result(nnodes)
        integer(HID_T), intent(in)  :: domain_id

        integer(HID_T)              :: did, sid, grid_id
        integer(HSIZE_T)            :: rank_one_dims(1), maxdims(2)
        character(:), allocatable   :: user_msg
        logical                     :: grid_exists
        integer(ik)                 :: ierr, nnodes


        !
        ! Check file has a node set 
        !
        grid_exists = check_link_exists_hdf(domain_id,'Grid/Coordinate1')
        user_msg = "get_domain_nelements_hdf: Trying to determine the number of elements in a &
                    domain without a 'Grid/Elements' group. This file probably doesn't contain &
                    a grid, so you will want to figure out how to remedy this."
        if (.not. grid_exists) call chidg_signal(FATAL,user_msg)


        !
        ! Open the Grid group
        !
        call h5gopen_f(domain_id, 'Grid', grid_id, ierr, H5P_DEFAULT_F)
        user_msg = "get_domain_nnodes_hdf: Domain/Grid group did not open properly."
        if (ierr /= 0) call chidg_signal(FATAL,user_msg)


        !
        ! Get number of nodes in the domain
        !
        call h5dopen_f(grid_id, 'Coordinate1', did, ierr, H5P_DEFAULT_F)
        user_msg = "get_domain_nnodes_hdf: Domain/Grid/Coordinate1 group did not open properly."
        if (ierr /= 0) call chidg_signal(FATAL,user_msg)


        !
        !  Get the dataspace id and dimensions
        !
        call h5dget_space_f(did, sid, ierr)
        user_msg = "get_domain_nnodes_hdf: h5dget_space_f did not return 'Coordinate1' dataspace properly."
        if (ierr /= 0) call chidg_signal(FATAL,user_msg)

        call h5sget_simple_extent_dims_f(sid, rank_one_dims, maxdims, ierr)
        user_msg = "get_domain_nnodes_hdf: h5sget_simple_extent_dims_f did not return extent propertly."
        if (ierr == -1) call chidg_signal(FATAL,user_msg)

        
        !
        ! Get node count from size of coordinate dataset in file
        !
        nnodes = int(rank_one_dims(1), ik)


        !
        ! Close identifiers
        !
        call h5sclose_f(sid,ierr)
        call h5dclose_f(did,ierr)
        call h5gclose_f(grid_id,ierr)


    end function get_domain_nnodes_hdf
    !***************************************************************************************



















    !>  Set boundary condition patch face indices for a block boundary.
    !!
    !!  /D_domainname/BoundaryConditions/"face"/Faces
    !!
    !!  "face" is XI_MIN, XI_MAX, ETA_MIN, etc.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/15/2016
    !!
    !!  @param[in]  dom_id      HDF identifier for the block to be written to
    !!  @param[in]  faces       Face indices to be set for the boundary condition patch
    !!  @param[in]  bcface      Integer specifying which boundary of the block to write to
    !!
    !---------------------------------------------------------------------------------------
    subroutine set_bc_patch_hdf(dom_id,faces,bcface)
        integer(HID_T), intent(in)  :: dom_id
        integer(ik),    intent(in)  :: faces(:,:)
        integer(ik),    intent(in)  :: bcface

        character(len=8)            :: bc_face_strings(6)
        character(:), allocatable   :: bc_face_string
        integer                     :: ierr
        integer(HID_T)              :: bc_id, face_id, face_space_id, face_set_id
        integer(HSIZE_T)            :: dims(2)
            

        !
        ! Create a boundary condition-group within the current block domain
        !
        call h5gopen_f(dom_id, "BoundaryConditions", bc_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_bc_patch_hdf: h5gcreate_f")


        !
        ! Get boundary condition string
        !
        bc_face_strings = ["XI_MIN  ","XI_MAX  ","ETA_MIN ","ETA_MAX ","ZETA_MIN","ZETA_MAX"]
        bc_face_string  = trim(bc_face_strings(bcface))


        !
        ! Create empty group for boundary condition
        !
        call h5gcreate_f(bc_id,bc_face_string,face_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_bc_patch_hdf: h5gcreate_f")


        !
        ! Create dataspaces for boundary condition connectivity
        !
        dims(1) = size(faces,1)
        dims(2) = size(faces,2)
        call h5screate_simple_f(2, dims, face_space_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_bc_patch_hdf: h5screate_simple_f")


        !
        ! Create datasets for boundary condition connectivity
        !
        call h5dcreate_f(face_id,"Faces", H5T_NATIVE_INTEGER, face_space_id, face_set_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_bc_patch_hdf: h5dcreate_f")


        !
        ! Write bc faces
        !
        call h5dwrite_f(face_set_id, H5T_NATIVE_INTEGER, faces, dims, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_bc_patch_hdf: h5dwrite_f")


        !
        ! Close groups
        !
        call h5dclose_f(face_set_id, ierr)
        call h5sclose_f(face_space_id, ierr)
        call h5gclose_f(face_id, ierr) 
        call h5gclose_f(bc_id, ierr)


    end subroutine set_bc_patch_hdf
    !****************************************************************************************












    !>  Return the bc_patch connectivity information for a boundary condition face group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_bc_patch_hdf(bcface_id) result(bc_patch)
        use iso_c_binding,  only: c_ptr, c_loc
        integer(HID_T), intent(in)  :: bcface_id

        integer(HID_T)      :: faces_did, faces_sid
        integer(HSIZE_T)    :: dims(2), maxdims(2)
        integer(ik)         :: nbcfaces, npts_face, ierr

        integer(ik), allocatable, target    :: bc_patch(:,:)
        type(c_ptr)                         :: bc_patch_p
        

        ! Open Faces patch data
        ! TODO: WARNING, should replace with XI_MIN, XI_MAX, etc. somehow. Maybe not...
        call h5dopen_f(bcface_id, "Faces", faces_did, ierr, H5P_DEFAULT_F)


        !  Get the dataspace id and dimensions
        call h5dget_space_f(faces_did, faces_sid, ierr)
        call h5sget_simple_extent_dims_f(faces_sid, dims, maxdims, ierr)
        nbcfaces  = dims(1)
        npts_face = dims(2)


        ! Read boundary condition patch connectivity
        allocate(bc_patch(nbcfaces,npts_face),stat=ierr)
        if (ierr /= 0) call AllocationError
        bc_patch_p = c_loc(bc_patch(1,1))
        call h5dread_f(faces_did, H5T_NATIVE_INTEGER, bc_patch_p, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_bc_patch_hdf: h5dread_f")


        call h5dclose_f(faces_did,ierr)
        call h5sclose_f(faces_sid,ierr)

    end function get_bc_patch_hdf
    !***************************************************************************************







    

    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine set_bc_patch_group_hdf(patch_id,group)
        integer(HID_T), intent(in)  :: patch_id
        character(*),   intent(in)  :: group

        integer(ik) :: ierr

        ! Set 'Boundary State Group'
        call h5ltset_attribute_string_f(patch_id, ".", "Boundary State Group", trim(group), ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_bc_patch_group_hdf: error setting the attribute 'Boundary State Group'")


    end subroutine set_bc_patch_group_hdf
    !***************************************************************************************

    
    !>  Return 'Boundary State Group' attribute for a given patch.
    !!
    !!
    !!  If found, returns the group attribute.
    !!  If not found, returns 'empty'.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    function get_bc_patch_group_hdf(patch_id) result(group_trim)
        integer(HID_T), intent(in)  :: patch_id

        character(1024)             :: group
        character(:),   allocatable :: group_trim
        integer(ik)                 :: ierr
        integer(ik)                 :: exists


        call check_attribute_exists_hdf(patch_id,"Boundary State Group","Soft Fail",exists)

        ! Get 'Boundary State Group'
        if (exists==0) then
            call h5ltget_attribute_string_f(patch_id, ".", "Boundary State Group", group, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_bc_patch_group_hdf: error setting the attribute 'Boundary State Group'")
            group_trim = trim(group)
        else
            group_trim = 'empty'
        end if
            
    end function get_bc_patch_group_hdf
    !***************************************************************************************





    !>  Add a bc_state group to the ChiDG HDF file.
    !!
    !!  A bc_state group is a group of bc_state functions. An example might be a fluid inlet.
    !!  The fluid inlet group might be composed of the Total Inlet and Spalart-Allmaras Inlet
    !!  bc_state functions. The Total Inlet bc_state defines the boundary state for
    !!  the Euler/Navier-Stokes equations. The Spalart-Allmaras Inlet defines the boundary state
    !!  for the extra PDE defined in the Spalart-Allmaras turbulence model.
    !!  
    !!      
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid             HDF file identifier
    !!  @param[in]  group_name      Unique name for the new boundary condition state group.
    !!
    !----------------------------------------------------------------------------------------
    subroutine create_bc_group_hdf(fid,group_name)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: group_name

        character(:),   allocatable :: user_msg
        integer(HID_T)              :: bcgroup_id
        integer(ik)                 :: ierr
        logical                     :: group_exists


        ! Check if bc_state group exists
        group_exists = check_link_exists_hdf(fid,"BCSG_"//trim(group_name))

        user_msg = "create_bc_state_group_hdf: Boundary condition state group already exists. &
                    Cannot have two groups with the same name"
        if (group_exists) call chidg_signal_one(FATAL,user_msg,trim(group_name))


        !
        ! Create a new group for the bc_state_t
        !
        call h5gcreate_f(fid, "BCSG_"//trim(group_name), bcgroup_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'create_bc_state_group_hdf: error creating new group for bc_state.')


        ! Set 'Family'
        call h5ltset_attribute_string_f(bcgroup_id, '.', 'Family', 'none', ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"create_bc_state_group_hdf: error setting the attribute 'Family'")


        call h5gclose_f(bcgroup_id,ierr)

    end subroutine create_bc_group_hdf
    !****************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_bc_state_group_family_hdf(bcgroup_id) result(family_trimmed)
        integer(HID_T),     intent(in)  :: bcgroup_id

        integer(ik)                 :: ierr
        character(1024)             :: family
        character(:),   allocatable :: family_trimmed

        call h5ltget_attribute_string_f(bcgroup_id, ".", "Family", family, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_bc_state_group_family_hdf: h5ltget_attribute_int_f")
        
        family_trimmed = trim(family)

    end function get_bc_state_group_family_hdf
    !*****************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine set_bc_state_group_family_hdf(bcgroup_id,group_family)
        integer(HID_T),     intent(in)  :: bcgroup_id
        character(*),       intent(in)  :: group_family

        integer(ik)                 :: ierr


        call h5ltset_attribute_string_f(bcgroup_id, ".", "Family", trim(adjustl(group_family)), ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_bc_state_group_family_hdf: h5ltget_attribute_int_f")
        

    end subroutine set_bc_state_group_family_hdf
    !*****************************************************************************************









    !>  Return the number of boundary condition state groups are in the HDF file.
    !!
    !!  Boundary condition state groups: 'BCSG_'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_nbc_state_groups_hdf(fid) result(ngroups)
        integer(HID_T)  :: fid
        
        integer(ik)     :: ngroups
        type(svector_t) :: bc_state_group_names


        bc_state_group_names = get_bc_state_group_names_hdf(fid)

        ngroups = bc_state_group_names%size()

    end function get_nbc_state_groups_hdf
    !****************************************************************************************





    
    !>  Return a vector of the names for each boundary condition state group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_bc_state_group_names_hdf(fid) result(bc_state_group_names)
        integer(HID_T), intent(in)  :: fid

        integer(ik)     :: nmembers, ierr, igrp, type
        character(1024) :: gname
        type(svector_t) :: bc_state_group_names


        !  Get number of groups linked to the current bc_face
        call h5gn_members_f(fid, ".", nmembers, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_bc_state_group_names_hdf: error h5gn_members_f")

        if ( nmembers > 0 ) then
            do igrp = 0,nmembers-1
                ! Get group name
                call h5gget_obj_info_idx_f(fid, ".", igrp, gname, type, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_bc_state_group_names_hdf: error h5gget_obj_info_idx_f")

                ! Test if group is a boundary condition state. 'BCSG_'
                if (gname(1:5) == 'BCSG_') then
                    call bc_state_group_names%push_back(string_t(trim(gname(6:))))
                end if
            end do  ! igrp
        end if

    end function get_bc_state_group_names_hdf
    !***************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/11/2016
    !!
    !---------------------------------------------------------------------------------------
    function open_bc_group_hdf(fid,group_name) result(bcgroup_id)
        integer(HID_T),     intent(in)  :: fid
        character(*),       intent(in)  :: group_name

        integer(ik)     :: ierr
        integer(HID_T)  :: bcgroup_id

        call h5gopen_f(fid,"BCSG_"//trim(group_name),bcgroup_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"open_bc_group_hdf: Error opening boundary condition group")

    end function open_bc_group_hdf
    !***************************************************************************************


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/11/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine close_bc_group_hdf(bcgroup_id)
        integer(HID_T), intent(in)  :: bcgroup_id

        integer(ik)     :: ierr

        call h5gclose_f(bcgroup_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"close_bc_group_hdf: Error closing bc_group")

    end subroutine close_bc_group_hdf
    !***************************************************************************************



    !>  Add bc_state function to a group of bc_states.
    !!
    !!  /BCSG_name/BCS_bc_state
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!  @date   11/8/2016   moved from domain to groups
    !!
    !----------------------------------------------------------------------------------------
    subroutine add_bc_state_hdf(bcgroup_id,bc_state)
        integer(HID_T),     intent(in)  :: bcgroup_id
        class(bc_state_t),  intent(in)  :: bc_state

        integer(ik)                         :: ierr
        integer(HID_T)                      :: state_id
        character(:),   allocatable         :: current_family, user_msg
        logical                             :: link_exists, state_found


        if ( bc_state%get_name() == 'empty' ) then
            !
            ! If 'empty' do not allocate new bc
            !
            
        else

            ! Check to make sure the bc_state wasn't previously added
            link_exists = check_link_exists_hdf(bcgroup_id,"BCS_"//bc_state%get_name())


            if (.not. link_exists) then

                ! Check bc_state exists in the register. 
                ! If not, user probably entered the wrong string, so do nothing
                state_found = check_bc_state_registered(bc_state%get_name())

                if (state_found) then

                    ! Get bcgroup family
                    current_family = get_bc_state_group_family_hdf(bcgroup_id)

                    !
                    ! Check if new bc_state is of same family
                    !
                    if ( (trim(current_family) == 'none') .or. &
                         (trim(current_family) == trim(bc_state%get_family())) ) then

                        ! Set group 'Family'
                        call set_bc_state_group_family_hdf(bcgroup_id, bc_state%get_family())

                        ! Create a new group for the bc_state_t
                        call h5gcreate_f(bcgroup_id, "BCS_"//bc_state%get_name(), state_id, ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,"add_bc_state_hdf: error creating new group for bc_state")

                        ! Add bc_state properties to the group that was created
                        call add_bc_properties_hdf(state_id,bc_state)

                        ! Close function group
                        call h5gclose_f(state_id,ierr)

                    else
                        user_msg = "add_bc_state_hdf: Boundary condition state functions in a group &
                                    must be of the same family"
                        call chidg_signal_one(FATAL,user_msg,bc_state%get_family())
                    end if
                end if

            end if

        end if


    end subroutine add_bc_state_hdf
    !*****************************************************************************************











    !>  Return the number of bc_state's attached to a boundary condition face group.
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   10/17/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_nbc_states_hdf(bcgroup_id) result(nbc_states)
        integer(HID_T), intent(in)  :: bcgroup_id

        type(svector_t) :: bc_states
        integer(ik)     :: nbc_states

        bc_states = get_bc_state_names_hdf(bcgroup_id)
        nbc_states = bc_states%size()

    end function get_nbc_states_hdf
    !****************************************************************************************







    

    !>  Return the names of bc_state's attached to a boundary condition state group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    function get_bc_state_names_hdf(bcgroup_id) result(bc_state_names)
        integer(HID_T), intent(in)  :: bcgroup_id

        integer(ik)     :: nmembers, ierr, igrp, type
        character(1024) :: gname
        type(svector_t) :: bc_state_names


        !  Get number of groups linked to the current bc_face
        call h5gn_members_f(bcgroup_id, ".", nmembers, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_bc_states_names_hdf: error h5gn_members_f")

        if ( nmembers > 0 ) then
            do igrp = 0,nmembers-1
                ! Get group name
                call h5gget_obj_info_idx_f(bcgroup_id, ".", igrp, gname, type, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_bc_state_names_hdf: error h5gget_obj_info_idx_f")

                ! Test if group is a boundary condition state. 'BCS_'
                if (gname(1:4) == 'BCS_') then
                    call bc_state_names%push_back(string_t(trim(gname(5:))))
                end if
            end do  ! igrp
        end if

    end function get_bc_state_names_hdf
    !****************************************************************************************







    

    !>  Given the name of a bc_state on a face, return an initialized bc_state instance.
    !!
    !!  You may consider calling 'get_bc_state_names_hdf' first to get a list of 
    !!  available bc_state's on a face. Then the names could be passed into this routine
    !!  to return the bc_state instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_bc_state_hdf(bcgroup_id,bcstate_name) result(bc_state)
        integer(HID_T), intent(in)  :: bcgroup_id
        character(*),   intent(in)  :: bcstate_name


        class(bc_state_t),  allocatable :: bc_state 
        character(:),       allocatable :: bcname, pname, oname
        character(1024)                 :: fname
        integer(HID_T)                  :: bcstate_id, bcprop_id
        integer(ik)                     :: ierr, iprop, nprop, iopt, noptions
        real(rdouble), dimension(1)     :: buf
        real(rk)                        :: ovalue


        ! Open bc_state group
        call h5gopen_f(bcgroup_id, "BCS_"//trim(bcstate_name), bcstate_id, ierr)
        if (ierr /= 0) call chidg_signal_one(FATAL,"get_bc_state_hdf: error opening bc_state group.",trim(bcstate_name))

        
        ! Get boundary condition name string
        if (bcstate_name(1:4) == "BCS_") then
            bcname = trim(bcstate_name(5:))
        else
            bcname = trim(bcstate_name)
        end if


        ! Create boundary condition state and get number of properties
        call create_bc(bcname,bc_state)
        nprop = bc_state%get_nproperties()

        
        ! Loop through properties
        do iprop = 1,nprop


            ! Get property name + open HDF group
            pname = bc_state%get_property_name(iprop)
            call h5gopen_f(bcstate_id, "BCP_"//trim(pname), bcprop_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_bc_state_hdf: error opening bcproperty group.")


            ! Read the function name set for the property.
            call h5ltget_attribute_string_f(bcprop_id, ".", "Function", fname, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_bc_state_hdf: error getting function name.")

            
            ! Set/Create the function for the current property
            call bc_state%set_fcn(trim(pname), trim(fname))

            
            ! Get number of options for the function
            noptions = bc_state%get_noptions(iprop)



            ! Get each option value
            do iopt = 1,noptions
                ! Get option name
                oname = bc_state%get_option_key(iprop,iopt)

                ! Get option value from file
                call h5ltget_attribute_double_f(bcprop_id, ".", trim(oname), buf, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_bc_state_hdf: error getting option value")
                ovalue = real(buf(1),rk)

                ! Set boundary condition option
                call bc_state%set_fcn_option(trim(pname), trim(oname), ovalue)
            end do ! iopt



            ! Close current property group
            call h5gclose_f(bcprop_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_bc_state_hdf: h5gclose")



        end do !iprop



        ! Close boundary condition state group
        call h5gclose_f(bcstate_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_bc_state_hdf: h5gclose")


    end function get_bc_state_hdf
    !*****************************************************************************************










    !>  Add properties to a bc_state on a boundary for a particular domain.
    !!
    !!  /D_domainname/BoundaryConditions/"face"/BCS_bc_state/BCP_bc_property
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !!  @param[in]      bcstate_id      HDF identifier of the bc_state group in the file to 
    !!                                  be modified.
    !!  @param[inout]   bc_state        bc_state class that can be queried for properties to 
    !!                                  be set.
    !!
    !-----------------------------------------------------------------------------------------
    subroutine add_bc_properties_hdf(bcstate_id, bc_state)
        integer(HID_T),     intent(in)  :: bcstate_id
        class(bc_state_t),  intent(in)  :: bc_state

        integer(HID_T)                  :: prop_id
        integer(HSIZE_T)                :: adim
        integer(ik)                     :: iprop, nprop, iopt, nopt, ierr
        character(len=1024)             :: pstring
        character(len=:),   allocatable :: option_key, fcn_name
        real(rk)                        :: option_value


        !
        ! Get number of functions in the boundary condition
        !
        nprop = bc_state%get_nproperties()


        !
        ! Loop through and add properties
        !
        do iprop = 1,nprop

            ! Get string the property is associated with
            pstring = bc_state%get_property_name(iprop)

            ! Create a new group for the property
            call h5gcreate_f(bcstate_id, "BCP_"//trim(adjustl(pstring)), prop_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"add_bcproperties_hdf: error creating new group for bcfunction")

            ! Set property function attribute
            call set_bc_property_function_hdf(prop_id, bc_state%bcproperties%bcprop(iprop)%fcn)


            !
            ! Get number of options available for the current property
            !
            nopt = bc_state%get_noptions(iprop)

            if (nopt > 0 ) then
                do iopt = 1,nopt

                    ! Get the current option and default value.
                    option_key   = bc_state%get_option_key(iprop,iopt)
                    option_value = bc_state%get_option_value(iprop,option_key)

                    ! Set the option as a real attribute
                    adim = 1
                    call h5ltset_attribute_double_f(prop_id, ".", option_key, [real(option_value,rdouble)], adim, ierr)

                end do
            end if


            ! Close function group
            call h5gclose_f(prop_id,ierr)

        end do !ifcn


    end subroutine add_bc_properties_hdf
    !******************************************************************************************









    !>  Set a function for a bc_state property.
    !!
    !!  /D_domainname/BoundaryConditions/"face"/BCS_bcstatename/BCP_bcpropertyname/
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !--------------------------------------------------------------------------------------------
    subroutine set_bc_property_function_hdf(bcprop_id, fcn)
        integer(HID_T),     intent(in)  :: bcprop_id
        class(function_t),  intent(in)  :: fcn

        integer(HSIZE_T)                :: adim
        character(len=:),   allocatable :: option
        real(rk)                        :: val
        integer(ik)                     :: nopt, iopt
        integer                         :: ierr
        

        !
        ! Delete bcproperty attributes
        !
        call delete_group_attributes_hdf(bcprop_id)


        !
        ! Set 'Function' attribute
        !
        call h5ltset_attribute_string_f(bcprop_id, ".", "Function", trim(fcn%get_name()), ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_bc_property_function_hdf: error setting function name")


        !
        ! Set function options
        !
        nopt = fcn%get_noptions()

        do iopt = 1,nopt

            option = fcn%get_option_key(iopt)
            val    = fcn%get_option_value(option)
            !
            ! Set option
            !
            adim = 1
            call h5ltset_attribute_double_f(bcprop_id, ".", trim(option), [real(val,rdouble)], adim, ierr)

        end do ! iopt


    end subroutine set_bc_property_function_hdf
    !***********************************************************************************************







    !>  Remove a bc_state group from the HDF file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine remove_bc_state_group_hdf(fid,group_name)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: group_name

        integer(HID_T)  :: bcgroup_id
        integer(ik)     :: istate, ierr
        logical         :: group_exists
        type(svector_t) :: bc_state_names
        type(string_t)  :: state_name


        group_exists = check_link_exists_hdf(fid,"BCSG_"//trim(group_name))

        if (group_exists) then

            ! Open boundary condition state group
            call h5gopen_f(fid,"BCSG_"//trim(group_name),bcgroup_id,ierr)

            ! Get the names of all bc_states in the group
            bc_state_names = get_bc_state_names_hdf(bcgroup_id)

            ! Call remove for each one
            do istate = 1,bc_state_names%size()
                state_name = bc_state_names%at(istate)
                call remove_bc_state_hdf(bcgroup_id,state_name%get())
            end do

            ! Close the bc_state group
            call h5gclose_f(bcgroup_id,ierr)

            ! Unlink the bc_state group
            call h5gunlink_f(fid,"BCSG_"//trim(group_name),ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"remove_bc_state_group_hdf: error unlinking bc_state group")
        end if

    end subroutine remove_bc_state_group_hdf
    !***********************************************************************************************










    !>  Remove a bc_state from the HDF file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!  @note   Modified to include bc_states
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine remove_bc_state_hdf(bcgroup_id,state_string)
        integer(HID_T),     intent(in)  :: bcgroup_id
        character(len=*),   intent(in)  :: state_string

        integer(HID_T)                          :: bc_state

        integer(HSIZE_T)                        :: iattr, idx
        integer(ik)                             :: nattr
        integer                                 :: nmembers, igrp, type, ierr, iter, iprop, nprop
        character(len=10)                       :: faces(NFACES)
        character(len=1024),    allocatable     :: anames(:), pnames(:)
        character(len=1024)                     :: gname
        type(h5o_info_t), target                :: h5_info



        ! Open the bc_state group
        call h5gopen_f(bcgroup_id, "BCS_"//trim(state_string), bc_state, ierr)


        ! Delete overall boundary condition face attributes
        call delete_group_attributes_hdf(bc_state)


        !  Get number of groups linked to the current bc_state
        call h5gn_members_f(bc_state, ".", nmembers, ierr)


        !
        !  Loop through groups and delete properties
        !
        if ( nmembers > 0 ) then

            !
            ! First get number of states. This could be different than number of groups.
            !
            nprop = 0
            do igrp = 0,nmembers-1

                ! Get group name
                call h5gget_obj_info_idx_f(bc_state, ".", igrp, gname, type, ierr)

                ! Test if group is a boundary condition function. 'BCP_'
                if (gname(1:4) == 'BCP_') then
                    ! increment nprop
                    nprop = nprop + 1
                end if

            end do  ! igrp


            !
            ! Second, get all state names
            !
            allocate(pnames(nprop), stat=ierr)
            if (ierr /= 0) call AllocationError
            iprop = 1
            do igrp = 0,nmembers-1

                ! Get group name
                call h5gget_obj_info_idx_f(bc_state, ".", igrp, gname, type, ierr)

                ! Test if group is a boundary condition function. 'BCP_'
                if (gname(1:4) == 'BCP_') then
                    ! Store name
                    pnames(iprop) = gname
                    iprop = iprop + 1
                end if

            end do ! igrp



            !
            ! Now, go about deleting them all.
            ! Previously, we were deleting them one at a time, but then the index
            ! traversal call get_obj_info_idx was failing for more than one property
            ! because the index was screwed up.
            !
            do iprop = 1,nprop
                call remove_bc_property_hdf(bc_state,pnames(iprop))
            end do


        end if ! nmembers


        !
        ! Close the bc_state group
        !
        call h5gclose_f(bc_state,ierr)

        !
        ! Unlink the bc_state group
        !
        call h5gunlink_f(bcgroup_id,"BCS_"//trim(state_string),ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"delete_bc_state_hdf: error unlinking bc_state group")



        !
        ! If no bc_state's are left attached, clear group family.
        !
        if (get_nbc_states_hdf(bcgroup_id) == 0) then
            call set_bc_state_group_family_hdf(bcgroup_id,'none')
        end if



    end subroutine remove_bc_state_hdf
    !****************************************************************************************************











    !>  Remove the properties from a bc_state group in the HDF file
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!  @note   Modified to include bc_states
    !!
    !!
    !------------------------------------------------------------------------------------------------------
    subroutine remove_bc_property_hdf(bc_state,pname)
        integer(HID_T),     intent(in)      :: bc_state
        character(*),       intent(in)      :: pname

        integer(HID_T)  :: bcprop
        integer(ik)     :: ierr

        !
        ! Open bcproperty group
        !
        call h5gopen_f(bc_state, trim(adjustl(pname)), bcprop, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"delete_bc_property_hdf: error opening bcproperty group")


        !
        ! Delete bcproperty attributes
        !
        call delete_group_attributes_hdf(bcprop)


        !
        ! Close bcproperty group
        !
        call h5gclose_f(bcprop, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"delete_bc_property_hdf: error closing bcproperty group")


        !
        ! Now that the data in bcproperty has been removed, unlink the bcproperty group.
        !
        call h5gunlink_f(bc_state,trim(pname),ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"delete_bcfunction_hdf: error unlinking bcproperty group")


    end subroutine remove_bc_property_hdf
    !******************************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/15/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    function check_bc_state_exists_hdf(bcface_id,bc_state) result(exist_status)
        integer(HID_T),     intent(in)  :: bcface_id
        character(len=*),   intent(in)  :: bc_state

        integer(ik) :: ierr
        logical     :: exist_status

        ! Check if face contains the bc_state
        call h5lexists_f(bcface_id, "BCS_"//trim(bc_state), exist_status, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"check_bc_state_exists: Error in call to h5lexists_f")


    end function check_bc_state_exists_hdf
    !*********************************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!  @note   Modified to include bc_states
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function check_bc_property_exists_hdf(bcface_id,pname) result(exist_status)
        integer(HID_T),     intent(in)  :: bcface_id
        character(*),       intent(in)  :: pname

        integer(HID_T)          :: bc_state
        integer                 :: ierr, nmembers, igrp, type, iop
        character(len=1024)     :: gname
        logical                 :: exist_status
        type(svector_t)         :: bc_state_strings
        type(string_t)          :: string

        
        !
        !  Loop through groups and detect bc_state's that could contain property
        !
        call h5gn_members_f(bcface_id, ".", nmembers, ierr)
        if ( nmembers > 0 ) then

            ! First get number of states. This could be different than number of groups.
            do igrp = 0,nmembers-1

                ! Get group name
                call h5gget_obj_info_idx_f(bcface_id, ".", igrp, gname, type, ierr)

                ! Test if group is a boundary condition function. 'BCS_'
                if (gname(1:4) == 'BCS_') then
                    call bc_state_strings%push_back(string_t(trim(gname)))
                end if

            end do  ! igrp

        end if



        !
        ! Find the state with the property
        !
        exist_status = .false.
        do iop = 1,bc_state_strings%size()

            ! Open the state group
            string = bc_state_strings%at(iop)
            call h5gopen_f(bcface_id, string%get(), bc_state, ierr)

            ! Check if it contains a link to the property group
            call h5lexists_f(bc_state, "BCP_"//trim(pname), exist_status, ierr)



            if (exist_status) then
                ! Close state
                call h5gclose_f(bc_state,ierr)
                exit
            end if

            ! Close state
            call h5gclose_f(bc_state,ierr)

        end do !iop



    end function check_bc_property_exists_hdf
    !********************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !--------------------------------------------------------------------------------------------
    function check_domain_exists_hdf(fid,domain_name) result(exist_status)
        integer(HID_T),     intent(in)  :: fid
        character(len=*),   intent(in)  :: domain_name

        integer(ik) :: ierr
        logical     :: exist_status

        ! Check if face contains the bc_state
        call h5lexists_f(fid, "D_"//trim(domain_name), exist_status, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"check_domain_exists_hdf: Error in call to h5lexists_f")


    end function check_domain_exists_hdf
    !*********************************************************************************************








    !>  Delete all the attributes attached to a specified group identifier.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @param[in]  gid     HDF5 group identifier.
    !!
    !---------------------------------------------------------------------------------------
    subroutine delete_group_attributes_hdf(gid)
        integer(HID_T),         intent(in)      :: gid

        integer(ik)                             :: nattr, ierr
        integer(HSIZE_T)                        :: iattr, idx
        type(h5o_info_t), target                :: h5_info


        !
        ! Get number of attributes attached to the group id
        !
        call h5oget_info_f(gid, h5_info, ierr)
        nattr = h5_info%num_attrs
        if (ierr /= 0) call chidg_signal(FATAL,"delete_group_attributes_hdf: error getting current number of attributes.")


        !
        ! Delete any existing attributes
        !
        if ( nattr > 0 ) then

            !
            ! Delete by index. h5adelete_by_idx_f returns idx with the next index so it doesn't need manually updated.
            !
            idx = 0
            do iattr = 1,nattr
                call h5adelete_by_idx_f(gid, ".", H5_INDEX_CRT_ORDER_F, H5_ITER_NATIVE_F, idx, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"delete_group_attributes_hdf: error deleting attribute")
            end do


        end if ! nattr


    end subroutine delete_group_attributes_hdf
    !****************************************************************************************







    !>  Check if an attribute exists on an HDF group.
    !!
    !!  fail_type
    !!      = 'Soft Fail'   routine exists and returns a failure status to the caller.
    !!      = 'Hard Fail'   routine calls a FATAL exit routine and brings down the program.
    !!
    !!  fail_status
    !!      = 1             Attribute was not found, does not exists.
    !!      = 0             Attribute found, no failure.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/25/2016
    !!
    !----------------------------------------------------------------------------------------
    subroutine check_attribute_exists_hdf(id,attribute,fail_type,fail_status)
        integer(HID_T), intent(in)                  :: id
        character(*),   intent(in)                  :: attribute
        character(*),   intent(in),     optional    :: fail_type
        integer(ik),    intent(inout),  optional    :: fail_status

        logical                         :: attribute_exists, hard_stop
        character(len=:), allocatable   :: msg
        integer(ik)                     :: ierr


        !
        ! Query attribute existence
        !
        call h5aexists_f(id,trim(attribute),attribute_exists,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"check_attribute_exists_hdf: Error checking if attribute exists")

        !
        ! Default error mode
        !
        hard_stop = .true.

        !
        ! Check if failure type was specified
        !
        if (present(fail_type)) then

            if (trim(fail_type) == "Soft Fail") then
                hard_stop = .false.
            else if (trim(fail_type) == "Hard Fail") then
                hard_stop = .true.
            else
                call chidg_signal(FATAL,"check_attribute_exists_hdf: Didn't recognize the specified 'fail_status'.")
            end if

        end if

        !
        ! Handle error if necessary
        !
        msg = "Attribute "//trim(attribute)//" not found in the file. Maybe the file was generated &
               with an old version of the ChiDG library. Try regenerating the HDF grid file with an &
               updated version of the ChiDG library to make sure the file is formatted properly"

        if (.not. attribute_exists) then

            if (present(fail_status)) then
                fail_status = 1     ! Didnt find attribute
            end if

            if (hard_stop) then
                call chidg_signal(FATAL,msg)
            end if

        else
            if (present(fail_status)) then
                fail_status = 0     ! Found attribute, no failure
            end if
        end if


    end subroutine check_attribute_exists_hdf
    !****************************************************************************************





    !>  Check if group or dataset exists on an HDF identifier.
    !!
    !!  This calls h5lexists_f. Groups and Datasets don't actually have names, but, the
    !!  links to those objects in the file do have names. So, we can check if there is 
    !!  a link with the desired name for either a Group or a Dataset.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/15/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    function check_link_exists_hdf(id,linkname) result(exist_status)
        integer(HID_T),     intent(in)  :: id
        character(*),       intent(in)  :: linkname

        integer(ik) :: ierr
        logical     :: exist_status

        ! Check if group exists
        call h5lexists_f(id, trim(linkname), exist_status, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"check_link_exists_hdf: Error in call to h5lexists_f")


    end function check_link_exists_hdf
    !*********************************************************************************************








    !>  Check if a ChiDG file exists.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/23/2016
    !!
    !----------------------------------------------------------------------------------------------
    function check_file_exists_hdf(filename) result(exist_status)
        character(*),   intent(in)  :: filename

        logical :: exist_status


        ! Check file exists
        inquire(file=filename, exist=exist_status)


    end function check_file_exists_hdf
    !*********************************************************************************************




end module mod_hdf_utilities
