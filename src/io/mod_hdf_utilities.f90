module mod_hdf_utilities
#include <messenger.h>
    use mod_kinds,              only: rk, ik, rdouble
    use mod_constants,          only: NFACES, TWO_DIM, THREE_DIM
    use mod_file_utilities,     only: delete_file
    use mod_bc,                 only: check_bc_state_registered
    use mod_string,             only: string_t
    use mod_function,           only: create_function
    use type_function,          only: function_t
    use type_svector,           only: svector_t
    use type_bc_state,          only: bc_state_t
    use type_file_properties,   only: file_properties_t
    use hdf5
    use h5lt
    implicit none

    

    !
    ! HDF5 storage format
    !
    integer, parameter :: STORAGE_FORMAT_MAJOR = 1
    integer, parameter :: STORAGE_FORMAT_MINOR = 0


    ! Attribute sizes
    integer(HSIZE_T), parameter :: SIZE_ONE = 1



contains

    !----------------------------------------------------------------------------------------
    !!
    !!  ChiDG HDF File Format API
    !!
    !!  Procedures:
    !!  -----------
    !!
    !!  initialize_file_hdf
    !!  check_file_storage_version_hdf
    !!
    !!  set_storage_version_major_hdf
    !!  set_storage_version_minor_hdf
    !!  get_storage_version_major_hdf
    !!  get_storage_version_minor_hdf
    !!
    !!  get_properties_hdf
    !!
    !!  set_ndomains_hdf
    !!  get_ndomains_hdf
    !!
    !!  set_contains_grid_hdf
    !!  get_contains_grid_hdf
    !!
    !!  set_contains_solution_hdf
    !!  get_contains_solution_hdf
    !!
    !!  get_domain_name_hdf
    !!  get_domain_names_hdf
    !!
    !!  set_domain_index_hdf
    !!  get_domain_index_hdf
    !!  get_domain_indices_hdf
    !!
    !!  set_domain_coordinates_hdf
    !!  set_domain_elements_hdf
    !!
    !!
    !!  set_coordinate_order_hdf
    !!  get_coordinate_order_hdf
    !!  get_coordinate_orders_hdf
    !!
    !!  set_solution_order_hdf
    !!  get_solution_order_hdf
    !!  get_solution_orders_hdf
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
    !!  set_bc_patch_hdf
    !!  add_bc_state_hdf
    !!  get_nbc_states_hdf
    !!  add_bc_properties_hdf
    !!  remove_bc_state_hdf
    !!  remove_bc_property_hdf
    !!
    !!  get_bcnames_hdf
    !!
    !!  check_bc_state_exists_hdf
    !!  check_bc_property_exists_hdf
    !!  delete_group_attributes_hdf
    !!  check_attribute_exists_hdf
    !!      
    !!
    !****************************************************************************************




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
    function initialize_file_hdf(filename) result(fid)
        character(*),   intent(in)  :: filename

        logical         :: file_exists
        integer(HID_T)  :: fid
        integer(ik)     :: ierr
        

        !
        ! Check if input file already exists
        !
        inquire(file=filename, exist=file_exists)
        if (file_exists) then
            call write_line("Found "//trim(filename)//" that already exists. Deleting it to create new file...")
            call delete_file(trim(filename))
        end if


        !
        ! Create file
        !
        call h5fcreate_f(trim(filename)//".h5", H5F_ACC_TRUNC_F, fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"initialize_file_hdf: Error h5fcreate_f")
        call write_line("File created: "//trim(filename)//".h5")


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


    end function initialize_file_hdf
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

            call write_line(msg)

            read_user_input = .true.
            do while(read_user_input)
                
                read(*,*) user_option

                if (user_option == 1) then
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
            prop%order_s = get_solution_orders_hdf(fid,prop%domain_names)
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

        character(len=1024), allocatable    :: names(:)
        character(len=1024)                 :: gname
        integer                             :: ndomains, nmembers, type
        integer                             :: igrp, idom, ierr

        !
        ! Get number of domains
        !
        ndomains = get_ndomains_hdf(fid)


        !
        ! Allocate names
        !
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
            !
            ! Get group name
            !
            call h5gget_obj_info_idx_f(fid,"/", igrp, gname, type, ierr)

            !
            ! Test if group is a 'Domain'
            !
            if (gname(1:2) == 'D_') then

                !
                ! Store name
                !
                names(idom) = trim(gname)

                idom = idom + 1
            end if
        end do


    end function get_domain_names_hdf
    !****************************************************************************************









    !>  Return a domain name given a domain index.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid         HDF file identifier
    !!  @param[in]  idom_hdf    A specified domain index to be queried. This is an attribute of each domain in 
    !!                          the HDF file, per the ChiDG convention.
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_name_hdf(fid,idom_hdf) result(dname)
        integer(HID_T),     intent(in)  :: fid
        integer(ik),        intent(in)  :: idom_hdf

        character(len=1024)                 :: dname
        character(len=1024), allocatable    :: dnames(:)
        integer(ik),         allocatable    :: dindices(:)
        integer                             :: iind, ndomains

        !
        ! Get number of domains, domain names, and domain indices
        !
        ndomains = get_ndomains_hdf(fid)
        dnames   = get_domain_names_hdf(fid)
        dindices = get_domain_indices_hdf(fid)



        do iind = 1,ndomains

            if ( dindices(iind) == idom_hdf ) then
               dname = dnames(iind) 
            end if

        end do



    end function get_domain_name_hdf
    !****************************************************************************************









    !> Return a list of domain indices from an HDF5 file identifier. This is because, the current method of detecting
    !! domains by name can change the order they are detected in. So, each domain is given an idomain attribute that 
    !! is independent of the order of discovery from the file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_index_hdf(block_id,domain_index)
        integer(HID_T),     intent(in)  :: block_id
        integer(ik),        intent(in)  :: domain_index

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(block_id,".","Domain Index",[domain_index],SIZE_ONE,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_index_hdf: Error h5ltset_attribute_int_f")


    end subroutine set_domain_index_hdf
    !****************************************************************************************





    !> Return a list of domain indices from an HDF5 file identifier. This is because, the current method of detecting
    !! domains by name can change the order they are detected in. So, each domain is given an idomain attribute that 
    !! is independent of the order of discovery from the file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_index_hdf(block_id) result(domain_index)
        integer(HID_T),     intent(in)  :: block_id

        integer(ik) :: domain_index, ierr
        integer, dimension(1) :: buf

        call h5ltget_attribute_int_f(block_id,".","Domain Index",buf,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_index_hdf: Error h5ltget_attribute_int_f")

        domain_index = int(buf(1), kind=ik)

    end function get_domain_index_hdf
    !****************************************************************************************







    !> Return a list of domain indices from an HDF5 file identifier. This is because, 
    !! the current method of detecting domains by name can change the order they are 
    !! detected in. So, each domain is given an idomain attribute that is independent 
    !! of the order of discovery from the file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier
    !!
    !---------------------------------------------------------------------------------------
    function get_domain_indices_hdf(fid) result(indices)
        integer(HID_T),     intent(in)  :: fid

        integer(HID_T)                          :: did
        integer(ik),            allocatable     :: indices(:)
        character(len=1024),    allocatable     :: names(:)
        integer(ik)                             :: idom, ndomains, ierr
        integer, dimension(1)                   :: buf
        integer(HSIZE_T)                        :: adim
        logical                                 :: attribute_exists


        !
        ! Get number of domains
        !
        ndomains = get_ndomains_hdf(fid)
        names    = get_domain_names_hdf(fid)

        !
        ! Allocate indices
        !
        allocate(indices(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        !  Loop through groups and read domain names
        !
        idom = 1
        do idom = 1,ndomains
            !
            ! Open domain group
            !
            call h5gopen_f(fid,trim(adjustl(names(idom))), did, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_indices_hdf: error opening domain group")

            !
            ! Get idomain attribute from fid/domain/idomain
            !
            call h5aexists_f(did, 'Domain Index', attribute_exists, ierr)

            
            !
            ! If it doesn't exist, set to the current value of idom
            !
            adim = 1
            if ( .not. attribute_exists ) then

                ! Set value.
                call h5ltset_attribute_int_f(fid, trim(adjustl(names(idom))), 'Domain Index', [idom], adim, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_domain_indices_hdf: error writing an initial domain index")

            end if


            !
            ! Get value that was just set to be sure. 
            !
            call h5ltget_attribute_int_f(fid, trim(adjustl(names(idom))), 'Domain Index', buf, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_indices_hdf: error retrieving domain indices")

            !
            ! Set value detected to indices array that will be passed back from the function
            !
            indices(idom) = buf(1)


            !
            ! Close domain
            !
            call h5gclose_f(did,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_indices_hdf: h5gclose")

        end do ! idom


    end function get_domain_indices_hdf
    !****************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  block_id        HDF file identifier of a block-domain group
    !!  @param[in]  domain_mapping  Integer specifying the block-domain mapping 
    !!                              1-linear, 2-quadratic, etc.
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_mapping_hdf(block_id,domain_mapping)
        integer(HID_T),     intent(in)  :: block_id
        integer(ik),        intent(in)  :: domain_mapping

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(block_id,".","Domain Mapping",[domain_mapping],SIZE_ONE,ierr)
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
    function get_domain_mapping_hdf(block_id) result(domain_mapping)
        integer(HID_T),     intent(in)  :: block_id

        integer(ik) :: domain_mapping, ierr
        integer, dimension(1) :: buf

        call h5ltget_attribute_int_f(block_id,".","Domain Mapping",buf,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_mapping_hdf: Error h5ltget_attribute_int_f")

        domain_mapping = int(buf(1), kind=ik)

    end function get_domain_mapping_hdf
    !****************************************************************************************










    !>  For a block-domain, set the domain coordinates
    !!
    !!  /D_domainname/Grid/Coordinates
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   10/15/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_coordinates_hdf(block_id,xcoords,ycoords,zcoords)
        integer(HID_T), intent(in)  :: block_id
        real(rk),       intent(in)  :: xcoords(:,:,:)
        real(rk),       intent(in)  :: ycoords(:,:,:)
        real(rk),       intent(in)  :: zcoords(:,:,:)

        integer(HID_T)      :: grid_id, xspace_id, yspace_id, zspace_id, xset_id, yset_id, zset_id
        integer(HSIZE_T)    :: dims_rank_one(1), dims_rank_two(2)
        integer(ik)         :: ierr, ipt, ipt_i, ipt_j, ipt_k, npts, npt_i, npt_j, npt_k
        real(rk), allocatable, dimension(:) :: xcoords_linear, ycoords_linear, zcoords_linear

        !
        ! Create a grid-group within the current block domain
        !
        call h5gcreate_f(block_id, "Grid", grid_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_hdf: h5gcreate_f")




        !
        ! Re-order coordinates to be linear arrays
        !
        npt_i = size(xcoords,1)
        npt_j = size(xcoords,2)
        npt_k = size(xcoords,3)
        npts = npt_i*npt_j*npt_k
        dims_rank_one = npts
        allocate(xcoords_linear(npts), ycoords_linear(npts), zcoords_linear(npts), stat=ierr)

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
        call h5dcreate_f(grid_id, "CoordinateX", H5T_NATIVE_DOUBLE, xspace_id, xset_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_hdf: h5dcreate_f")
        call h5dcreate_f(grid_id, "CoordinateY", H5T_NATIVE_DOUBLE, yspace_id, yset_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_hdf: h5dcreate_f")
        call h5dcreate_f(grid_id, "CoordinateZ", H5T_NATIVE_DOUBLE, zspace_id, zset_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_hdf: h5dcreate_f")

        !
        ! Write coordinates to datasets
        !
        call h5dwrite_f(xset_id, H5T_NATIVE_DOUBLE, xcoords_linear, dims_rank_one, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5dwrite_f")
        call h5dwrite_f(yset_id, H5T_NATIVE_DOUBLE, ycoords_linear, dims_rank_one, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5dwrite_f")
        call h5dwrite_f(zset_id, H5T_NATIVE_DOUBLE, zcoords_linear, dims_rank_one, ierr)
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
    subroutine set_domain_elements_hdf(block_id,elements)
        integer(HID_T), intent(in)  :: block_id
        integer(ik),    intent(in)  :: elements(:,:)

        integer(ik)         :: ierr
        integer(HID_T)      :: element_set_id, element_space_id, grid_id
        integer(HSIZE_T)    :: dims_rank_two(2)

        !
        ! Size element connectivities
        !
        dims_rank_two(1) = size(elements,1)
        dims_rank_two(2) = size(elements,2)


        !
        ! Create dataspace for element connectivity
        !
        call h5screate_simple_f(2, dims_rank_two, element_space_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_elements_hdf: h5screate_simple_f")


        !
        ! Create a grid-group within the current block domain
        !
        call h5gopen_f(block_id, "Grid", grid_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_elements_hdf: h5gcreate_f")


        !
        ! Create dataset for element connectivity
        !
        call h5dcreate_f(grid_id, "Elements", H5T_NATIVE_INTEGER, element_space_id, element_set_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_elements_hdf: h5dcreate_f")



        !
        ! Write element connectivities
        !
        call h5dwrite_f(element_set_id, H5T_NATIVE_INTEGER, elements, dims_rank_two, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_elements_hdf: h5dwrite_f")


        !
        ! Close dataset, dataspace, Grid group
        !
        call h5dclose_f(element_set_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_elements_hdf: h5dclose_f")
        call h5sclose_f(element_space_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_elements_hdf: h5sclose_f")
        call h5gclose_f(grid_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_elements_hdf: h5gclose_f")



    end subroutine set_domain_elements_hdf
    !****************************************************************************************










    !>  Set boundary condition patch face indices for a block boundary.
    !!
    !!  /D_domainname/BoundaryConditions/"face"/Faces
    !!
    !!  "face" is XI_MIN, XI_MAX, ETA_MIN, etc.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/15/2016
    !!
    !!  @param[in]  block_id    HDF identifier for the block to be written to
    !!  @param[in]  faces       Face indices to be set for the boundary condition patch
    !!  @param[in]  bcface      Integer specifying which boundary of the block to write to
    !!
    !---------------------------------------------------------------------------------------
    subroutine set_bc_patch_hdf(block_id,faces,bcface)
        integer(HID_T), intent(in)  :: block_id
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
        call h5gopen_f(block_id, "BoundaryConditions", bc_id, ierr)
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














    !>  Add bc_state function for a block boundary condition.
    !!
    !!  /D_domainname/BoundaryConditions/"face"/BCS_bc_state
    !!
    !!  "face" is XI_MIN, XI_MAX, ETA_MIN, etc.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine add_bc_state_hdf(bcface_id,bc_state)
        integer(HID_T),     intent(in)      :: bcface_id
        class(bc_state_t),  intent(inout)   :: bc_state

        integer(ik)                         :: ierr
        integer(HID_T)                      :: state_id
        logical                             :: link_exists, state_found


        if ( bc_state%get_name() == 'empty' ) then
            !
            ! If 'empty' do not allocate new bc
            !
            
        else

            ! Check to make sure the bc_state wasn't previously added
            call h5lexists_f(bcface_id, "BCS_"//bc_state%get_name(), link_exists, ierr)


            if (.not. link_exists) then

                ! Check bc_state exists in the register. 
                ! If not, user probably entered the wrong string, so do nothing
                state_found = check_bc_state_registered(bc_state%get_name())

                if (state_found) then
                    ! Create a new group for the bc_state_t
                    call h5gcreate_f(bcface_id, "BCS_"//bc_state%get_name(), state_id, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"add_domain_bc_state_hdf: error creating new group for bc_state")

                    ! Add bc_state properties to the group that was created
                    call add_bc_properties_hdf(state_id,bc_state)

                    ! Close function group
                    call h5gclose_f(state_id,ierr)
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
    function get_nbc_states_hdf(bcface_id) result(nbc_states)
        integer(HID_T), intent(in)  :: bcface_id

        integer(ik)     :: nbc_states, nmembers, igrp, type, ierr
        character(1024) :: gname

        nbc_states = 0

        ! Loop over groups and detect bc_state's; Increment nbc_states.
        call h5gn_members_f(bcface_id, ".", nmembers, ierr)
        if ( nmembers > 0 ) then
            do igrp = 0,nmembers-1

                ! Get group name
                call h5gget_obj_info_idx_f(bcface_id, ".", igrp, gname, type, ierr)

                ! Test if group is a boundary condition state. 'BCS_'
                if (gname(1:4) == 'BCS_') then
                    nbc_states = nbc_states + 1
                end if

            end do  ! igrp
        end if


    end function get_nbc_states_hdf
    !****************************************************************************************







    

    !>  Return the names of bc_state's attached to a boundary conditions face group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    function get_bc_state_names_hdf(bcface_id) result(bc_state_names)
        integer(HID_T), intent(in)  :: bcface_id

        integer(ik)     :: nmembers, ierr, igrp, type
        character(1024) :: gname
        type(svector_t) :: bc_state_names



        !  Get number of groups linked to the current bc_face
        call h5gn_members_f(bcface_id, ".", nmembers, ierr)
        if ( nmembers > 0 ) then
            do igrp = 0,nmembers-1
                ! Get group name
                call h5gget_obj_info_idx_f(bcface_id, ".", igrp, gname, type, ierr)

                ! Test if group is a boundary condition state. 'BCS_'
                if (gname(1:4) == 'BCS_') then
                    call bc_state_names%push_back(string_t(trim(gname)))
                end if
            end do  ! igrp
        end if

    end function get_bc_state_names_hdf
    !****************************************************************************************










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
        integer(HID_T),     intent(in)      :: bcstate_id
        class(bc_state_t),  intent(inout)   :: bc_state

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

            !
            ! Get string the property is associated with
            !
            pstring = bc_state%get_property_name(iprop)


            !
            ! Create a new group for the property
            !
            call h5gcreate_f(bcstate_id, "BCP_"//trim(adjustl(pstring)), prop_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"add_bcproperties_hdf: error creating new group for bcfunction")


            !
            ! Print property function attribute
            !
            fcn_name = bc_state%bcproperties%bcprop(iprop)%fcn%get_name()
            call h5ltset_attribute_string_f(prop_id, ".", "function", fcn_name, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"add_bcproperties_hdf: error setting function attribute")


            !
            ! Get number of options available for the current property
            !
            nopt = bc_state%get_noptions(iprop)


            if (nopt > 0 ) then
                do iopt = 1,nopt

                    !
                    ! Get the current option and default value.
                    !
                    option_key   = bc_state%get_option_key(iprop,iopt)
                    option_value = bc_state%get_option_value(iprop,option_key)

                    !
                    ! Set the option as a real attribute
                    !
                    adim = 1
                    call h5ltset_attribute_double_f(prop_id, ".", option_key, [real(option_value,rdouble)], adim, ierr)

                end do
            end if


            !
            ! Close function group
            !
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
    subroutine set_bc_property_function_hdf(bcprop_id, func)
        integer(HID_T),     intent(in)  :: bcprop_id
        class(function_t),  intent(in)  :: func

        integer(HSIZE_T)                :: adim
        class(function_t),  allocatable :: fcn
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
        call h5ltset_attribute_string_f(bcprop_id, ".", "Function", trim(func%get_name()), ierr)
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
    !********************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!  @note   Modified to include bc_states
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine remove_bc_state_hdf(bcface_id,state_string)
        integer(HID_T),     intent(in)  :: bcface_id
        character(len=*),   intent(in)  :: state_string

        integer(HID_T)                          :: bc_state

        integer(HSIZE_T)                        :: iattr, idx
        integer(ik)                             :: nattr
        integer                                 :: nmembers, igrp, type, ierr, iter, iprop, nprop
        character(len=10)                       :: faces(NFACES)
        character(len=1024),    allocatable     :: anames(:), pnames(:)
        character(len=1024)                     :: gname
        type(h5o_info_t), target                :: h5_info



        !
        ! Open the bc_state group
        !
        call h5gopen_f(bcface_id, "BCS_"//trim(state_string), bc_state, ierr)


        !
        ! Delete overall boundary condition face attributes
        !
        call delete_group_attributes_hdf(bc_state)


        !
        !  Get number of groups linked to the current bc_state
        !
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
        call h5gunlink_f(bcface_id,"BCS_"//trim(state_string),ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"delete_bc_state_hdf: error unlinking bc_state group")


    end subroutine remove_bc_state_hdf
    !****************************************************************************************************











    !>
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
    !!  @author Nathan A. Wukie
    !!  @date   9/25/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_coordinate_order_hdf(block_id,order)
        integer(HID_T), intent(in)  :: block_id
        integer(ik),    intent(in)  :: order

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(block_id,".","Coordinate Order", [order],SIZE_ONE,ierr)
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
    subroutine get_coordinate_order_hdf(block_id,order)
        integer(HID_T), intent(in)  :: block_id
        integer(ik),    intent(in)  :: order

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(block_id,".","Coordinate Order", [order],SIZE_ONE,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_coordinate_order_hdf: Error getting 'Coordinate Order' attribute")

    end subroutine get_coordinate_order_hdf
    !****************************************************************************************








    !> Returns an array of integers that specifies the order of the coordinate expansion for every domain.
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
            call h5ltget_attribute_int_f(fid, trim(dnames(idom)), "Domain Mapping", buf, ierr)
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
    subroutine set_solution_order_hdf(block_id,order)
        integer(HID_T), intent(in)  :: block_id
        integer(ik),    intent(in)  :: order

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(block_id,".","Solution Order", [order],SIZE_ONE,ierr)
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
    function get_solution_order_hdf(block_id) result(order)
        integer(HID_T), intent(in)  :: block_id

        integer(ik) :: order, ierr
        integer, dimension(1) :: buf

        call h5ltget_attribute_int_f(block_id,".","Solution Order",buf,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_solution_order_hdf: Error getting 'Solution Order' attribute")

        order = int(buf(1), kind=ik)

    end function get_solution_order_hdf
    !****************************************************************************************










    !> Returns an array of integer that specifies the order of the solution expansion for every domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier.
    !!  @param[in]  dnames  List of domain names to be interrogated.
    !!
    !----------------------------------------------------------------------------------------
    function get_solution_orders_hdf(fid, dnames) result(orders)
        integer(HID_T),         intent(in)  :: fid
        character(len=1024),    intent(in)  :: dnames(:)

        integer(HID_T)              :: did
        integer(HSIZE_T)            :: adim
        logical                     :: order_exists
        integer(ik), allocatable    :: orders(:)
        integer                     :: ierr, idom, order
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
            ! Open domain group
            !
            call h5gopen_f(fid, trim(dnames(idom)), did, ierr)
            if (ierr /= 0) call chidg_signal_one(FATAL,"get_solution_orders_hdf: error opening domain group.", trim(dnames(idom)) )

            
            !
            ! Check 'order_solution' attribute exists
            !
            call h5aexists_f(did, 'Solution Order', order_exists, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_solution_orders_hdf: error check attribue exists.")



            !
            ! Handle attribute does not exist
            !
            if ( .not. order_exists ) then
                adim = 1
                buf = 0
                call h5ltset_attribute_int_f(did, ".", 'Solution Order', buf, adim, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_solution_orders_hdf: error setting attribute.")
            end if



            !
            !  Get solution order
            !
            call h5ltget_attribute_int_f(did, ".", 'Solution Order', buf, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_solution_orders_hdf: error getting 'Solution Order' attribute.")


            !
            ! Compute number of terms in coordinate expansion
            !
            order = buf(1)
            orders(idom) = int(order, kind=ik)


            call h5gclose_f(did,ierr)

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
    subroutine set_domain_dimensionality_hdf(block_id,dimensionality)
        integer(HID_T), intent(in)  :: block_id
        integer(ik),    intent(in)  :: dimensionality

        integer(ik)         :: ierr

        !  Get coordinate mapping
        call h5ltset_attribute_int_f(block_id,".","Domain Dimensionality",[dimensionality],SIZE_ONE,ierr)
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
    function get_domain_dimensionality_hdf(block_id) result(dimensionality)
        integer(HID_T), intent(in)  :: block_id

        integer(ik) :: dimensionality, ierr
        integer, dimension(1)   :: buf

        !  Get coordinate mapping
        call h5ltget_attribute_int_f(block_id,".","Domain Dimensionality",buf,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_dimensionality_hdf: h5ltget_attribute_int_f")

        dimensionality = int(buf(1),kind=ik)

    end function get_domain_dimensionality_hdf
    !****************************************************************************************










    !>  Returns an array of integers that specifies the number of spatial dimensions to use for every domain.
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

            !
            !  Get coordinate mapping
            !
            call h5ltget_attribute_int_f(fid, trim(dnames(idom)), "Domain Dimensionality", dimensionality, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_dimensionalities_hdf: Error h5ltget_attribute_int_f")

            dimensionalities(idom) = dimensionality(1)

        end do

    end function get_domain_dimensionalities_hdf
    !****************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_equation_set_hdf(block_id,equation_set)
        integer(HID_T), intent(in)  :: block_id
        character(*),   intent(in)  :: equation_set

        integer(ik) :: ierr

        !  Get coordinate mapping
        call h5ltset_attribute_string_f(block_id,".","Equation Set",equation_set,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_equation_set_hdf: h5ltset_attribute_string_f")

    end subroutine set_domain_equation_set_hdf
    !****************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_equation_set_hdf(block_id) result(equation_set)
        integer(HID_T), intent(in)  :: block_id

        character(1024) :: equation_set
        integer(ik)     :: ierr

        !  Get coordinate mapping
        call h5ltget_attribute_string_f(block_id,".","Equation Set",equation_set,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_equation_set_hdf: h5ltset_attribute_string_f")

    end function get_domain_equation_set_hdf
    !****************************************************************************************






    !> Returns an array of integer that specifies the order of the solution expansion for every domain.
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
            call h5gopen_f(fid, trim(dnames(idom)), did, ierr)
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










    !>  Return the boundary condition names from the HDF file that are set for each face of a particular domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid         ChiDG HDF5 file identifier.
    !!  @param[in]  dname       String indicating the boundary condition to query.
    !!  @result     bcnames     Array of boundary condition names for the specified domain. bcnames(ibc)
    !!
    !----------------------------------------------------------------------------------------
    function get_bcnames_hdf(fid,dname) result(bcnames)
        integer(HID_T),         intent(in)  :: fid
        character(len=1024),    intent(in)  :: dname

        integer(HID_T)                      :: bc_id, bcface_id
        integer(ik)                         :: ndom, ierr, idom, iface
        integer                             :: igroup, nmembers, type
        type(svector_t),    allocatable     :: bcnames(:)
        character(len=1024)                 :: gname
        character(len=10)                   :: faces(NFACES)
        logical                             :: exists, bcname_exists



        !
        ! Get file information
        !
        ndom   = get_ndomains_hdf(fid)

        allocate(bcnames(NFACES), stat=ierr)
        if (ierr /= 0) call AllocationError




        !
        ! Check if 'BoundaryConditions' group exists
        !
        call h5lexists_f(fid, "/"//trim(dname)//"/BoundaryConditions", exists, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5lexists - checking BoundaryConditions")



        !
        ! Open the Domain/BoundaryConditions group
        !
        if (exists) then
            call h5gopen_f(fid, "/"//trim(dname)//"/BoundaryConditions", bc_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5gopen - BoundaryConditions")
        else
            call h5gcreate_f(fid, '/'//trim(dname)//"/BoundaryConditions", bc_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5gcreate - BoundaryConditions")
        end if


        !
        ! Loop through faces
        !
        faces = ["  XI_MIN","  XI_MAX"," ETA_MIN"," ETA_MAX","ZETA_MIN","ZETA_MAX"]
        do iface = 1,NFACES

            !
            ! Open boundary condition face group
            !
            call h5gopen_f(bc_id, trim(adjustl(faces(iface))), bcface_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5gopen_f - bcface (ex. XI_MIN)")


            !
            !  Get number of groups linked to the current bc_op
            !
            call h5gn_members_f(bcface_id, ".", nmembers, ierr)


            !
            !  Loop through members and add bc_operators if any exist
            !
            if (nmembers > 0) then

                do igroup = 0,nmembers-1

                    ! Get group name
                    call h5gget_obj_info_idx_f(bcface_id, ".", igroup, gname, type, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: error getting boundary condition group name")

                    ! Test if group is a boundary condition operator. 'BCS_'
                    if (gname(1:4) == 'BCS_') then
                        call bcnames(iface)%push_back(string_t(trim(gname(5:))))
                    end if

                end do

            end if


            ! Set empty string if none detected
            if (bcnames(iface)%size() == 0) then
                call bcnames(iface)%push_back(string_t('empty'))
            end if

            
            ! Close boundary condition face group
            call h5gclose_f(bcface_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5gclose.")


        end do !iface


        ! Close the boundary condition group
        call h5gclose_f(bc_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5gclose.")


    end function get_bcnames_hdf
    !****************************************************************************************










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







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/25/2016
    !!
    !!
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
        msg = "Attribute "//trim(attribute)//" not found in the file. Maybe the file was generated with an old &
               version of the ChiDG library. Try regenerating the HDF grid file with an updated version &
               of the ChiDG library to make sure the file is formatted properly"

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








end module mod_hdf_utilities
