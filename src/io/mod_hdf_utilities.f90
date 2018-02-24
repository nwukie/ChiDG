!----------------------------------------------------------------------------------------
!!
!!  ChiDG HDF File Format API
!!
!!  =====================================================================
!!  HDF:
!!  =====================================================================
!!  open_hdf
!!  close_hdf
!!
!!
!!  =====================================================================
!!  File:
!!      A: Storage Version Major    (Integer)
!!      A: Storage Version Minor    (Integer)
!!      A: Contains Grid            (Logical)
!!      A: Contains Solution        (Logical)
!!  
!!      G: BCSG_name1 ; ...  ('BCSG_' = Boundary condition state groups)
!!      G: EQN_name1  ; ...  ('EQN_'  = Equations)
!!      G: D_name1    ; ...  ('D_'    = Domains)
!!
!!  =====================================================================
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
!!  set_time_step_hdf
!!  get_time_step_hdf
!!
!!  set_nsteps_hdf
!!  get_nsteps_hdf
!!
!!  set_nwrite_hdf
!!  get_nwrite_hdf
!!
!!  set_frequencies_hdf
!!  get_frequencies_hdf
!!
!!  !set_time_levels_hdf
!!  !get_time_levels_hdf
!!
!!  
!!
!!  =====================================================================
!!  Domain-level routines:
!!  "D_"
!!  ---------------------------------------------------------------------
!!      A: Name                     (String)
!!      A: Equation Set             (String)
!!      A: Coordinate System        (String)
!!      A: Coordinate Order         (Integer)
!!      A: Field Order              (Integer)
!!      A: Dimensionality           (Integer)
!!      G: Grid
!!      G: Patches
!!      G: Fields
!!  =====================================================================
!!  add_domain_hdf
!!  create_domain_hdf
!!  open_domain_hdf
!!  close_domain_hdf
!!  check_domain_exists_hdf
!!
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
!!  set_domain_coordinate_displacements_hdf
!!  get_domain_coordinate_displacements_hdf
!!  set_domain_coordinate_velocities_hdf
!!  get_domain_coordinate_velocities_hdf
!!
!!  set_domain_coordinate_system_hdf
!!  get_domain_coordinate_system_hdf
!!
!!  set_domain_connectivity_hdf
!!  get_domain_connectivity_hdf
!!  set_domain_connectivity_partition_hdf
!!
!!  get_domain_nelements_hdf
!!  get_domain_nnodes_hdf
!!
!!  set_domain_coordinate_order_hdf
!!  get_domain_coordinate_order_hdf
!!  get_domain_coordinate_orders_hdf
!!
!!  set_domain_field_order_hdf
!!  get_domain_field_order_hdf
!!  get_domain_field_orders_hdf
!!
!!  TODO: set_domain_field_order_hdf
!!  TODO: get_domain_field_order_hdf
!!  TODO: get_domain_field_orders_hdf
!!
!!  set_domain_dimensionality_hdf   DEPRECATED
!!  get_domain_dimensionality_hdf   DEPRECATED
!!  get_domain_dimensionalities_hdf DEPRECATED
!!
!!
!!  set_domain_equation_set_hdf
!!  get_domain_equation_set_hdf
!!  get_domain_equation_sets_hdf
!!
!!
!!  create_patch_hdf
!!  open_patch_hdf
!!  close_patch_hdf
!!  set_patch_hdf
!!  get_patch_hdf
!!  set_patch_group_hdf
!!  get_patch_group_hdf
!!  
!!
!!  get_npatches_hdf
!!  copy_patches_attributes_hdf
!!
!!
!!
!!  =====================================================================
!!  Boundary Conditions: 
!!  "BCSG_", "BCS_", "BCP_"
!!  ---------------------------------------------------------------------
!!  =====================================================================
!!  "BCSG_"
!!  create_bc_state_group_hdf
!!  open_bc_state_group_hdf
!!  close_bc_state_group_hdf
!!
!!  get_nbc_state_groups_hdf
!!  get_bc_state_group_names_hdf
!!
!!  copy_bc_state_groups_hdf
!!  remove_bc_state_group_hdf
!!
!!
!!  "BCS_"
!!  open_bc_state_hdf
!!  close_bc_state_hdf
!!  add_bc_state_hdf
!!
!!  get_bc_states_hdf
!!  get_nbc_states_hdf
!!  get_bc_state_names_hdf
!!  remove_bc_state_hdf
!!  check_bc_state_exists_hdf
!!
!!
!!  "BCP_"
!!  add_bc_properties_hdf
!!  get_bc_properties_hdf
!!  remove_bc_property_hdf
!!  check_bc_property_exists_hdf
!!
!!
!!  get_bcnames_hdf
!!
!!
!!  =====================================================================
!!  Equations:
!!  ---------------------------------------------------------------------
!!  "EQN_"
!!  =====================================================================
!!  create_eqn_group_hdf
!!  remove_eqn_group_hdf
!!  get_eqn_group_names_hdf
!!  check_eqn_group_exists_hdf
!!  prune_eqn_groups_hdf
!!
!!
!!  Time Integrators:
!!  ---------------------------------------------------------------------
!!  set_time_integrator_hdf
!!  get_time_integrator_hdf
!!
!!  !set_ntimes_hdf
!!  !get_ntimes_hdf
!!
!!  set_times_hdf
!!  get_times_hdf
!!
!!  set_frequencies_hdf
!!  get_frequencies_hdf
!!
!!
!!  Utilities:
!!  ---------------------------------------------------------------------
!!  delete_group_attributes_hdf
!!  check_attribute_exists_hdf
!!  check_link_exists_hdf
!!  check_file_exists_hdf
!!  check_file_has_extension_hdf
!!  check_file_structure_hdf      
!!  check_file_domains_hdf
!!  check_file_ntimes_hdf
!!
!****************************************************************************************
module mod_hdf_utilities
#include <messenger.h>
    use mod_kinds,              only: rk, ik, rdouble
    use mod_constants,          only: ZERO, NFACES, TWO_DIM, THREE_DIM
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


    use type_prescribed_mesh_motion,                only: prescribed_mesh_motion_t
    use type_prescribed_mesh_motion_group,          only: prescribed_mesh_motion_group_t
    use type_prescribed_mesh_motion_group_wrapper,  only: prescribed_mesh_motion_group_wrapper_t
    implicit none

    

    !
    ! HDF5 storage format
    !
    integer, parameter :: STORAGE_VERSION_MAJOR = 1
    integer, parameter :: STORAGE_VERSION_MINOR = 6


    ! Attribute sizes
    integer(HSIZE_T), parameter :: SIZE_ONE = 1

    
    ! HDF library status
    logical    :: HDF_is_open     = .false.
    integer    :: HDF_nfiles_open = 0


contains





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


        filename_init = trim(filename)
        

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
        call set_storage_version_major_hdf(fid,STORAGE_VERSION_MAJOR)
        call set_storage_version_minor_hdf(fid,STORAGE_VERSION_MINOR)


        !
        ! Set contains status for grid/solution 
        !
        call set_contains_grid_hdf(fid,"False")
        call set_contains_solution_hdf(fid,"False")



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

        integer(ik)                 :: idom, eqn_ID
        integer(HID_T)              :: domain_id
        character(:),   allocatable :: domain_name

        do idom = 1,data%mesh%ndomains()

            ! Create domain group
            domain_name = data%mesh%domain(idom)%name
            call create_domain_hdf(fid,domain_name)


            ! Set additional attributes
            eqn_ID    = data%mesh%domain(idom)%elems(1)%eqn_ID !assume each element has same eqn_ID
            domain_id = open_domain_hdf(fid,trim(domain_name))
            call set_domain_equation_set_hdf(domain_id,data%eqnset(eqn_ID)%get_name())
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
    function open_file_hdf(filename,silent) result(fid)
        character(*),   intent(in)              :: filename
        logical,        intent(in), optional    :: silent

        character(:),   allocatable :: filename_open, user_msg
        integer(HID_T)  :: fid
        integer         :: ierr, loc
        logical         :: file_exists, file_is_hdf5

        ! Trim if we need to
        filename_open = trim(filename)

        ! Check file exists
        inquire(file=filename_open, exist=file_exists)
        if (.not. file_exists) then
            call chidg_signal_one(FATAL,"open_file_hdf: Could not find file.",trim(filename_open))
        end if

        ! Check file is HDF5/ChiDG
        call h5fis_hdf5_f(filename_open,file_is_hdf5,ierr)
        user_msg = "open_file_hdf: Specified file is not an HDF5 file."
        if (.not. file_is_hdf5) call chidg_signal_one(FATAL,user_msg,trim(filename_open))

        !
        !  Open input file using default properties.
        !
        call open_hdf()
        call write_line('   Opening file: ', filename_open, io_proc=GLOBAL_MASTER, ltrim=.false., silence=silent)
        call h5fopen_f(filename_open, H5F_ACC_RDWR_F, fid, ierr)
        user_msg = "open_file_hdf: There was an error opening the file."
        if (ierr /= 0) call chidg_signal_one(FATAL,user_msg,trim(filename_open))


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
    subroutine close_file_hdf(fid,silent)
        integer(HID_T), intent(in)              :: fid
        logical,        intent(in), optional    :: silent

        character(1024) :: buf
        integer(SIZE_T) :: name_size
        integer         :: ierr


        ! Get file name
        call h5fget_name_f(fid, buf, name_size, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"close_file_hdf: error getting file name from identifier.")
        call write_line('   Closing file: ', trim(buf), io_proc=GLOBAL_MASTER, ltrim=.false., silence=silent)


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
        if ( (file_major_version /= STORAGE_VERSION_MAJOR) .or. &
             (file_minor_version /= STORAGE_VERSION_MINOR) ) then

            msg = "The storage format of the file being worked with &
                   does not match the storage format in the ChiDG library being used. This &
                   probably means the file was generated with another version of the ChiDG &
                   library and may not be compatible with the library currently being used. &
                   You could try a few things here. 1: regenerate the with with the ChiDG &
                   library being used. 2: Use a different version of the ChiDG library &
                   that uses a storage format for the file being used. 3: Full-speed ahead! &
                   Proceed anyways and try your luck!"//NEW_LINE('A')//" &
                   Options: Exit(1), Continue(2)."
                   

            call chidg_signal_two(MSG,msg,file_minor_version,STORAGE_VERSION_MINOR)

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











!    !>  Check if an already existing HDF5 file has the correct structure to be use to write 
!    !!  a new solution
!    !!
!    !!  Errors if incompatible.
!    !!
!    !!  @author Matteo Ugolotti
!    !!  @date   02/20/2017
!    !!
!    !!  TODO: add check to the domains' attribute (dimensionality, equation-set)
!    !!        Many other check can be implemented. This is the first step
!    !-----------------------------------------------------------------------------------------
!    subroutine check_file_structure_hdf(fid,data)
!    
!        integer(HID_T)      , intent(in)  :: fid
!        type(chidg_data_t)  , intent(in)  :: data
!        
!        ! Check if the version is correct
!        call check_file_storage_version_hdf(fid)
!    
!        ! Check if the number and names of the domains are correct
!        call check_file_domains_hdf(fid,data)
!        
!        ! Check the number of time levels, and update it in case.
!        call check_file_ntimes_hdf(fid,data)
!
!
!    end subroutine check_file_structure_hdf
!    !*****************************************************************************************









!    !>  Check check that the ntimes attirbute of the existing file matches the  
!    !!  a new data solution
!    !!
!    !!  Errors if incompatible.
!    !!
!    !!  @author Matteo Ugolotti
!    !!  @date   02/20/2017
!    !!
!    !-----------------------------------------------------------------------------------------
!    subroutine check_file_ntimes_hdf(fid,data)
!        
!        integer(HID_T)      , intent(in)    :: fid
!        type(chidg_data_t)  , intent(in)    :: data
!
!        integer(ik)                         :: time_lev, ntimes
!        character(len=1024),    allocatable :: msg
!
!        time_lev = get_ntimes_hdf(fid)
!        ntimes   = data%ntime() - 1
!    
!
!        if ( time_lev /= ntimes ) then
!            
!            msg = "The attribute ntimes in the existing  HDF5 fils is not up-to-date. It will automatically be updated"
!            call chidg_signal(WARN, msg)
!
!            !update ntimes
!            call set_ntimes_hdf(fid,ntimes)
!
!        end if
!
!
!    end subroutine check_file_ntimes_hdf
!    !*****************************************************************************************










!    !>  Check if the number of domains and name in an old HDF5 file are coherent with the
!    !!  new solution
!    !!
!    !!  Errors if incompatible.
!    !!
!    !!  @author Matteo Ugolotti
!    !!  @date   02/20/2017
!    !!
!    !!
!    !!
!    !-----------------------------------------------------------------------------------------
!    subroutine check_file_domains_hdf(id,data)
!
!        integer(HID_T)      ,intent(in) :: id
!        type(chidg_data_t)  ,intent(in) :: data
!
!        integer(ik)                        :: domain_number, idom, domain_check, dimension
!        character(len=1024),   allocatable :: dom_names(:), domain_name, msg
!        
!        !
!        ! check whether the number of domians is correct or not. IF it is, then 
!        ! compare domains' names
!        !
!        
!        domain_number = data%ndomains()
!        domain_check  = get_ndomains_hdf(id)
!
!        if (domain_number /= domain_check) then
!            msg = "The number of domains in the existing file and the number of domains that need to be stored &
!                    do not match. Please, delete the existing HDF5 file."
!            call chidg_signal(FATAL,msg)
!        else
!
!            !
!            ! check if the name of the domainis are correct
!            !
!            dom_names = get_domain_names_hdf(id)
!
!            do idom = 1, domain_number
!                
!                domain_name = "D_"//trim(data%info(idom)%name)
!                        
!                msg = "There exists a mismatch in a domain's name due to the overwriting of and existing HDF5 &
!                       file. To avoid this delete the existing HDF5 file."
!
!                if ( domain_name /= dom_names(idom) ) call chidg_signal(FATAL,msg)
!
!            end do
!        
!        end if
!        
!
!
!    end subroutine check_file_domains_hdf
!    !*****************************************************************************************





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

        call h5ltset_attribute_int_f(fid, "/", 'Storage Version Major', [major_version], SIZE_ONE, ierr)
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

        call h5ltset_attribute_int_f(fid, "/", 'Storage Version Minor', [minor_version], SIZE_ONE, ierr)
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


        call check_attribute_exists_hdf(fid,"Storage Version Major","Soft Fail",attr_status)

        if (attr_status == 0) then
            call h5ltget_attribute_int_f(fid, "/", 'Storage Version Major', buf, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_storage_version_major_hdf: error getting the major storage version.")
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


        call check_attribute_exists_hdf(fid,"Storage Version Minor","Soft Fail",attr_status)

        if (attr_status == 0) then
            call h5ltget_attribute_int_f(fid, "/", 'Storage Version Minor', buf, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_storage_version_minor_hdf: error getting the minor storage version.")
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
        integer(HSIZE_T)            :: ntime, nfreq

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
        call open_hdf()



        !
        !  Open input file using default properties.
        !
        fid = open_file_hdf(filename)



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
        ! Get time_integrator name
        !
        prop%time_integrator = get_time_integrator_hdf(fid)


        !
        ! Get order of coordinate and solution expansions
        !
        if ( prop%contains_grid ) then
            prop%order_c = get_domain_coordinate_orders_hdf(fid,prop%domain_names)
        end if

        if ( prop%contains_solution ) then
            prop%order_s = get_domain_field_orders_hdf(fid)
        end if




        !
        ! Compute number of terms in the polynomial expansions for each domain
        !
        do idom = 1,prop%ndomains
            

            nterms_1d = (prop%order_c(idom) + 1)
            prop%nterms_c(idom) = nterms_1d * nterms_1d * nterms_1d
            !if ( prop%spacedim(idom) == THREE_DIM ) then
            !    prop%nterms_c(idom) = nterms_1d * nterms_1d * nterms_1d
            !else if ( prop%spacedim(idom) == TWO_DIM ) then
            !    prop%nterms_c(idom) = nterms_1d * nterms_1d
            !end if


 
            nterms_1d = (prop%order_s(idom) + 1)
            prop%nterms_s(idom) = nterms_1d * nterms_1d * nterms_1d
            !if ( prop%spacedim(idom) == THREE_DIM ) then
            !    prop%nterms_s(idom) = nterms_1d * nterms_1d * nterms_1d
            !else if ( prop%spacedim(idom) == TWO_DIM ) then
            !    prop%nterms_s(idom) = nterms_1d * nterms_1d
            !end if


        end do ! idom


        !
        ! Get equation set for each domain
        !
        prop%eqnset = get_domain_equation_sets_hdf(fid)

        

        !
        ! Close file
        !
        call close_file_hdf(fid)
        call close_hdf()

    end function get_properties_hdf
    !****************************************************************************************










    !>  Add a domain to the file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine add_domain_hdf(fid,domain_name,nodes,elements,coord_system,equation_set)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: domain_name
        real(rk),       intent(in)  :: nodes(:,:)
        integer(ik),    intent(in)  :: elements(:,:)
        character(*),   intent(in)  :: coord_system
        character(*),   intent(in)  :: equation_set


        real(rk),   allocatable :: dnodes(:,:), vnodes(:,:)
        integer(HID_T)          :: dom_id, grid_id, bc_id, var_id
        integer(ik)             :: mapping, ierr


        !
        ! Create new domain. /D_domain_name
        !
        call create_domain_hdf(fid,domain_name)


        !
        ! Size perturbation array and set to default=zero.
        !  
        dnodes = nodes
        vnodes = nodes
        dnodes = ZERO
        vnodes = ZERO


        !
        ! Open Domain
        !
        dom_id = open_domain_hdf(fid,domain_name)



        ! Set 'Coordinate System', 'Coordinates', 'DCoordinates', 'VCoordinates', and 'Elements'
        mapping  = elements(1,3)
        call set_domain_coordinate_order_hdf(dom_id,mapping)
        call set_domain_coordinate_system_hdf(dom_id,coord_system)
        call set_domain_coordinates_hdf(dom_id,nodes)
        call set_domain_coordinate_displacements_hdf(dom_id,dnodes)
        call set_domain_coordinate_velocities_hdf(dom_id,vnodes)
        call set_domain_connectivity_hdf(dom_id,elements)
        call set_domain_equation_set_hdf(dom_id,trim(equation_set))



        !
        ! Close Domain
        !
        call close_domain_hdf(dom_id)


        !
        ! Write equation set group to file root. /EQN_equation_set
        !
        call create_eqn_group_hdf(fid, trim(equation_set))


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

        integer(ik)     :: ierr
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

            call h5gcreate_f(domain_id, "Patches", bc_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"create_domain_hdf: h5gcreate_f")
            call h5gclose_f(bc_id,ierr)

            call h5gcreate_f(domain_id, "Fields", var_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"create_domain_hdf: h5gcreate_f")
            call h5gclose_f(var_id,ierr)


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





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/15/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    function check_domain_exists_hdf(fid,domain_name) result(exist_status)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: domain_name

        integer(ik) :: ierr
        logical     :: exist_status

        ! Check if face contains the bc_state
        call h5lexists_f(fid, "D_"//trim(domain_name), exist_status, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"check_domain_exists_hdf: Error in call to h5lexists_f")


    end function check_domain_exists_hdf
    !*********************************************************************************************






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

        integer                     :: nmembers, ierr
        integer(ik)                 :: ndomains
        integer(HSIZE_T)            :: igrp
        character(1024)             :: gname
        character(:),   allocatable :: user_msg


        !  Get number of groups in the file root
        call h5gn_members_f(fid, "/", nmembers, ierr)


        !  Loop through groups and count those starting with "D_"
        ndomains = 0
        do igrp = 0,nmembers-1

            ! Get group name
            call h5lget_name_by_idx_f(fid,".",H5_INDEX_NAME_F,H5_ITER_INC_F,igrp,gname,ierr)
            user_msg = "get_ndomains_hdf: Error iterating through links to detect domain groups."
            if (ierr /= 0) call chidg_signal(FATAL,user_msg)

            ! Test if group is a 'Domain'
            if (gname(1:2) == 'D_') then
                ndomains = ndomains + 1
            end if

        end do


    end function get_ndomains_hdf
    !****************************************************************************************









!    !>  Given a file identifier, return the time integrator name in a hdf5 file.
!    !!
!    !!  @author Mayank Sharma
!    !!  @date   3/18/2017
!    !!
!    !!  @param[in]  fid     HDF file identifier
!    !!
!    !----------------------------------------------------------------------------------------
!    function get_time_integrator_hdf(fid) result(time_string)
!        integer(HID_T), intent(in)  :: fid
!        
!        character(1024)             :: temp_string
!        character(:),   allocatable :: time_string
!        integer(ik)                 :: ierr
!
!        call h5ltget_attribute_string_f(fid, "/", "time_integrator", temp_string, ierr)
!        if (ierr /= 0) call chidg_signal(FATAL,"get_time_integrator_hdf: h5ltget_attribute_string_f & 
!                                         had a problem getting the time integrator name")
!        time_string = trim(temp_string)
!
!    end function get_time_integrator_hdf
!    !****************************************************************************************









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

        call h5ltset_attribute_string_f(domain_id, ".", "Name", trim(domain_name), ierr)
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

        call h5ltget_attribute_string_f(domain_id, ".", "Name", temp_string, ierr)
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

        character(:),           allocatable :: user_msg
        character(len=1024),    allocatable :: names(:)
        character(len=1024)                 :: gname
        integer(HSIZE_T)                    :: igrp
        integer                             :: ndomains, nmembers, type
        integer                             :: idom, ierr

        ! Get number of domains
        ndomains = get_ndomains_hdf(fid)
        allocate(names(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError


        !  Get number of groups in the file root
        call h5gn_members_f(fid, "/", nmembers, ierr)


        !  Loop through groups and read domain names
        idom = 1
        do igrp = 0,nmembers-1

            ! Get group name
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
        real(rk),       intent(in)  :: nodes(:,:)

        integer(HID_T)      :: grid_id, xspace_id, yspace_id, zspace_id, xset_id, yset_id, zset_id
        integer(HSIZE_T)    :: dims_rank_one(1)
        integer(ik)         :: ipt, ierr
        logical             :: exists

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
        dims_rank_one = size(nodes,1)


        exists = check_link_exists_hdf(grid_id,"Coordinate1")
        if (exists) then
            !
            ! Open 'Coordinate' data sets
            !
            call h5dopen_f(grid_id, "Coordinate1", xset_id, ierr, H5P_DEFAULT_F)
            call h5dopen_f(grid_id, "Coordinate2", yset_id, ierr, H5P_DEFAULT_F)
            call h5dopen_f(grid_id, "Coordinate3", zset_id, ierr, H5P_DEFAULT_F)

        else

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
            ! Close dataspaces
            !
            call h5sclose_f(xspace_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5sclose_f")
            call h5sclose_f(yspace_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5sclose_f")
            call h5sclose_f(zspace_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5sclose_f")



        end if

        !
        ! Write coordinates to datasets
        !
        call h5dwrite_f(xset_id, H5T_NATIVE_DOUBLE, nodes(:,1), dims_rank_one, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5dwrite_f")
        call h5dwrite_f(yset_id, H5T_NATIVE_DOUBLE, nodes(:,2), dims_rank_one, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5dwrite_f")
        call h5dwrite_f(zset_id, H5T_NATIVE_DOUBLE, nodes(:,3), dims_rank_one, ierr)
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
        ! Close Grid group
        !
        call h5gclose_f(grid_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinates_her: h5gclose_f")


    end subroutine set_domain_coordinates_hdf
    !***************************************************************************************










    !>  For a domain, set the displacements for the domain coordinates.
    !!
    !!  /D_domainname/Grid/DCoordinate1
    !!  /D_domainname/Grid/DCoordinate2
    !!  /D_domainname/Grid/DCoordinate3
    !!
    !!  This is used in the context of the ALE formulation, which is formulated using a map
    !!  from original to deformed elements. 
    !!
    !!  The datasets DCoordinate1, DCoordinate2, DCoordinate3, represent the displacements
    !!  of the original nodes to the deformed nodes. By default, the displacement values
    !!  are present in the file, but are zero.
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   6/15/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_coordinate_displacements_hdf(dom_id,nodes)
        integer(HID_T), intent(in)  :: dom_id
        real(rk),       intent(in)  :: nodes(:,:)

        integer(HID_T)      :: grid_id, xspace_id, yspace_id, zspace_id, xset_id, yset_id, zset_id
        integer(HSIZE_T)    :: dims_rank_one(1)
        integer(ik)         :: ipt, ierr
        logical             :: exists

        !
        ! Create a grid-group within the current block domain
        !
        exists = check_link_exists_hdf(dom_id,"Grid")
        if (exists) then
            call h5gopen_f(dom_id,"Grid", grid_id,ierr)
        else
            call h5gcreate_f(dom_id, "Grid", grid_id, ierr)
        end if
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_displacements_hdf: h5gcreate_f")




        !
        ! Re-order coordinates to be linear arrays
        !
        dims_rank_one = size(nodes,1)


        exists = check_link_exists_hdf(grid_id,"DCoordinate1")
        if (exists) then
            !
            ! Open 'Coordinate' data sets
            !
            call h5dopen_f(grid_id, "DCoordinate1", xset_id, ierr, H5P_DEFAULT_F)
            call h5dopen_f(grid_id, "DCoordinate2", yset_id, ierr, H5P_DEFAULT_F)
            call h5dopen_f(grid_id, "DCoordinate3", zset_id, ierr, H5P_DEFAULT_F)

        else

            !
            ! Create dataspaces for grid coordinates
            !
            call h5screate_simple_f(1, dims_rank_one, xspace_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_displacements_hdf: h5screate_simple_f")
            call h5screate_simple_f(1, dims_rank_one, yspace_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_displacements_hdf: h5screate_simple_f")
            call h5screate_simple_f(1, dims_rank_one, zspace_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_displacements_hdf: h5screate_simple_f")


            !
            ! Create datasets for grid coordinates
            !
            call h5dcreate_f(grid_id, "DCoordinate1", H5T_NATIVE_DOUBLE, xspace_id, xset_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_displacements_hdf: h5dcreate_f")
            call h5dcreate_f(grid_id, "DCoordinate2", H5T_NATIVE_DOUBLE, yspace_id, yset_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_displacements_hdf: h5dcreate_f")
            call h5dcreate_f(grid_id, "DCoordinate3", H5T_NATIVE_DOUBLE, zspace_id, zset_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_displacements_hdf: h5dcreate_f")


            !
            ! Close dataspaces
            !
            call h5sclose_f(xspace_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_displacements_hdf: h5sclose_f")
            call h5sclose_f(yspace_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_displacements_hdf: h5sclose_f")
            call h5sclose_f(zspace_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_displacements_hdf: h5sclose_f")



        end if

        !
        ! Write coordinates to datasets
        !
        call h5dwrite_f(xset_id, H5T_NATIVE_DOUBLE, nodes(:,1), dims_rank_one, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_displacements_hdf: h5dwrite_f")
        call h5dwrite_f(yset_id, H5T_NATIVE_DOUBLE, nodes(:,2), dims_rank_one, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_displacements_hdf: h5dwrite_f")
        call h5dwrite_f(zset_id, H5T_NATIVE_DOUBLE, nodes(:,3), dims_rank_one, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_displacements_hdf: h5dwrite_f")


        !
        ! Close datasets
        !
        call h5dclose_f(xset_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinate_displacements_hdf: h5dclose_f")
        call h5dclose_f(yset_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinate_displacements_hdf: h5dclose_f")
        call h5dclose_f(zset_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinate_displacements_hdf: h5dclose_f")


        !
        ! Close Grid group
        !
        call h5gclose_f(grid_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinate_displacements_hdf: h5gclose_f")


    end subroutine set_domain_coordinate_displacements_hdf
    !***************************************************************************************








    !>  For a domain, set the velocities for the domain coordinates.
    !!
    !!  /D_domainname/Grid/VCoordinate1
    !!  /D_domainname/Grid/VCoordinate2
    !!  /D_domainname/Grid/VCoordinate3
    !!
    !!  This is used in the context of the ALE formulation, which is formulated using a map
    !!  from original to deformed elements. 
    !!
    !!  The datasets VCoordinate1, VCoordinate2, VCoordinate3, represent the velocities
    !!  of the deformed nodes. By default, the velocity values are present in the file, 
    !!  but are zero.
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   6/15/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_coordinate_velocities_hdf(dom_id,nodes)
        integer(HID_T), intent(in)  :: dom_id
        real(rk),       intent(in)  :: nodes(:,:)

        integer(HID_T)      :: grid_id, xspace_id, yspace_id, zspace_id, xset_id, yset_id, zset_id
        integer(HSIZE_T)    :: dims_rank_one(1)
        integer(ik)         :: ipt, ierr
        logical             :: exists

        !
        ! Create a grid-group within the current block domain
        !
        exists = check_link_exists_hdf(dom_id,"Grid")
        if (exists) then
            call h5gopen_f(dom_id,"Grid", grid_id,ierr)
        else
            call h5gcreate_f(dom_id, "Grid", grid_id, ierr)
        end if
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_velocities_hdf: h5gcreate_f")




        !
        ! Re-order coordinates to be linear arrays
        !
        dims_rank_one = size(nodes,1)


        exists = check_link_exists_hdf(grid_id,"VCoordinate1")
        if (exists) then
            !
            ! Open 'Coordinate' data sets
            !
            call h5dopen_f(grid_id, "VCoordinate1", xset_id, ierr, H5P_DEFAULT_F)
            call h5dopen_f(grid_id, "VCoordinate2", yset_id, ierr, H5P_DEFAULT_F)
            call h5dopen_f(grid_id, "VCoordinate3", zset_id, ierr, H5P_DEFAULT_F)

        else

            !
            ! Create dataspaces for grid coordinates
            !
            call h5screate_simple_f(1, dims_rank_one, xspace_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_velocities_hdf: h5screate_simple_f")
            call h5screate_simple_f(1, dims_rank_one, yspace_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_velocities_hdf: h5screate_simple_f")
            call h5screate_simple_f(1, dims_rank_one, zspace_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_velocities_hdf: h5screate_simple_f")


            !
            ! Create datasets for grid coordinates
            !
            call h5dcreate_f(grid_id, "VCoordinate1", H5T_NATIVE_DOUBLE, xspace_id, xset_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_velocities_hdf: h5dcreate_f")
            call h5dcreate_f(grid_id, "VCoordinate2", H5T_NATIVE_DOUBLE, yspace_id, yset_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_velocities_hdf: h5dcreate_f")
            call h5dcreate_f(grid_id, "VCoordinate3", H5T_NATIVE_DOUBLE, zspace_id, zset_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_velocities_hdf: h5dcreate_f")


            !
            ! Close dataspaces
            !
            call h5sclose_f(xspace_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_velocities_hdf: h5sclose_f")
            call h5sclose_f(yspace_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_velocities_hdf: h5sclose_f")
            call h5sclose_f(zspace_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_velocities_hdf: h5sclose_f")



        end if

        !
        ! Write coordinates to datasets
        !
        call h5dwrite_f(xset_id, H5T_NATIVE_DOUBLE, nodes(:,1), dims_rank_one, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_velocities_hdf: h5dwrite_f")
        call h5dwrite_f(yset_id, H5T_NATIVE_DOUBLE, nodes(:,2), dims_rank_one, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_velocities_hdf: h5dwrite_f")
        call h5dwrite_f(zset_id, H5T_NATIVE_DOUBLE, nodes(:,3), dims_rank_one, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dcoordinate_velocities_hdf: h5dwrite_f")


        !
        ! Close datasets
        !
        call h5dclose_f(xset_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinate_velocities_hdf: h5dclose_f")
        call h5dclose_f(yset_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinate_velocities_hdf: h5dclose_f")
        call h5dclose_f(zset_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinate_velocities_hdf: h5dclose_f")


        !
        ! Close Grid group
        !
        call h5gclose_f(grid_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinate_velocities_hdf: h5gclose_f")


    end subroutine set_domain_coordinate_velocities_hdf
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

        real(rk),      allocatable  :: nodes(:,:)

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
        allocate(nodes(npts,3), stat=ierr)
        if (ierr /= 0) call AllocationError
            

        do ipt = 1,rank_one_dims(1)
            nodes(ipt,1) = real(pts1(ipt),rk)
            nodes(ipt,2) = real(pts2(ipt),rk)
            nodes(ipt,3) = real(pts3(ipt),rk)
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






    !>  For a domain, get the displacements for the domain coordinates.
    !!
    !!  /D_domainname/Grid/DCoordinate1
    !!  /D_domainname/Grid/DCoordinate2
    !!  /D_domainname/Grid/DCoordinate3
    !!
    !!  This is used in the context of the ALE formulation, which is formulated using a map
    !!  from original to deformed elements. 
    !!
    !!  The datasets DCoordinate1, DCoordinate2, DCoordinate3, represent the displacements
    !!  of the original nodes to the deformed nodes. By default, the displacement values
    !!  are present in the file, but are zero.
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   6/15/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_coordinate_displacements_hdf(dom_id) result(nodes)
        integer(HID_T), intent(in)  :: dom_id

        real(rk),      allocatable  :: nodes(:,:)

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
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_displacements_hdf: h5gopen_f")
        else
            call chidg_signal(FATAL,"get_domain_coordinate_displacements_hdf: The current domain does not contain a 'Grid'.")
        end if


        !
        !  Open the Coordinate datasets
        !
        call h5dopen_f(grid_id, "DCoordinate1", did_1, ierr, H5P_DEFAULT_F)
        call h5dopen_f(grid_id, "DCoordinate2", did_2, ierr, H5P_DEFAULT_F)
        call h5dopen_f(grid_id, "DCoordinate3", did_3, ierr, H5P_DEFAULT_F)


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
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_displacements_hdf: h5dread_f")

        cp_pts = c_loc(pts2(1))
        call h5dread_f(did_2, H5T_NATIVE_DOUBLE, cp_pts, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_displacements_hdf: h5dread_f")

        cp_pts = c_loc(pts3(1))
        call h5dread_f(did_3, H5T_NATIVE_DOUBLE, cp_pts, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_displacements_hdf: h5dread_f")



        !
        !  Accumulate pts into a single points_t matrix to initialize domain
        !
        allocate(nodes(npts,3), stat=ierr)
        if (ierr /= 0) call AllocationError
            

        do ipt = 1,rank_one_dims(1)
            nodes(ipt,1) = real(pts1(ipt),rk)
            nodes(ipt,2) = real(pts2(ipt),rk)
            nodes(ipt,3) = real(pts3(ipt),rk)
        end do


        ! Close the Coordinate datasets
        call h5dclose_f(did_1,ierr)
        call h5dclose_f(did_2,ierr)
        call h5dclose_f(did_3,ierr)

        ! Close the dataspace id
        call h5sclose_f(sid,ierr)

        ! Close Grid group
        call h5gclose_f(grid_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_displacements_hdf: h5gclose_f")


    end function get_domain_coordinate_displacements_hdf
    !***************************************************************************************





    !>  For a domain, get the velocities for the domain coordinates.
    !!
    !!  /D_domainname/Grid/VCoordinate1
    !!  /D_domainname/Grid/VCoordinate2
    !!  /D_domainname/Grid/VCoordinate3
    !!
    !!  This is used in the context of the ALE formulation, which is formulated using a map
    !!  from original to deformed elements. 
    !!
    !!  The datasets VCoordinate1, VCoordinate2, VCoordinate3, represent the velocities
    !!  of the the deformed nodes. By default, the velocity values are present in the file, 
    !!  but are zero.
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   6/15/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_coordinate_velocities_hdf(dom_id) result(nodes)
        integer(HID_T), intent(in)  :: dom_id

        real(rk),      allocatable  :: nodes(:,:)

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
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_velocities_hdf: h5gopen_f")
        else
            call chidg_signal(FATAL,"get_domain_coordinate_velocities_hdf: The current domain does not contain a 'Grid'.")
        end if


        !
        !  Open the Coordinate datasets
        !
        call h5dopen_f(grid_id, "VCoordinate1", did_1, ierr, H5P_DEFAULT_F)
        call h5dopen_f(grid_id, "VCoordinate2", did_2, ierr, H5P_DEFAULT_F)
        call h5dopen_f(grid_id, "VCoordinate3", did_3, ierr, H5P_DEFAULT_F)


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
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_velocities_hdf: h5dread_f")

        cp_pts = c_loc(pts2(1))
        call h5dread_f(did_2, H5T_NATIVE_DOUBLE, cp_pts, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_velocities_hdf: h5dread_f")

        cp_pts = c_loc(pts3(1))
        call h5dread_f(did_3, H5T_NATIVE_DOUBLE, cp_pts, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_velocities_hdf: h5dread_f")



        !
        !  Accumulate pts into a single points_t matrix to initialize domain
        !
        allocate(nodes(npts,3), stat=ierr)
        if (ierr /= 0) call AllocationError
            

        do ipt = 1,rank_one_dims(1)
            nodes(ipt,1) = real(pts1(ipt),rk)
            nodes(ipt,2) = real(pts2(ipt),rk)
            nodes(ipt,3) = real(pts3(ipt),rk)
        end do


        ! Close the Coordinate datasets
        call h5dclose_f(did_1,ierr)
        call h5dclose_f(did_2,ierr)
        call h5dclose_f(did_3,ierr)

        ! Close the dataspace id
        call h5sclose_f(sid,ierr)

        ! Close Grid group
        call h5gclose_f(grid_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_velocities_hdf: h5gclose_f")


    end function get_domain_coordinate_velocities_hdf
    !***************************************************************************************













    !>  For a domain, set the domain coordinate system.
    !!
    !!  Attribute: /D_domainname/Coordinate System
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
        integer(ik)                 :: ierr
        logical                     :: cartesian, cylindrical



        !
        ! Check valid input
        !
        cartesian   = (coord_system == 'Cartesian')
        cylindrical = (coord_system == 'Cylindrical')
        user_msg = "set_domain_coordinate_system_hdf: Invalid coordinate system. Valid systems are 'Cartesian' or 'Cylindrical'."
        if (.not. (cartesian .or. cylindrical)) call chidg_signal_one(FATAL,user_msg,coord_system)


        !
        ! Set 'Coordinate System' attribute
        !
        call h5ltset_attribute_string_f(dom_id,".","Coordinate System",coord_system,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinate_system_hdf: Error setting 'Coordinate System' attribute.")



    end subroutine set_domain_coordinate_system_hdf
    !***************************************************************************************




    !>  For a domain, get the domain coordinate system.
    !!
    !!  Attribute: /D_domainname/Coordinate System
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
        integer(ik)                 :: ierr

        call h5ltget_attribute_string_f(dom_id,".","Coordinate System",coord_system,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_system_hdf: Error getting 'Coordinate System' attribute.")
        
        ! Trim blanks
        coord_system_trim = trim(coord_system)

    end function get_domain_coordinate_system_hdf
    !***************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/25/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_coordinate_order_hdf(dom_id,order)
        integer(HID_T), intent(in)  :: dom_id
        integer(ik),    intent(in)  :: order

        integer(ik) :: ierr

        call h5ltset_attribute_int_f(dom_id,".","Coordinate Order", [order],SIZE_ONE,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_coordinate_order_hdf: Error setting 'Coordinate Order' attribute")

    end subroutine set_domain_coordinate_order_hdf
    !****************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/25/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_domain_coordinate_order_hdf(dom_id) result(order)
        integer(HID_T), intent(in)  :: dom_id

        integer                 :: ierr
        integer, dimension(1)   :: buf
        integer(ik)             :: order

        call h5ltget_attribute_int_f(dom_id,".","Coordinate Order", buf, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_order_hdf: Error getting 'Coordinate Order' attribute")

        order = int(buf(1),kind=ik)

    end function get_domain_coordinate_order_hdf
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
    function get_domain_coordinate_orders_hdf(fid, dnames) result(orders)
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

            !  Get coordinate mapping
            call h5ltget_attribute_int_f(fid, "D_"//trim(dnames(idom)), "Coordinate Order", buf, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_coordinate_orders_hdf: h5ltget_attribute_int_f")

            ! Compute number of terms in coordinate expansion
            mapping = buf(1)
            orders(idom) = int(mapping, kind=ik)

        end do

    end function get_domain_coordinate_orders_hdf
    !****************************************************************************************







    !>  Set the 'Field Order' attribute on a domain.
    !!
    !!  'Field Order' indicates the order of the polynomial expansion representing
    !!  the solution on the domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/25/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_field_order_hdf(dom_id,order)
        integer(HID_T), intent(in)  :: dom_id
        integer(ik),    intent(in)  :: order

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(dom_id,".","Field Order", [order],SIZE_ONE,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_field_order_hdf: Error setting 'Field Order' attribute")

    end subroutine set_domain_field_order_hdf
    !****************************************************************************************





    !>  Get the 'Field Order' attribute on a domain.
    !!
    !!  'Field Order' indicates the order of the polynomial expansion representing
    !!  the solution on the domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/25/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    function get_domain_field_order_hdf(dom_id) result(order)
        integer(HID_T), intent(in)  :: dom_id

        integer(ik) :: order, ierr
        integer, dimension(1) :: buf

        call h5ltget_attribute_int_f(dom_id,".","Field Order",buf,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_field_order_hdf: Error getting 'Field Order' attribute")

        order = int(buf(1), kind=ik)

    end function get_domain_field_order_hdf
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
    function get_domain_field_orders_hdf(fid) result(orders)
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

            orders(idom) = get_domain_field_order_hdf(did)

            call close_domain_hdf(did)

        end do

    end function get_domain_field_orders_hdf
    !****************************************************************************************






!    !>  
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   4/11/2016
!    !!
!    !!  DEPRECATED: Nathan A. Wukie
!    !!
!    !---------------------------------------------------------------------------------------
!    subroutine set_domain_dimensionality_hdf(dom_id,dimensionality)
!        integer(HID_T), intent(in)  :: dom_id
!        integer(ik),    intent(in)  :: dimensionality
!
!        integer(ik)         :: ierr
!
!        !  Get coordinate mapping
!        call h5ltset_attribute_int_f(dom_id,".","Dimensionality",[dimensionality],SIZE_ONE,ierr)
!        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_dimensionality_hdf: h5ltset_attribute_int_f")
!
!    end subroutine set_domain_dimensionality_hdf
!    !****************************************************************************************
!
!
!
!
!    !>  
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   4/11/2016
!    !!
!    !!  DEPRECATED: Nathan A. Wukie
!    !!
!    !----------------------------------------------------------------------------------------
!    function get_domain_dimensionality_hdf(dom_id) result(dimensionality)
!        integer(HID_T), intent(in)  :: dom_id
!
!        integer(ik) :: dimensionality, ierr
!        integer, dimension(1)   :: buf
!
!        !  Get coordinate mapping
!        call h5ltget_attribute_int_f(dom_id,".","Dimensionality",buf,ierr)
!        if (ierr /= 0) call chidg_signal(FATAL,"get_domain_dimensionality_hdf: h5ltget_attribute_int_f")
!
!        dimensionality = int(buf(1),kind=ik)
!
!    end function get_domain_dimensionality_hdf
!    !****************************************************************************************
!
!
!
!
!
!
!
!
!
!
!    !>  Returns an array of integers that specifies the number of spatial dimensions to use 
!    !!  for every domain.
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   4/11/2016
!    !!
!    !!  DEPRECATED: Nathan A. Wukie
!    !!
!    !!  @param[in]  fid         HDF file identifier.
!    !!  @param[in]  dnames(:)   List of domain names to be interrogated. 
!    !!
!    !----------------------------------------------------------------------------------------
!    function get_domain_dimensionalities_hdf(fid, dnames) result(dimensionalities)
!        integer(HID_T),         intent(in)  :: fid
!        character(len=1024),    intent(in)  :: dnames(:)
!
!        integer(ik), allocatable    :: dimensionalities(:)
!        integer(ik)                 :: ierr, idom
!        integer, dimension(1)       :: dimensionality
!
!
!        !
!        ! Allocate storage for orders
!        !
!        allocate(dimensionalities(size(dnames)), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!        !
!        !  Loop through groups and read domains
!        !
!        do idom = 1,size(dnames)
!
!            call h5ltget_attribute_int_f(fid, "D_"//trim(dnames(idom)), "Dimensionality", dimensionality, ierr)
!            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_dimensionalities_hdf: Error h5ltget_attribute_int_f")
!
!            dimensionalities(idom) = dimensionality(1)
!
!        end do
!
!    end function get_domain_dimensionalities_hdf
!    !****************************************************************************************







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
    function get_domain_equation_sets_hdf(fid) result(eqnsets)
        integer(HID_T),         intent(in)  :: fid

        character(len=1024), allocatable    :: dnames(:)
        integer(HID_T)                      :: did
        logical                             :: eqnset_exists
        character(len=1024), allocatable    :: eqnsets(:)
        integer                             :: ierr, idom


        !
        ! Get domain names
        !
        dnames = get_domain_names_hdf(fid)


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









    !>  Set the element connectivities for a block, from a partition.
    !!
    !!  Accepts array of element connectivities as
    !!      nelements_g                         ! indicates the total number of elements 
    !!                                          ! in the unpartitioned domain
    !!  
    !!      elements(nelem, size_connectivity)  ! element connectivities on the local partition
    !!
    !!  /D_domainname/Grid/Elements
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/15/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_domain_connectivity_partition_hdf(dom_id,nelements_g,elements)
        integer(HID_T), intent(in)  :: dom_id
        integer(ik),    intent(in)  :: nelements_g
        integer(ik),    intent(in)  :: elements(:,:)

        integer(ik)         :: ierr, connectivity_size, ndims, ielem, ielement_g
        integer(HID_T)      :: element_set_id, element_space_id, grid_id, sid, memspace
        integer(HSIZE_T)    :: dims_rank_two(2), start(2), count(2), dimsm(2)
        logical             :: exists

        integer(ik),  allocatable, target :: var(:,:)
        type(c_ptr)                       :: cp_var


        !
        ! Get number of nodes in element connectivity
        !
        connectivity_size = size(elements,2)

        !
        ! Size element connectivities
        !
        dims_rank_two(1) = nelements_g
        dims_rank_two(2) = connectivity_size


        !
        ! Create a grid-group within the current block domain
        !
        exists = check_link_exists_hdf(dom_id,"Grid")
        if (exists) then
            call h5gopen_f(dom_id, "Grid", grid_id, ierr)
        else
            call h5gcreate_f(dom_id,"Grid", grid_id, ierr)
        end if
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_partition_hdf: h5gcreate_f/h5gopen_f")


        !
        ! Open/Create dataset for element connectivity: element_set_id
        !
        exists = check_link_exists_hdf(grid_id,"Elements")
        if (exists) then
            call h5dopen_f(grid_id,"Elements", element_set_id, ierr, H5P_DEFAULT_F)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_partition_hdf: h5dopen_f")

        else
            call h5screate_simple_f(2, dims_rank_two, element_space_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_partition_hdf: h5screate_simple_f")

            call h5dcreate_f(grid_id, "Elements", H5T_NATIVE_INTEGER, element_space_id, element_set_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_partition_hdf: h5dcreate_f")

            call h5sclose_f(element_space_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_partition_hdf: h5sclose_f")

        end if


        !
        ! Get data space
        !
        call h5dget_space_f(element_set_id, element_space_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL, "set_domain_connectivity_partition_hdf: h5dget_space_f.")



        !
        ! Allocate buffer
        !
        allocate(var(1,connectivity_size), stat=ierr)
        if (ierr /= 0) call AllocationError

        !
        ! Write element connectivities:
        !   - one at a time taking ielement_g from the connectivity since they might 
        !     not be in order after partitioning
        !
        do ielem = 1,size(elements,1)


            !
            ! get domain-global element index
            !
            ielement_g = elements(ielem,2)
            start = [ielement_g-1,1-1]   ! 0-based
            count = [1, connectivity_size]


            !
            ! Select subset of dataspace - sid
            !
            call h5sselect_hyperslab_f(element_space_id, H5S_SELECT_SET_F, start, count, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'set_domain_connectivity_partition_hdf: h5sselect_hyperslab_f')

            !
            ! Create a memory dataspace
            !
            ndims = 2
            dimsm(1) = 1
            dimsm(2) = connectivity_size
            call h5screate_simple_f(ndims,dimsm,memspace,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'set_domain_connectivity_partition_hdf: h5screate_simple_f')


            !
            ! Write modes
            !
            var(1,:) = elements(ielem,:)
            cp_var = c_loc(var(1,1))
            call h5dwrite_f(element_set_id, H5T_NATIVE_INTEGER, cp_var, ierr, memspace, element_space_id)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_partition_hdf: h5dwrite_f")


            call h5sclose_f(memspace,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_partition_hdf: h5sclose_f")



        end do


        !
        ! Close dataset, dataspace, Grid group
        !
        call h5dclose_f(element_set_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_partition_hdf: h5dclose_f")
        call h5sclose_f(element_space_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_partition_hdf: h5sclose_f")
        call h5gclose_f(grid_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_domain_connectivity_partition_hdf: h5gclose_f")



    end subroutine set_domain_connectivity_partition_hdf
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
    !!  TODO: Switch this to read an attribute. That way we can use this in 
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













    !>  Open a boundary patch group on domain 'dom_id' and 'patch'.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   06/01/2017
    !!
    !!
    !--------------------------------------------------------------------------------------
    function open_patch_hdf(dom_id,patch) result(patch_id)
        integer(HID_T), intent(in)  :: dom_id
        character(*),   intent(in)  :: patch

        integer(HID_T)              :: patch_id
        integer(ik)                 :: ierr
        character(:),   allocatable :: msg
        logical                     :: exists

        ! Check exists
        exists = check_link_exists_hdf(dom_id,"Patches/"//"P_"//trim(adjustl(patch)))
        msg    = "open_patch_hdf: Couldn't find bc patch "//trim(adjustl(patch))//" on domain."
        if (.not. exists) call chidg_signal(FATAL,msg)


        ! Open face boundary condition group
        call h5gopen_f(dom_id, "Patches/"//"P_"//trim(adjustl(patch)), patch_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"open_patch_hdf: error opening boundary face group")


    end function open_patch_hdf
    !**************************************************************************************




    !>  Close a boundary patch group, patch_id.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   06/01/2017
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine close_patch_hdf(patch_id)
        integer(HID_T), intent(in)  :: patch_id

        integer(ik) :: ierr

        ! Close face boundary condition group
        call h5gclose_f(patch_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"close_patch_hdf: h5gclose")


    end subroutine close_patch_hdf
    !***************************************************************************************






    !>  Create a bc patch HDF group. Return HDF identifier.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   06/01/2017
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function create_patch_hdf(dom_id,patch_name) result(patch_id)
        integer(HID_T), intent(in)  :: dom_id
        character(*),   intent(in)  :: patch_name

        integer(HID_T)              :: patch_id
        integer(ik)                 :: ierr


        !
        ! Create empty group for boundary condition
        !
        call h5gcreate_f(dom_id,"Patches/"//"P_"//patch_name,patch_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"create_patch_hdf: h5gcreate_f")
        

    end function create_patch_hdf
    !*****************************************************************************************

    


    !>  Return the number of patches on a domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/13/2017
    !!
    !--------------------------------------------------------------------------------------
    function get_npatches_hdf(dom_id) result(npatches)
        integer(HID_T), intent(in)  :: dom_id

        integer(ik)     :: nmembers, ierr, igrp, type, npatches
        character(1024) :: gname


        !  Get number of groups linked to the current bc_face
        call h5gn_members_f(dom_id, "Patches", nmembers, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_npatches_hdf: error h5gn_members_f")

        npatches = 0
        if ( nmembers > 0 ) then
            do igrp = 0,nmembers-1
                ! Get group name
                call h5gget_obj_info_idx_f(dom_id, "Patches", igrp, gname, type, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_npatches_hdf: error h5gget_obj_info_idx_f")

                ! Test if group is a patch: 'P_'
                if (gname(1:2) == 'P_') npatches = npatches + 1
            end do  ! igrp
        end if
        

    end function get_npatches_hdf
    !**************************************************************************************







    !>  Copy bc group association for all patches from file_a to file_b.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/13/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine copy_patches_attributes_hdf(fid_a,fid_b)
        integer(HID_T), intent(in)  :: fid_a
        integer(HID_T), intent(in)  :: fid_b

        integer(HID_T)                  :: dom_id_a, dom_id_b, patch_id_a, patch_id_b
        integer(ik)                     :: idom, ipatch
        character(1024),    allocatable :: domain_names(:), patch_names(:)
        character(:),       allocatable :: patch_bc_group
        logical                         :: exists


        !
        ! Copy Domain patch configuration
        !   For each domain in file_a:
        !       - loop through patches
        !       - open patch and get bc group name
        !       - set bc group name for same patch on file_b
        !       - close current patch
        !       - repeat for next patch
        !
        call write_line("Copying patch boundary group associations...")
        domain_names = get_domain_names_hdf(fid_a) 
        do idom = 1,size(domain_names)

            ! Check if there is a domain of the same name in 'fid_b'
            exists = check_domain_exists_hdf(fid_b,trim(domain_names(idom)))

            if (exists) then
                
                dom_id_a = open_domain_hdf(fid_a, trim(domain_names(idom)))
                dom_id_b = open_domain_hdf(fid_b, trim(domain_names(idom)))

                patch_names = get_patch_names_hdf(dom_id_a)
                do ipatch = 1,size(patch_names)
           
                    ! Check if there is a patch of the same name in fid_b
                    patch_id_a = open_patch_hdf(dom_id_a, trim(patch_names(ipatch)))
                    patch_id_b = open_patch_hdf(dom_id_b, trim(patch_names(ipatch)))

                    patch_bc_group = get_patch_group_hdf(patch_id_a)
                    call set_patch_group_hdf(patch_id_b,patch_bc_group)

                    call close_patch_hdf(patch_id_a)
                    call close_patch_hdf(patch_id_b)

                end do !ipatch

                call close_domain_hdf(dom_id_a)
                call close_domain_hdf(dom_id_b)

            end if !exists

        end do !idom


    end subroutine copy_patches_attributes_hdf
    !*************************************************************************************










    !>  Return a vector of the names for each boundary condition state group.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !----------------------------------------------------------------------------------------
    function get_patch_names_hdf(dom_id) result(patch_names)
        integer(HID_T), intent(in)  :: dom_id

        integer(ik)                     :: nmembers, ierr, igrp, type, npatches, ipatch
        character(1024)                 :: gname
        character(1024),    allocatable :: patch_names(:)


        ! Get number of patches, allocate names
        npatches = get_npatches_hdf(dom_id)
        allocate(patch_names(npatches), stat=ierr)
        if (ierr /= 0) call AllocationError

        !  Get number of groups linked to the "Patches"
        call h5gn_members_f(dom_id, "Patches", nmembers, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_patch_names_hdf: error h5gn_members_f")

        ! Iterate over them and accumulate patch names
        ipatch = 1
        if ( nmembers > 0 ) then
            do igrp = 0,nmembers-1
                ! Get group name
                call h5gget_obj_info_idx_f(dom_id, "Patches", igrp, gname, type, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_patch_names_hdf: error h5gget_obj_info_idx_f")

                ! Test if group is a boundary condition state. 'P_'
                if (gname(1:2) == 'P_') then
                    patch_names(ipatch) = trim(gname(3:))
                    ipatch = ipatch + 1
                end if
            end do  ! igrp
        end if

    end function get_patch_names_hdf
    !***************************************************************************************







    !>  Set patch face connectivity indices for a block boundary.
    !!
    !!  /D_domainname/Patches/P_patchname/Faces
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/15/2016
    !!
    !!  @param[in]  patch_id    HDF identifier for the patch group to be written to
    !!  @param[in]  patch       Face connectivity indices to be set for the patch.
    !!
    !---------------------------------------------------------------------------------------
    subroutine set_patch_hdf(patch_id,patch)
        integer(HID_T), intent(in)  :: patch_id
        integer(ik),    intent(in)  :: patch(:,:)

        integer                     :: ierr
        integer(HID_T)              :: patch_space_id, patch_set_id
        integer(HSIZE_T)            :: dims(2)
        logical                     :: exists
            
        exists = check_link_exists_hdf(patch_id,"Faces")

        !
        ! : already exists, just open data set.
        ! : doesn't exists, create data set
        !
        if (exists) then
            call h5dopen_f(patch_id,"Faces",patch_set_id, ierr, H5P_DEFAULT_F)

        else
            ! Create dataspaces for boundary condition connectivity
            dims(1) = size(patch,1)
            dims(2) = size(patch,2)
            call h5screate_simple_f(2, dims, patch_space_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_patch_hdf: h5screate_simple_f")

            ! Create datasets for boundary condition connectivity
            call h5dcreate_f(patch_id,"Faces", H5T_NATIVE_INTEGER, patch_space_id, patch_set_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_patch_hdf: h5dcreate_f")

        end if


        !
        ! Write patch connectivity
        !
        call h5dwrite_f(patch_set_id, H5T_NATIVE_INTEGER, patch, dims, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_patch_hdf: h5dwrite_f")


        !
        ! Close groups
        !
        call h5dclose_f(patch_set_id, ierr)

        ! If dataset didn't already exists, then a dataspace was created and needs closed also.
        if (.not. exists) call h5sclose_f(patch_space_id, ierr)

    end subroutine set_patch_hdf
    !****************************************************************************************








    !>  Return the patch connectivity information for a boundary condition face group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !----------------------------------------------------------------------------------------
    function get_patch_hdf(patch_id) result(patch)
        use iso_c_binding,  only: c_ptr, c_loc
        integer(HID_T), intent(in)  :: patch_id

        integer(HID_T)      :: faces_did, faces_sid
        integer(HSIZE_T)    :: dims(2), maxdims(2)
        integer(ik)         :: nbcfaces, npts_face, ierr

        integer(ik), allocatable, target    :: patch(:,:)
        type(c_ptr)                         :: patch_p
        

        ! Open Faces patch data
        call h5dopen_f(patch_id, "Faces", faces_did, ierr, H5P_DEFAULT_F)


        !  Get the dataspace id and dimensions
        call h5dget_space_f(faces_did, faces_sid, ierr)
        call h5sget_simple_extent_dims_f(faces_sid, dims, maxdims, ierr)
        nbcfaces  = dims(1)
        npts_face = dims(2)


        ! Read boundary condition patch connectivity
        allocate(patch(nbcfaces,npts_face),stat=ierr)
        if (ierr /= 0) call AllocationError
        patch_p = c_loc(patch(1,1))
        call h5dread_f(faces_did, H5T_NATIVE_INTEGER, patch_p, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_patch_hdf: h5dread_f")


        call h5dclose_f(faces_did,ierr)
        call h5sclose_f(faces_sid,ierr)

    end function get_patch_hdf
    !***************************************************************************************





    

    !>  Set patch attribute 'Boundary State Group'
    !!
    !!  Exists as: P_patchname/Boundary State Group
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine set_patch_group_hdf(patch_id,group)
        integer(HID_T), intent(in)  :: patch_id
        character(*),   intent(in)  :: group

        integer(ik) :: ierr

        ! Set 'Boundary State Group'
        call h5ltset_attribute_string_f(patch_id, ".", "Boundary State Group", trim(group), ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_patch_group_hdf: error setting the attribute 'Boundary State Group'")


    end subroutine set_patch_group_hdf
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
    function get_patch_group_hdf(patch_id) result(group_trim)
        integer(HID_T), intent(in)  :: patch_id

        character(1024)             :: group
        character(:),   allocatable :: group_trim
        integer(ik)                 :: ierr
        integer(ik)                 :: exists


        call check_attribute_exists_hdf(patch_id,"Boundary State Group","Soft Fail",exists)

        ! Get 'Boundary State Group'
        if (exists==0) then
            call h5ltget_attribute_string_f(patch_id, ".", "Boundary State Group", group, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_patch_group_hdf: error setting the attribute 'Boundary State Group'")
            group_trim = trim(group)
        else
            group_trim = 'empty'
        end if
            
    end function get_patch_group_hdf
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
    subroutine create_bc_state_group_hdf(fid,group_name)
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

    end subroutine create_bc_state_group_hdf
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
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
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





    !>  Copy a boundary condition state group from one file to another file.
    !!
    !!  NOTE: if BCSG_ of same name already exists on fid_b, it is removed first.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/14/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine copy_bc_state_groups_hdf(fid_a,fid_b)
        integer(HID_T), intent(in)  :: fid_a
        integer(HID_T), intent(in)  :: fid_b

        integer(ik) :: igroup, ierr

        type(svector_t) :: bc_state_group_names
        type(string_t)  :: group_name


        !
        ! Copy BCSG configuration
        !   - copy entire BCSG groups from file_a to file_b
        !
        call write_line("Copying boundary condition state groups...")
        bc_state_group_names = get_bc_state_group_names_hdf(fid_a)
        do igroup = 1,bc_state_group_names%size()
            ! Get string of current group
            group_name = bc_state_group_names%at(igroup)

            ! Remove group in target if already exists
            call remove_bc_state_group_hdf(fid_b,group_name%get())

            ! Copy group from fid_a to fid_b
            call h5ocopy_f(fid_a,"BCSG_"//trim(group_name%get()),fid_b,"BCSG_"//trim(group_name%get()),ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"copy_bc_state_groups_hdf: error copying boundary state groups.")

        end do !igroup


    end subroutine copy_bc_state_groups_hdf
    !***************************************************************************************






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






    !>  Open a boundary condition state group, return HDF identifier.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/11/2016
    !!
    !---------------------------------------------------------------------------------------
    function open_bc_state_group_hdf(fid,group_name) result(bcgroup_id)
        integer(HID_T),     intent(in)  :: fid
        character(*),       intent(in)  :: group_name

        integer(ik)     :: ierr
        integer(HID_T)  :: bcgroup_id

        call h5gopen_f(fid,"BCSG_"//trim(group_name),bcgroup_id,ierr)
        if (ierr /= 0) call chidg_signal_one(FATAL,"open_bc_state_group_hdf: Error opening boundary condition group",trim(group_name))

    end function open_bc_state_group_hdf
    !***************************************************************************************


    !>  Close a boundary condition state group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/11/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine close_bc_state_group_hdf(bcgroup_id)
        integer(HID_T), intent(in)  :: bcgroup_id

        integer(ik)     :: ierr

        call h5gclose_f(bcgroup_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"close_bc_state_group_hdf: Error closing bc_group")

    end subroutine close_bc_state_group_hdf
    !***************************************************************************************






    !>  Open a boundary condition state, return identifier. "BCS_statename"
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/21/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    function open_bc_state_hdf(group_id,state_name) result(state_id)
        integer(HID_T), intent(in)  :: group_id
        character(*),   intent(in)  :: state_name

        integer(HID_T)  :: state_id
        integer(ik)     :: ierr


        call h5gopen_f(group_id,"BCS_"//trim(state_name), state_id, ierr)
        if (ierr /= 0) call chidg_signal_one(FATAL,"open_bc_state_hdf: Error opening boundary state.",trim(state_name))


    end function open_bc_state_hdf
    !****************************************************************************************



    !>  Close a boundary condition state.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/21/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine close_bc_state_hdf(state_id)
        integer(HID_T), intent(in)  :: state_id

        integer(ik)     :: ierr

        call h5gclose_f(state_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"close_bc_state_hdf: Error closing boundary state.")

    end subroutine close_bc_state_hdf
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


        if ( (bc_state%get_name() == 'empty') .or. &
             (bc_state%get_name() == 'Empty') ) then
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

                    ! Set group 'Family'
                    call set_bc_state_group_family_hdf(bcgroup_id, bc_state%get_family())

                    ! Create a new group for the bc_state_t
                    call h5gcreate_f(bcgroup_id, "BCS_"//bc_state%get_name(), state_id, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"add_bc_state_hdf: error creating new group for bc_state")

                    ! Add bc_state properties to the group that was created
                    call add_bc_properties_hdf(state_id,bc_state)

                    ! Close function group
                    call h5gclose_f(state_id,ierr)


                    !! Get bcgroup family
                    !current_family = get_bc_state_group_family_hdf(bcgroup_id)

                    !!
                    !! Check if new bc_state is of same family
                    !!
                    !if ( (trim(current_family) == 'none') .or. &
                    !     (trim(current_family) == trim(bc_state%get_family())) ) then

                    !    ! Set group 'Family'
                    !    call set_bc_state_group_family_hdf(bcgroup_id, bc_state%get_family())

                    !    ! Create a new group for the bc_state_t
                    !    call h5gcreate_f(bcgroup_id, "BCS_"//bc_state%get_name(), state_id, ierr)
                    !    if (ierr /= 0) call chidg_signal(FATAL,"add_bc_state_hdf: error creating new group for bc_state")

                    !    ! Add bc_state properties to the group that was created
                    !    call add_bc_properties_hdf(state_id,bc_state)

                    !    ! Close function group
                    !    call h5gclose_f(state_id,ierr)

                    !else
                    !    user_msg = "add_bc_state_hdf: Boundary condition state functions in a group &
                    !                must be of the same family"
                    !    call chidg_signal_one(FATAL,user_msg,bc_state%get_family())
                    !end if
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
    !!  /D_domainname/Patches/"face"/BCS_bc_state/BCP_bc_property
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
    !!  /D_domainname/Patches/"face"/BCS_bcstatename/BCP_bcpropertyname/
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









    !>  Set file attribute "/Time Integrator".
    !!
    !!
    !!  @author Mayank Sharma
    !!  @author Nathan A. Wukie
    !!  @date   4/4/2017
    !!
    !!  @param[in]  fid     HDF5 file identifier.
    !!  @param[in]  string  String indicating the time integrator
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_time_integrator_hdf(fid,string)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: string

        integer(ik) :: ierr

        call h5ltset_attribute_string_f(fid,"/","Time Integrator",trim(string), ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_time_integrator_hdf: h5ltset_attribute_string_f")

    end subroutine set_time_integrator_hdf
    !****************************************************************************************









    !>  Return file attribute "/Time Integrator".
    !!
    !!
    !!  @author Mayank Sharma
    !!  @author Nathan A. Wukie
    !!  @date   4/4/2017
    !!
    !!  @param[in]  fid     HDF5 file identifier.
    !!  @result     string  String indicating the time integrator
    !!
    !----------------------------------------------------------------------------------------
    function get_time_integrator_hdf(fid) result(string)
        integer(HID_T), intent(in)  :: fid

        character(:),   allocatable :: string
        character(100)              :: string_buffer

        integer(ik) :: ierr

        call check_attribute_exists_hdf(fid,"Time Integrator")

        call h5ltget_attribute_string_f(fid,"/","Time Integrator",string_buffer, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_time_integrator_hdf: h5ltget_attribute_string_f")

        ! Trim buffer for result
        string = trim(string_buffer)

    end function get_time_integrator_hdf
    !****************************************************************************************







    !>  Set solution times stored in the time.
    !!
    !!  Sets the array:  /Times
    !!
    !!  For time-steady:   size(times)  = 1
    !!  For time-marching: size(times)  = 1
    !!  For time-spectral: size(times) >= 1
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/21/2017
    !!
    !!  TODO: TEST
    !!
    !---------------------------------------------------------------------------------------
    subroutine set_times_hdf(fid,times)
        integer(HID_T), intent(in)  :: fid
        real(rk),       intent(in)  :: times(:)

        integer(HSIZE_T)    :: rank_dims(1), dimsc(1)
        integer(HID_T)      :: space_id, time_id, crp_list
        integer(ik)         :: ntime, ierr, nrank
        logical             :: exists


        ntime     = size(times)
        nrank     = 1
        rank_dims = [ntime]
        if (ntime == 0) call chidg_signal(FATAL,"set_times_hdf: ntimes == 0.")


        !
        ! Modify dataset creation properties, i.e. enable chunking in order to append
        ! dataspace, if needed.
        !
        dimsc = [1]  ! Chunk size

        call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, ierr)
        if (ierr /= 0) call chidg_signal(FATAL, "write_field_domain_hdf: h5pcreate_f error enabling chunking.")

        call h5pset_chunk_f(crp_list, nrank, dimsc, ierr)
        if (ierr /= 0) call chidg_signal(FATAL, "write_field_domain_hdf: h5pset_chunk_f error setting chunk properties.")




        !
        ! : Open 'Times'
        ! : If 'Times' doesn't exist, create new 'Times' data set
        !
        exists = check_link_exists_hdf(fid,'Times')
        if (exists) then
            call h5dopen_f(fid,"Times",time_id, ierr, H5P_DEFAULT_F)
            if (ierr /= 0) call chidg_signal(FATAL,"set_times_hdf: h5dopen_f.")

            ! Extend dataset if necessary
            call h5dset_extent_f(time_id, rank_dims, ierr)
            if (ierr /= 0) call chidg_signal(FATAL, "set_times_hdf: h5dset_extent_f.")

            ! Update existing dataspace ID since it may have been expanded
            call h5dget_space_f(time_id, space_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL, "set_times_hdf: h5dget_space_f.")



        else
            call h5screate_simple_f(nrank, rank_dims, space_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_times_hdf: h5screate_simple_f.")

            call h5dcreate_f(fid, "Times", H5T_NATIVE_DOUBLE, space_id, time_id, ierr, crp_list)
            if (ierr /= 0) call chidg_signal(FATAL,"set_times_hdf: h5dcreate_f.")

        end if

        ! Write times to dataset
        call h5dwrite_f(time_id, H5T_NATIVE_DOUBLE, times, rank_dims, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_times_hdf: h5dwrite_f")

        ! Close datasets
        call h5pclose_f(crp_list,ierr)
        call h5dclose_f(time_id,ierr)
        call h5sclose_f(space_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_times_hdf: h5sclose_f.")
        if (ierr /= 0) call chidg_signal(FATAL,"set_times_hdf: h5dclose_f")


    end subroutine set_times_hdf
    !****************************************************************************************








    !>  Return array of times stored in the file. /Times
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/21/2017
    !!
    !!  TODO: TEST
    !!
    !----------------------------------------------------------------------------------------
    function get_times_hdf(fid) result(times)
        integer(HID_T),             intent(in)      :: fid

        real(rk),   allocatable :: times(:)
        integer(HID_T)          :: time_id, space_id
        integer(HSIZE_T)        :: rank_dims(1), maxdims(3)
        integer(ik)             :: ierr, ntime, rank

        real(rdouble), dimension(:), allocatable, target    :: read_times
        type(c_ptr)                                         :: cp_times

        ! Open dataset
        call h5dopen_f(fid, "Times", time_id, ierr, H5P_DEFAULT_F)

        ! Get dataspace id and dimensions
        call h5dget_space_f(time_id, space_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_times_hdf: h5dget_space_f.")
        call h5sget_simple_extent_dims_f(space_id, rank_dims, maxdims, rank)
        if (rank == -1) call chidg_signal(FATAL,"get_times_hdf: h5sget_simple_extent_dims_f.")
        ntime = rank_dims(1)


        ! Allocat 'times' with storage
        allocate(read_times(ntime), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Read 'Times' to buffer cp_times, (read_times)
        cp_times = c_loc(read_times(1))
        call h5dread_f(time_id, H5T_NATIVE_DOUBLE, cp_times, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_times_hdf: h5dread_f")

        ! Close identifiers
        call h5dclose_f(time_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_times_hdf: h5dclose_f.")
        call h5sclose_f(space_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_times_hdf: h5sclose_f.")
        


        ! Set out-going array using buffer
        times = real(read_times,kind=rk)


    end function get_times_hdf
    !****************************************************************************************







    !>  Set solution times stored in the time.
    !!
    !!  Sets the array:  /Frequencies
    !!
    !!  For time-steady:   size(Frequencies)  = N/A
    !!  For time-marching: size(Frequencies)  = N/A
    !!  For time-spectral: size(Frequencies) >= 1
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/21/2017
    !!
    !!  TODO: TEST
    !!
    !---------------------------------------------------------------------------------------
    subroutine set_frequencies_hdf(fid,freqs)
        integer(HID_T), intent(in)  :: fid
        real(rk),       intent(in)  :: freqs(:)

        integer(HSIZE_T)    :: rank_dims(1)
        integer(HID_T)      :: space_id, freq_id
        integer(ik)         :: nfreq, ierr, nrank
        logical             :: exists

        nfreq     = size(freqs)
        nrank     = 1
        rank_dims = [nfreq]
        if (nfreq == 0) call chidg_signal(FATAL,"set_frequencies_hdf: nfreq == 0.")

        !
        ! : Open 'Frequencies'
        ! : If 'Times' doesn't exists, create new 'Frequencies' data set
        !
        exists = check_link_exists_hdf(fid,'Frequencies')
        if (exists) then
            call h5dopen_f(fid,"Frequencies",freq_id, ierr, H5P_DEFAULT_F)
            if (ierr /= 0) call chidg_signal(FATAL,"set_frequencies_hdf: h5dopen_f.")
        else
            call h5screate_simple_f(nrank, rank_dims, space_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_frequencies_hdf: h5screate_simple_f.")

            call h5dcreate_f(fid, "Frequencies", H5T_NATIVE_DOUBLE, space_id, freq_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_frequencies_hdf: h5dcreate_f.")

            call h5sclose_f(space_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_frequencies_hdf: h5sclose_f.")
        end if

        ! Write coordinates to datasets
        call h5dwrite_f(freq_id, H5T_NATIVE_DOUBLE, freqs, rank_dims, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_frequencies_hdf: h5dwrite_f")

        ! Close datasets
        call h5dclose_f(freq_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_frequencies_hdf: h5dclose_f")


    end subroutine set_frequencies_hdf
    !****************************************************************************************






    !>  Return array of frequencies stored in the file. /Frequencies
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/21/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_frequencies_hdf(fid) result(freqs)
        integer(HID_T),             intent(in)      :: fid

        real(rk),   allocatable :: freqs(:)
        integer(HID_T)          :: freq_id, space_id
        integer(HSIZE_T)        :: rank_dims(1), maxdims(3)
        integer(ik)             :: ierr, nfreq, rank

        real(rdouble), dimension(:), allocatable, target    :: read_freqs
        type(c_ptr)                                         :: cp_freqs

        ! Open dataset
        call h5dopen_f(fid, "Frequencies", freq_id, ierr, H5P_DEFAULT_F)

        ! Get dataspace id and dimensions
        call h5dget_space_f(freq_id, space_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_frequencies_hdf: h5dget_space_f.")
        call h5sget_simple_extent_dims_f(space_id, rank_dims, maxdims, rank)
        if (rank == -1) call chidg_signal(FATAL,"get_frequencies_hdf: h5sget_simple_extend_dims_f.")
        nfreq = rank_dims(1)


        ! Allocate 'read_freqs' with storage
        allocate(read_freqs(nfreq), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Read 'Frequencies' to buffer cp_freqs, (read_freqs)
        cp_freqs = c_loc(read_freqs(1))
        call h5dread_f(freq_id, H5T_NATIVE_DOUBLE, cp_freqs, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_frequencies_hdf: h5dread_f")

        ! Close identifiers
        call h5dclose_f(freq_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_frequencies_hdf: h5dclose_f.")
        call h5sclose_f(space_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_frequencies_hdf: h5sclose_f.")
        


        ! Set out-going array using buffer
        freqs = real(read_freqs,kind=rk)


    end function get_frequencies_hdf
    !****************************************************************************************

















    !>  Given a file identifier, set time step in a hdf5 file
    !!  Used in type_time_integrator_marching
    !!
    !!  @author Mayank Sharma
    !!  @date   4/12/2017
    !!
    !!  @param[in]  fid     HDF file identifier
    !!  @param[in]  dt      Time step
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_time_step_hdf(fid,dt)
        integer(HID_T),     intent(in)  :: fid
        real(rk),           intent(in)  :: dt

        integer(ik)         :: ierr

        call h5ltset_attribute_double_f(fid, "/", "dt", [dt], SIZE_ONE, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_time_step_hdf: Error h5ltget_attribute_double_f")


    end subroutine set_time_step_hdf
    !****************************************************************************************










    !>  Given a file identifier, return time step from a hdf5 file
    !!  Used in type_time_integrator_marching
    !!
    !!  @author Mayank Sharma
    !!  @date   4/12/2017
    !!
    !!  @param[in]  fid     HDF file identifier
    !!
    !----------------------------------------------------------------------------------------
    function get_time_step_hdf(fid) result(dt)
        integer(HID_T),     intent(in)  :: fid
        
        integer                     :: ierr
        real(rk)                    :: dt
        real(rk),   dimension(1)    :: buffer

        call h5ltget_attribute_double_f(fid, "/", "dt", buffer, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_time_step_hdf: h5ltget_attribute_double_f had a &
                                        problem getting time step")

        dt = buffer(1)

    end function get_time_step_hdf
    !***************************************************************************************











    !>  Given a file identifier, set number of time steps in a hdf5 file
    !!  Used in type_time_integrator_marching
    !!
    !!  @author Mayank Sharma
    !!  @date   4/12/2017
    !!
    !!  @param[in]  fid     HDF file identifier
    !!  @param[in]  nsteps  Number of time steps
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_nsteps_hdf(fid,nsteps)
        integer(HID_T),     intent(in)  :: fid
        integer(ik),        intent(in)  :: nsteps

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(fid, "/", "nsteps", [nsteps], SIZE_ONE, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_nsteps_hdf: Error h5ltget_attribute_int_f")


    end subroutine set_nsteps_hdf
    !****************************************************************************************










    !>  Given a file identifier, return number of time steps from a hdf5 file
    !!  Used in type_time_integrator_marching
    !!
    !!  @author Mayank Sharma
    !!  @date   4/12/2017
    !!
    !!  @param[in]  fid     HDF file identifier
    !!
    !----------------------------------------------------------------------------------------
    function get_nsteps_hdf(fid) result(nsteps)
        integer(HID_T),     intent(in)  :: fid
        
        integer                         :: ierr
        integer(ik)                     :: nsteps
        integer(ik),    dimension(1)    :: buffer

        call h5ltget_attribute_int_f(fid, "/", "nsteps", buffer, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_nsteps_hdf: h5ltget_attribute_double_f had a &
                                        problem getting nsteps")

        nsteps = buffer(1)

    end function get_nsteps_hdf
    !***************************************************************************************











    !>  Given a file identifier, set nwrite in a hdf5 file
    !!  Used in type_time_integrator_marching
    !!
    !!  @author Mayank Sharma
    !!  @date   4/12/2017
    !!
    !!  @param[in]  fid     HDF file identifier
    !!  @param[in]  nwrite  Number of steps after which a .h5 file is written
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_nwrite_hdf(fid,nwrite)
        integer(HID_T),     intent(in)  :: fid
        integer(ik),        intent(in)  :: nwrite

        integer(ik)         :: ierr

        call h5ltset_attribute_int_f(fid, "/", "nwrite", [nwrite], SIZE_ONE, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_nwrite_hdf: Error h5ltget_attribute_int_f")


    end subroutine set_nwrite_hdf
    !****************************************************************************************










    !>  Given a file identifier, return nwrite from a hdf5 file
    !!  Used in type_time_integrator_marching
    !!
    !!  @author Mayank Sharma
    !!  @date   4/12/2017
    !!
    !!  @param[in]  fid     HDF file identifier
    !!
    !----------------------------------------------------------------------------------------
    function get_nwrite_hdf(fid) result(nwrite)
        integer(HID_T),     intent(in)  :: fid
        
        integer                         :: ierr
        integer(ik)                     :: nwrite
        integer(ik),    dimension(1)    :: buffer

        call h5ltget_attribute_int_f(fid, "/", "nwrite", buffer, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_nwrite_hdf: h5ltget_attribute_double_f had a &
                                        problem getting nwrite")

        nwrite = buffer(1)

    end function get_nwrite_hdf
    !***************************************************************************************











!    !>  Given a file identifier, set frequency data in a hdf5 file
!    !!  Used in type_time_integrator_spectral
!    !!
!    !!  @author Mayank Sharma
!    !!  @date   4/12/2017
!    !!
!    !!  @param[in]  fid     HDF file identifier
!    !!  @param[in]  freq    Frequency data
!    !!  @param[in]  nfreq   Number of frequencies
!    !1
!    !----------------------------------------------------------------------------------------
!    subroutine set_frequencies_hdf(fid,freq,nfreq)
!        integer(HID_T),     intent(in)  :: fid
!        real(rk),           intent(in)  :: freq(:)
!        integer(HSIZE_T),   intent(in)  :: nfreq
!
!        integer(ik)         :: ierr
!
!        call h5ltset_attribute_double_f(fid, "/", "Frequencies", freq, nfreq, ierr)
!        if (ierr /= 0) call chidg_signal(FATAL,"set_frequencies_hdf: Error h5ltget_attribute_double_f")
!
!
!    end subroutine set_frequencies_hdf
!    !****************************************************************************************
!
!
!
!
!
!
!
!
!
!
!    !>  Given a file identifier, return frequency data from a hdf5 file
!    !!  Used in type_time_integrator_spectral
!    !!
!    !!  @author Mayank Sharma
!    !!  @date   4/12/2017
!    !!
!    !!  @param[in]  fid     HDF file identifier
!    !!  @param[in]  nfreq   Number of frequencies
!    !!
!    !----------------------------------------------------------------------------------------
!    function get_frequencies_hdf(fid,nfreq) result(freq)
!        integer(HID_T),     intent(in)  :: fid
!        integer(HSIZE_T),   intent(in)  :: nfreq
!        
!        integer                 :: ierr
!        real(rk)                :: freq(nfreq)
!
!        call h5ltget_attribute_double_f(fid, "/", "Frequencies", freq, ierr)
!        if (ierr /= 0) call chidg_signal(FATAL,"get_frequencies_hdf: h5ltget_attribute_double_f had a &
!                                        problem getting frequencies")
!
!
!    end function get_frequencies_hdf
!    !***************************************************************************************






!    !>  Given a file identifier, set time level data in a hdf5 file
!    !!  Used in type_time_integrator_spectral
!    !!
!    !!  @author Mayank Sharma
!    !!  @date   4/12/2017
!    !!
!    !!  @param[in]  fid         HDF file identifier
!    !!  @param[in]  time_lev    Time level data
!    !!  @param[in]  ntime       Number of time levels
!    !1
!    !----------------------------------------------------------------------------------------
!    subroutine set_time_levels_hdf(fid,time_lev,ntime)
!        integer(HID_T),     intent(in)  :: fid
!        real(rk),           intent(in)  :: time_lev(:)
!        integer(HSIZE_T),   intent(in)  :: ntime
!
!        integer(ik)         :: ierr
!
!        call h5ltset_attribute_double_f(fid, "/", "Times", time_lev, ntime, ierr)
!        if (ierr /= 0) call chidg_signal(FATAL,"set_time_levels_hdf: Error h5ltget_attribute_double_f")
!
!
!    end subroutine set_time_levels_hdf
!    !****************************************************************************************
!
!
!
!
!
!
!
!
!
!
!    !>  Given a file identifier, return time level data from a hdf5 file
!    !!  Used in type_time_integrator_spectral
!    !!
!    !!  @author Mayank Sharma
!    !!  @date   4/12/2017
!    !!
!    !!  @param[in]  fid     HDF file identifier
!    !!  @param[in]  ntime   Number of time levels
!    !!
!    !----------------------------------------------------------------------------------------
!    function get_time_levels_hdf(fid,ntime) result(time_lev)
!        integer(HID_T),     intent(in)  :: fid
!        integer(HSIZE_T),   intent(in)  :: ntime
!        
!        integer                 :: ierr
!        real(rk)                :: time_lev(ntime)
!
!        call h5ltget_attribute_double_f(fid, "/", "Times", time_lev, ierr)
!        if (ierr /= 0) call chidg_signal(FATAL,"get_time_levels_hdf: h5ltget_attribute_double_f had a &
!                                        problem getting time levels")
!
!
!    end function get_time_levels_hdf
!    !***************************************************************************************











    !>  Create an equation group on the ChiDG HDF file root.
    !!
    !!  If group already exists, no need to do anything, exit routine.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/4/2017
    !!
    !!  @param[in]  fid             HDF file identifier
    !!  @param[in]  group_name      Unique name for the new boundary condition state group.
    !!
    !---------------------------------------------------------------------------------------
    subroutine create_eqn_group_hdf(fid,group_name)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: group_name

        character(:),   allocatable :: user_msg
        integer(HID_T)              :: eqn_id
        integer(ik)                 :: ierr
        logical                     :: group_exists

        !
        ! Check if bc_state group exists
        !
        group_exists = check_link_exists_hdf(fid,"EQN_"//trim(group_name))


        !
        ! Create a new group for the equation set
        !   - if already exists, do nothing.
        !
        if (.not. group_exists) then

            call h5gcreate_f(fid, "EQN_"//trim(group_name), eqn_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'create_eqn_group_hdf: error creating new group for equation set.')
            call h5gclose_f(eqn_id,ierr)

        end if

    end subroutine create_eqn_group_hdf
    !**************************************************************************************







    !>  Remove a equation group from the HDF file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/5/2017
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine remove_eqn_group_hdf(fid,group_name)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: group_name

        integer(HID_T)  :: eqn_id
        integer(ik)     :: ierr
        logical         :: group_exists

        group_exists = check_eqn_group_exists_hdf(fid,trim(group_name))

        if (group_exists) then

            ! Unlink the bc_state group
            call h5gunlink_f(fid,"EQN_"//trim(group_name),ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"remove_eqn_group_hdf: error unlinking eqn group")

        end if

    end subroutine remove_eqn_group_hdf
    !***********************************************************************************************









    !>  Return a vector of the names for each equation group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/5/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_eqn_group_names_hdf(fid) result(eqn_group_names)
        integer(HID_T), intent(in)  :: fid

        integer(ik)     :: nmembers, ierr, igrp, type
        character(1024) :: gname
        type(svector_t) :: eqn_group_names


        !
        !  Get number of groups linked to the file root:
        !
        call h5gn_members_f(fid, ".", nmembers, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_eqn_group_names_hdf: error h5gn_members_f")


        !
        ! Loop through groups and detect "EQN_" groups:
        !
        if ( nmembers > 0 ) then
            do igrp = 0,nmembers-1
                ! Get group name
                call h5gget_obj_info_idx_f(fid, ".", igrp, gname, type, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_eqn_group_names_hdf: error h5gget_obj_info_idx_f")

                ! Test if group is an equation group. 'EQN_'
                if (gname(1:4) == 'EQN_') then
                    call eqn_group_names%push_back(string_t(trim(gname(5:))))
                end if
            end do  ! igrp
        end if


    end function get_eqn_group_names_hdf
    !***************************************************************************************













    !>  Check if an equation set group exists on the file root. 
    !!
    !!  Checks: "/EQ_name"
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/4/2017
    !!
    !--------------------------------------------------------------------------------------
    function check_eqn_group_exists_hdf(fid,group_name) result(exist_status)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: group_name

        integer(ik) :: ierr
        logical     :: exist_status

        ! Check if face contains the bc_state
        call h5lexists_f(fid, "EQN_"//trim(group_name), exist_status, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"check_eqn_group_exists_hdf: Error in call to h5lexists_f")


    end function check_eqn_group_exists_hdf
    !***************************************************************************************








    !>  Remove any equation groups that do not have an associated domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/5/2017
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine prune_eqn_groups_hdf(fid)
        integer(HID_T),     intent(in)  :: fid

        character(1024),    allocatable :: domain_eqns(:)
        type(svector_t)                 :: eqn_groups
        type(string_t)                  :: group_string
        logical                         :: has_domain
        integer(ik)                     :: igroup, idom

        ! Get domain equation sets
        domain_eqns = get_domain_equation_sets_hdf(fid)

        ! Get equation groups
        eqn_groups = get_eqn_group_names_hdf(fid)
        
        !
        ! Loop through equation groups in the file, if they 
        !
        do igroup = 1,eqn_groups%size()
            

            ! Get group name
            group_string = eqn_groups%at(igroup)


            ! Check if any domain is associated with the group
            has_domain = .false.
            do idom = 1,size(domain_eqns)
                has_domain = group_string%get() == trim(domain_eqns(idom))
                if (has_domain) exit
            end do


            ! If group not associated with any domain, remove
            if (.not. has_domain) call remove_eqn_group_hdf(fid,group_string%get())


        end do


    end subroutine prune_eqn_groups_hdf
    !***************************************************************************************











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
        if (ierr /= 0) call chidg_signal_one(FATAL,"check_link_exists_hdf: Error in call to h5lexists_f", trim(linkname) )


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


    !
    !   Mesh Motion
    !

    !   Prescribed Mesh Motion

    !>  Return the number of pmm groups are in the HDF file.
    !!
    !!  Boundary condition state groups: 'PMM_'
    !!
    !!  @author Eric Wolf 
    !!  @date   4/4/2017 
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_npmm_groups_hdf(fid) result(ngroups)
        integer(HID_T)  :: fid
        
        integer(ik)     :: ngroups
        type(svector_t) :: pmm_group_names


        pmm_group_names = get_pmm_group_names_hdf(fid)

        ngroups = pmm_group_names%size()

    end function get_npmm_groups_hdf
    !****************************************************************************************





    
    !>  Return a vector of the names for each pmm group.
    !!
    !!  @author Eric Wolf 
    !!  @date   4/4/2017
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_pmm_group_names_hdf(fid) result(pmm_group_names)
        integer(HID_T), intent(in)  :: fid

        integer(ik)     :: nmembers, ierr, igrp, type
        character(1024) :: gname
        type(svector_t) :: pmm_group_names


        !  Get number of groups linked to the current bc_face
        call h5gn_members_f(fid, ".", nmembers, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_pmm_group_names_hdf: error h5gn_members_f")

        if ( nmembers > 0 ) then
            do igrp = 0,nmembers-1
                ! Get group name
                call h5gget_obj_info_idx_f(fid, ".", igrp, gname, type, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_pmm_group_names_hdf: error h5gget_obj_info_idx_f")

                ! Test if group is a boundary condition state. 'BCSG_'
                if (gname(1:4) == 'PMM_') then
                    call pmm_group_names%push_back(string_t(trim(gname(5:))))
                end if
            end do  ! igrp
        end if

    end function get_pmm_group_names_hdf
    !***************************************************************************************

    !>  Open a PMM group and return HDF group identifier.
    !!
    !!  @author Eric Wolf 
    !!  @date  4/25/2017 
    !!
    !!
    !----------------------------------------------------------------------------------------
    function open_pmm_hdf(fid,pmmname) result(pmmgroup_id)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: pmmname

        integer(HID_T)  :: pmmgroup_id
        integer(ik)     :: ierr
        logical         :: exists


        ! Check exists
        exists = check_link_exists_hdf(fid,"PMM_"//trim(pmmname))
        if (.not. exists) call chidg_signal_one(FATAL,"open_pmm_hdf: Couldn't find PMM in file.","PMM_"//trim(pmmname))


        ! If so, open.
        call h5gopen_f(fid,"PMM_"//trim(pmmname), pmmgroup_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"open_pmm_hdf: Error in h5gopen_f")


    end function open_pmm_hdf
    !****************************************************************************************


    subroutine get_pmm_hdf_test(pmmgroup_id,pmm_name, pmm)
        integer(HID_T), intent(in)  :: pmmgroup_id
        character(*),   intent(in)  :: pmm_name
        class(prescribed_mesh_motion_t), allocatable, intent(inout)  :: pmm


        character(:),       allocatable :: pmmname, pname, oname
        !character(1024)                 :: fname
        character(:),   allocatable     :: fname
        integer(HID_T)                  :: pmm_id, pmmfo_id
        integer(ik)                     :: ierr, iprop, nprop, iopt, noptions
        real(rdouble), dimension(1)     :: buf
        real(rk)                        :: ovalue
        logical                         :: group_exists, function_exists


        !
        !   Prescribed mesh motion group structure
        ! /PMM_pmm_name/...
        !   ATTRIBUTE "Function"    - name of a registered pmmf
        !   /PMMFO_oname/...        - oname is the name of an option for this pmmf
        !       ATTRIBUTE "val"     - value of this option
        !   ...                     - more options

!        ! Open pmmf group 
!               
        ! Get boundary condition name string
        if (pmm_name(1:4) == "PMM_") then
            pmmname = trim(pmm_name(5:))
        else
            pmmname = trim(pmm_name)
        end if


        ! Create boundary condition state and get number of properties

        allocate(pmm, stat=ierr)
        if (ierr /=0) call AllocationError
        call pmm%set_name(pmmname)
       

        ! Read the function name set for the property.
        !fname = ''
        allocate(character(1024) :: fname)
        call h5ltget_attribute_string_f(pmmgroup_id, ".", "Function", fname, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_pmm_hdf: error getting function name.")

        
        ! Set/Create the function for the current property
        call pmm%add_pmmf(trim(fname))

        
        ! Get number of options for the function
        noptions = pmm%pmmf%get_noptions()



        ! Get each option value
        if (noptions>0) then
        do iopt = 1,noptions
            ! Get option name
            oname = pmm%pmmf%get_option_key(iopt)

            ! Get option value from file
            group_exists = check_link_exists_hdf(pmmgroup_id,"PMMFO_"//trim(oname))
            if (group_exists) then
                call h5gopen_f(pmmgroup_id, "PMMFO_"//trim(oname), pmmfo_id, ierr)
                call h5ltget_attribute_double_f(pmmfo_id,".", "val", buf, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_pmm_hdf: error getting option value")
                ovalue = real(buf(1),rk)

                ! Set boundary condition option
                call pmm%pmmf%set_option(trim(oname), ovalue)

                ! Close current property group
                call h5gclose_f(pmmfo_id,ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_pmm_hdf: h5gclose")
            end if
        end do ! iopt
        end if


    end subroutine get_pmm_hdf_test
    !*****************************************************************************************




    

    !>  Given the name of a bc_state on a face, return an initialized bc_state instance.
    !!
    !!  You may consider calling 'get_bc_state_names_hdf' first to get a list of 
    !!  available bc_state's on a face. Then the names could be passed into this routine
    !!  to return the bc_state instance.
    !!
    !!  @author Eric Wolf
    !!  @date   3/30/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_pmm_hdf(pmmgroup_id,pmm_name) result(pmm)
        integer(HID_T), intent(in)  :: pmmgroup_id
        character(*),   intent(in)  :: pmm_name


        class(prescribed_mesh_motion_t),  allocatable :: pmm
        character(:),       allocatable :: pmmname, pname, oname
        character(:), allocatable                 :: fname
        integer(HID_T)                  :: pmm_id, pmmfo_id
        integer(ik)                     :: ierr, iprop, nprop, iopt, noptions
        real(rdouble), dimension(1)     :: buf
        real(rk)                        :: ovalue
        logical                         :: group_exists, function_exists


        !
        !   Prescribed mesh motion group structure
        ! /PMM_pmm_name/...
        !   ATTRIBUTE "Function"    - name of a registered pmmf
        !   /PMMFO_oname/...        - oname is the name of an option for this pmmf
        !       ATTRIBUTE "val"     - value of this option
        !   ...                     - more options

!        ! Open pmmf group 
!               
        ! Get boundary condition name string
        if (pmm_name(1:4) == "PMM_") then
            pmmname = trim(pmm_name(5:))
        else
            pmmname = trim(pmm_name)
        end if


        ! Create boundary condition state and get number of properties

        allocate(pmm, stat=ierr)
        if (ierr /=0) call AllocationError
        call pmm%set_name(pmmname)
       

        ! Read the function name set for the property.
        !fname = ''
        allocate(character(1024) :: fname) 
        call h5ltget_attribute_string_f(pmmgroup_id, ".", "Function", fname, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_pmm_hdf: error getting function name.")

        
        ! Set/Create the function for the current property
        call pmm%add_pmmf(trim(fname))

        
        ! Get number of options for the function
        noptions = pmm%pmmf%get_noptions()



        ! Get each option value
        if (noptions>0) then
        do iopt = 1,noptions
            ! Get option name
            oname = pmm%pmmf%get_option_key(iopt)

            ! Get option value from file
            group_exists = check_link_exists_hdf(pmmgroup_id,"PMMFO_"//trim(oname))
            if (group_exists) then
                call h5gopen_f(pmmgroup_id, "PMMFO_"//trim(oname), pmmfo_id, ierr)
                call h5ltget_attribute_double_f(pmmfo_id,".", "val", buf, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_pmm_hdf: error getting option value")
                ovalue = real(buf(1),rk)

                ! Set boundary condition option
                call pmm%pmmf%set_option(trim(oname), ovalue)

                ! Close current property group
                call h5gclose_f(pmmfo_id,ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_pmm_hdf: h5gclose")
            end if
        end do ! iopt
        end if


    end function get_pmm_hdf
    !*****************************************************************************************



    !>
    !!
    !!  @author Eric Wolf
    !!  @date  3/30/2017 
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine set_pmm_domain_group_hdf(patch_id,group)
        integer(HID_T), intent(in)  :: patch_id
        character(*),   intent(in)  :: group

        integer(ik) :: ierr

        ! Set 'Prescribed Mesh Motion Group'
        call h5ltset_attribute_string_f(patch_id, ".", "Prescribed Mesh Motion Group", trim(group), ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_pmm_domain_group_hdf: error setting the attribute 'Boundary State Group'")


    end subroutine set_pmm_domain_group_hdf
    !***************************************************************************************

    
    !>  Return 'Boundary State Group' attribute for a given patch.
    !!
    !!
    !!  If found, returns the group attribute.
    !!  If not found, returns 'empty'.
    !!
    !!
    !!  @author Eric Wolf 
    !!  @date   3/30/2017
    !!
    !!
    !--------------------------------------------------------------------------------------
    function get_pmm_domain_group_hdf(patch_id) result(group_trim)
        integer(HID_T), intent(in)  :: patch_id

        character(1024)             :: group
        character(:),   allocatable :: group_trim
        integer(ik)                 :: ierr
        integer(ik)                 :: exists


        call check_attribute_exists_hdf(patch_id,"Prescribed Mesh Motion Group","Soft Fail",exists)

        ! Get 'Prescribed Mesh Motion Group'
        if (exists==0) then
            call h5ltget_attribute_string_f(patch_id, ".", "Prescribed Mesh Motion Group", group, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"set_pmm_domain_group_hdf: error setting the attribute 'Prescribed Mesh Motion Group'")
            group_trim = trim(group)
        else
            group_trim = 'empty'
        end if
            
    end function get_pmm_domain_group_hdf
    !***************************************************************************************



    !>  Add a pmm group to the ChiDG HDF file.
    !!
    !!  
    !!      
    !!
    !!  @author Eric Wolf
    !!  @date   4/21/2017 
    !!
    !!  @param[in]  fid             HDF file identifier
    !!  @param[in]  group_name      Unique name for the new boundary condition state group.
    !!
    !----------------------------------------------------------------------------------------
    subroutine create_pmm_group_hdf(fid,group_name,fname)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: group_name
        character(*),   intent(in),optional  :: fname

        character(:),   allocatable :: user_msg
        integer(HID_T)              :: pmmgroup_id
        integer(ik)                 :: ierr
        logical                     :: group_exists


        ! Check if bc_state group exists
        group_exists = check_link_exists_hdf(fid,"PMM_"//trim(group_name))

        user_msg = "create_pmm_group_hdf: Boundary condition state group already exists. &
                    Cannot have two groups with the same name"
!        if (group_exists) call chidg_signal_one(FATAL,user_msg,trim(group_name))


        !
        ! Create a new group for the bc_state_t
        !
        if (.not. group_exists) then
            call h5gcreate_f(fid, "PMM_"//trim(group_name), pmmgroup_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'create_pmm_group_hdf: error creating new group for bc_state.')
        else
            call h5gopen_f(fid, "PMM_"//trim(group_name), pmmgroup_id, ierr)
        end if


        ! Set 'Family'
        if (present(fname)) then
            call h5ltset_attribute_string_f(pmmgroup_id, ".", "Function", fname, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"create_pmm_group_hdf: error setting the attribute 'Function'")
        end if


        call h5gclose_f(pmmgroup_id,ierr)

    end subroutine create_pmm_group_hdf
    !****************************************************************************************

    !>  Sets pmmf name attribute to the ChiDG HDF file.
    !!
    !!  
    !!      
    !!
    !!  @author Eric Wolf
    !!  @date   4/25/2017 
    !!
    !!  @param[in]  fid             HDF file identifier
    !!  @param[in]  pmmname         Unique name for the pmm group.
    !!  @param[in]  fname           Unique name for the pmmf name - must be the name of a registered pmmf!.
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_pmmf_name_hdf(fid,pmmname,fname)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: pmmname 
        character(*),   intent(in)  :: fname

        character(:),   allocatable :: user_msg
        integer(HID_T)              :: pmmgroup_id
        integer(ik)                 :: ierr
        real(rdouble), dimension(1)     :: buf
        logical                     :: group_exists

        group_exists = check_link_exists_hdf(fid,"PMM_"//trim(pmmname))

        user_msg = "set_pmmf_name_hdf: PMM group does not exist. &
                    Create the pmm before trying to set the function name."
        if (.not. group_exists) call chidg_signal_one(FATAL,user_msg,trim(fname))



        pmmgroup_id = open_pmm_hdf(fid, pmmname)

        call h5ltset_attribute_string_f(pmmgroup_id, ".", "Function", fname, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_pmmf_name_hdf: error setting the attribute 'Function'")


        call h5gclose_f(pmmgroup_id,ierr)

    end subroutine set_pmmf_name_hdf
    !****************************************************************************************




    !>  Add a pmm group to the ChiDG HDF file.
    !!
    !!  
    !!      
    !!
    !!  @author Eric Wolf
    !!  @date   4/25/2017 
    !!
    !!  @param[in]  fid             HDF file identifier
    !!  @param[in]  group_name      Unique name for the new boundary condition state group.
    !!
    !----------------------------------------------------------------------------------------
    subroutine create_pmmfo_group_hdf(fid,pmmname,oname,val)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: pmmname 
        character(*),   intent(in)  :: oname
        real(rk),       intent(in),optional  :: val

        character(:),   allocatable :: user_msg
        integer(HSIZE_T)                :: adim
        integer(HID_T)              :: pmmgroup_id, pmmfo_id
        integer(ik)                 :: ierr
        real(rdouble), dimension(1)     :: buf
        logical                     :: group_exists


        pmmgroup_id = open_pmm_hdf(fid, pmmname)

        ! Check if bc_state group exists
        group_exists = check_link_exists_hdf(pmmgroup_id,"PMMFO_"//trim(oname))

        user_msg = "create_pmmfo_group_hdf: PMMF option group already exists. &
                    Cannot have two groups with the same name"
        if (group_exists) call chidg_signal_one(FATAL,user_msg,trim(oname))

        !
        ! Create a new group for the option
        !
        call h5gcreate_f(pmmgroup_id, "PMMFO_"//trim(oname), pmmfo_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'create_pmmfo_group_hdf: error creating new group for PMMF option.')
        if (present(val)) then
            adim = 1
            call h5ltset_attribute_double_f(pmmfo_id, ".", "val", [real(val,rdouble)], adim, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"create_pmmfo_group_hdf: error setting option value")
        end if

        

        call h5gclose_f(pmmfo_id,ierr)
        call h5gclose_f(pmmgroup_id,ierr)

    end subroutine create_pmmfo_group_hdf
    !****************************************************************************************


    !>  Sets pmmf option value attribute to the ChiDG HDF file.
    !!
    !!  
    !!      
    !!
    !!  @author Eric Wolf
    !!  @date   4/25/2017 
    !!
    !!  @param[in]  fid             HDF file identifier
    !!  @param[in]  pmmname         Unique name for the pmm group.
    !!  @param[in]  oname           Unique name for the pmmfo group.
    !!  @param[in]  val             Option value
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_pmmfo_val_hdf(fid,pmmname,oname,val)
        integer(HID_T), intent(in)  :: fid
        character(*),   intent(in)  :: pmmname 
        character(*),   intent(in)  :: oname
        real(rk),       intent(in)  :: val

        character(:),   allocatable :: user_msg
        integer(HSIZE_T)                :: adim
        integer(HID_T)              :: pmmgroup_id, pmmfo_id
        integer(ik)                 :: ierr
        real(rdouble), dimension(1)     :: buf
        logical                     :: group_exists


        pmmgroup_id = open_pmm_hdf(fid, pmmname)

        ! Check if bc_state group exists
        group_exists = check_link_exists_hdf(pmmgroup_id,"PMMFO_"//trim(oname))

        user_msg = "set_pmmfo_val_hdf: PMMF option group does not exist. &
                    Create the pmmfo group before trying to set the value."
        if (.not. group_exists) call chidg_signal_one(FATAL,user_msg,trim(oname))

        !
        ! Create a new group for the option
        !
        call h5gopen_f(pmmgroup_id, "PMMFO_"//trim(oname), pmmfo_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'set_pmmfo_val_hdf: error opening group for PMMF option.')
        adim = 1
        call h5ltset_attribute_double_f(pmmfo_id, ".", "val", [real(val,rdouble)], adim, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_pmmfo_val_hdf: error setting option value")

        

        call h5gclose_f(pmmfo_id,ierr)
        call h5gclose_f(pmmgroup_id,ierr)

    end subroutine set_pmmfo_val_hdf
    !****************************************************************************************




end module mod_hdf_utilities
