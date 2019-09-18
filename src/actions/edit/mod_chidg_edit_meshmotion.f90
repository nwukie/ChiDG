module mod_chidg_edit_meshmotion
#include <messenger.h>
    use mod_kinds,          only: rk, ik, rdouble
    use mod_constants,      only: NFACES, XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use mod_prescribed_mesh_motion_function       
    use mod_radial_basis_function
    use mod_string,         only: string_t
    use type_prescribed_mesh_motion_function,      only: prescribed_mesh_motion_function_t
    use type_svector,       only: svector_t
    use type_function,      only: function_t

    use mod_hdf_utilities
    use hdf5
    use h5lt

    use mod_chidg_edit_printoverview,   only: print_overview, chidg_clear_screen
    use mod_chidg_edit_boundaryconditions,   only: print_bc_states
    implicit none



contains


    !-----------------------------------------------------------------------------------------
    !!
    !!  chidg_edit_meshmotion
    !!
    !!  chidg_edit_boundarycondition_states
    !!  chidg_edit_boundarycondition_state_group
    !!  chidg_edit_boundarycondition_property
    !!
    !!  chidg_edit_boundarycondition_domains
    !!  chidg_edit_boundarycondition_domain_patches
    !!  chidg_edit_boundarycondition_patch
    !!  print_bc_overview
    !!  print_bc_states
    !!  print_patches
    !!  print_bc_state_properties
    !!  print_bc_state_property_options
    !!
    !!
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_meshmotion(fid)
        integer(HID_T),     intent(in)  :: fid

        integer(ik)                     :: ierr, idom_hdf, ndom
        logical                         :: run_bc_edit, print_help
        character(len=:),   allocatable :: command, help_message
        character(1024)                 :: selection

        run_bc_edit = .true.
        print_help  = .false.
        do while ( run_bc_edit )

            ! Refresh display
            call chidg_clear_screen()
            call print_overview(fid)
            call print_mm_overview(fid)


            !
            ! Print command options, accept user input.
            !
            command = "1: Edit Mesh Motion Groups, 2: Edit Mesh Motion Domains, (0: Exit, ?: Help)"
            call write_line(' ')
            call write_line(command,color='blue')



            if (print_help) then
                call write_line(" ")
                call write_line(" ")
                call write_line("Select topic to edit: Mesh Motion Groups or Mesh Motion Domains.", color='red', width=100)
                call write_line(" ")
                call print_mesh_motion_group_help()
                call write_line(" ")
                call print_mesh_motion_domain_help()
                print_help = .false.
            end if


            read(*,'(A1024)', iostat=ierr) selection


            !
            ! Direct user-selection
            !
            select case(trim(selection))
                case('1')
                    call chidg_edit_mm_groups(fid)
                case('2')
                    call chidg_edit_mm_domains(fid)
                case('?')
                    print_help = .true.
                case default
                    run_bc_edit = .false.
            end select 


        end do  ! run_bc_edit

    end subroutine chidg_edit_meshmotion
    !******************************************************************************************


    



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine chidg_edit_mm_groups(fid)
        integer(HID_T), intent(in)  :: fid

        character(:),   allocatable :: command
        character(1024)             :: group_name, group_family
        integer(ik)                 :: selection, ierr
        logical                     :: run_states, edit_group, open_group

        run_states = .true.
        do while(run_states)

            !
            ! Refresh display
            !
            call chidg_clear_screen()
            call print_overview(fid)
            call print_mm_groups(fid,active_topic='Groups')


            command = "1: Create group, 2: Edit group, 3: Remove group, (0 to exit):"
            call write_line(' ')
            call write_line(command,color='blue')

            read(*,'(I8)', iostat=ierr) selection


            select case(selection)
                !
                ! Create group
                !
                case(1)
                    command = "Enter group name: "
                    call write_line(' ')
                    call write_line(command, color='blue')
                    read(*,'(A1024)') group_name

                    call chidg_create_mm_group(fid,group_name)

                !
                ! Edit group
                !
                case(2)
                    edit_group = .true.
                    open_group = .true.
                    do while(edit_group)
                        ! Refresh display
                        call chidg_clear_screen()
                        call print_overview(fid)
                        call print_mm_groups(fid,active_topic='Groups')

                        command = "Enter group name: "
                        call write_line(' ')
                        call write_line(command, color='blue')
                        read(*,'(A1024)') group_name

                        if (trim(group_name) == '') then
                            edit_group = .false.
                            open_group = .false.
                        else
                            open_group = .true.
                        end if

                        if (open_group) call chidg_edit_mm_group(fid,group_name)
                        edit_group = .false.
                    end do
                
                !
                ! Remove group
                !
                case(3)
                    command = "Enter group name: "
                    call write_line(' ')
                    call write_line(command, color='blue')
                    read(*,'(A1024)') group_name

                    call remove_mm_group_hdf(fid,trim(group_name))

                case default
                    run_states = .false.
            end select

        end do ! run_states

    end subroutine chidg_edit_mm_groups
    !******************************************************************************************

 



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine chidg_create_mm_group(fid, group_name)
        integer(HID_T), intent(in)  :: fid
        character(*),       intent(in)  :: group_name

        character(:),   allocatable :: command
        integer(ik)                 :: selection, ierr
        logical                     :: run_states, edit_group, open_group

        run_states = .true.
        do while(run_states)

            !
            ! Refresh display
            !
            call chidg_clear_screen()
            call print_overview(fid)
            call print_mm_groups(fid, active_topic='Groups')


            command = "1: Create RBF Mesh Motion Group, 2: Create Prescribed Mesh Motion Group, (0 to exit):"
            call write_line(' ')
            call write_line(command,color='blue')

            read(*,'(I8)', iostat=ierr) selection


            select case(selection)
                !
                ! Create group
                !
                case(1)
                    call create_mm_group_hdf(fid,group_name, "RBF")
                    !call chidg_edit_rbfmm_group(fid,group_name)

                    run_states = .false.
                case(2)
                    call create_mm_group_hdf(fid,group_name, "PMM")
                    !call chidg_edit_pmm_group(fid,group_name)

                    run_states = .false.
                case default
                    run_states = .false.
            end select

        end do ! run_states

    end subroutine chidg_create_mm_group
    !******************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_mm_group(fid,group_name)
        integer(HID_T),     intent(in)  :: fid
        character(*),       intent(in)  :: group_name

    
        integer(HID_T)              :: bcgroup_id
        integer(ik)                 :: ierr, selection
        character(1024)             :: state_name, bc_string, pname, fname
        character(:),   allocatable :: command
        logical                     :: group_exists, edit_state_group, print_bcs, set_bc, found_bc_state, &
                                       bc_state_exists, remove_bc, get_property, property_exists, edit_property
        class(bc_state_t), allocatable    :: bc_state


        group_exists = check_link_exists_hdf(fid,"MM_"//trim(group_name))

        if (group_exists) then


            call h5gopen_f(fid,"MM_"//trim(group_name),bcgroup_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_boundarycondition_state_group: h5gopen_f")


            !Get the MM Family and call the appropriate edit procedure 
            call h5ltget_attribute_string_f(bcgroup_id, ".", "Family", fname, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_mm_group: h5ltget_attribute_string_f")

            ! Close group
            call h5gclose_f(bcgroup_id,ierr)

            select case(trim(fname))

            case("RBF")

                call chidg_edit_rbfmm_group(fid,group_name)


            case("PMM")

                call chidg_edit_pmm_group(fid,group_name)

            case default

                call chidg_signal(FATAL,"chidg_edit_mm_group: MM Family not recognized! Must be one of: 'RBF', or 'PMM'.")

            end select

        end if !group_exists

    end subroutine chidg_edit_mm_group
    !******************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_rbfmm_group(fid,group_name)
        integer(HID_T),     intent(in)  :: fid
        character(*),       intent(in)  :: group_name

    
        integer(HID_T)              :: bcgroup_id
        integer(ik)                 :: ierr, selection
        character(1024)             :: state_name, bc_string, pname
        character(:),   allocatable :: command
        logical                     :: group_exists, edit_state_group, print_bcs, set_bc, found_bc_state, &
                                       bc_state_exists, remove_bc, get_property, property_exists, edit_property
        class(bc_state_t), allocatable    :: bc_state


        group_exists = check_link_exists_hdf(fid,"MM_"//trim(group_name))

        if (group_exists) then


            call h5gopen_f(fid,"MM_"//trim(group_name),bcgroup_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_rbfmm_group: h5gopen_f")

            edit_state_group = .true.
            do while(edit_state_group)

                !
                ! Refresh display
                !
                call chidg_clear_screen()
                call print_overview(fid)
                call print_mm_groups(fid,active_topic='Groups',active_group=trim(group_name))

                command = "1: Edit RBF parameters, 2: Edit RBF sources, (0 to exit):"
                call write_line(' ')
                call write_line(command,color='blue')
                read(*,'(I8)', iostat=ierr) selection


                select case(selection)
                    case(1)

                        call chidg_edit_rbf_parameters(fid, bcgroup_id, group_name)
                    !
                    ! Remove bc_state
                    !
                    case(2)

                        call chidg_edit_rbf_sources(fid, bcgroup_id, group_name)

                    case default
                        edit_state_group = .false.
                end select




            end do !edit_state_group


            ! Close group
            call h5gclose_f(bcgroup_id,ierr)

        end if !group_exists

    end subroutine chidg_edit_rbfmm_group
    !******************************************************************************************

    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_rbf_parameters(fid,sourcegroup_id, sourcegroup_name)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: sourcegroup_id
        character(*),       intent(in)  :: sourcegroup_name
    
        integer(HID_T)              :: bcgroup_id, pmmfgroup_id
        integer(ik)                 :: ierr, selection, confirm
        character(1024)             :: state_name, bc_string, pname
        character(:),   allocatable :: command
        logical                     :: group_exists, edit_state_group, print_bcs, set_bc, found_bc_state, &
                                       bc_state_exists, remove_bc, get_property, property_exists, edit_property, pmmf_exists
        class(prescribed_mesh_motion_function_t), allocatable    :: bc_state

        type(svector_t)             :: pmmf_group_names
        type(string_t)              :: group_name

        edit_state_group = .true.
        do while(edit_state_group)

            !
            ! Refresh display
            !
            call chidg_clear_screen()
            call print_overview(fid)
            call print_mm_groups(fid, active_topic='Groups', active_group=trim(sourcegroup_name))
            call print_rbf_parameters(sourcegroup_id)
            !call print_rbf_overview(fid, sourcegroup_id)

            command = "1: Edit base RBF parameters, 2: Edit explicit RBF parameters, (0 to exit):"
            call write_line(' ')
            call write_line(command,color='blue')
            read(*,'(I8)', iostat=ierr) selection


            select case(selection)
                !
                ! Add bc_state
                !
                case(1)

                    call chidg_edit_rbf_base_parameters(fid, sourcegroup_id, sourcegroup_name)
                !
                ! Remove bc_state
                !
                case(2)

                    call chidg_edit_rbf_explicit_parameters(fid, sourcegroup_id, sourcegroup_name)
                !
                case default
                    edit_state_group = .false.
            end select




        end do !edit_state_group



    end subroutine chidg_edit_rbf_parameters
    !******************************************************************************************


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_rbf_base_parameters(fid,sourcegroup_id, sourcegroup_name)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: sourcegroup_id
        character(*),       intent(in)  :: sourcegroup_name
    
        integer(HID_T)              :: bcgroup_id, pmmfgroup_id
        integer(ik)                 :: ierr, selection, confirm
        character(1024)             :: state_name, bc_string, pname, rbftype
        character(:),   allocatable :: command
        logical                     :: group_exists, edit_state_group, print_bcs, set_bc, found_bc_state, &
                                       bc_state_exists, remove_bc, get_property, property_exists, edit_property, pmmf_exists

        type(svector_t)             :: pmmf_group_names
        type(string_t)              :: group_name

        integer(ik)                 :: read_status
        logical                     :: read_val, print_rbfs, set_rbf_type
        real(rk)                    :: radius(3)
        edit_state_group = .true.
        do while(edit_state_group)

            !
            ! Refresh display
            !
            call chidg_clear_screen()
            call print_overview(fid)
            call print_mm_groups(fid, active_topic='Groups', active_group=trim(sourcegroup_name))
            call print_rbf_parameters(sourcegroup_id, active_rbf='Base')
            !call print_rbf_overview(fid, sourcegroup_id)

            command = "1: Edit base RBF Type, 2: Edit base RBF Radius, (0 to exit):"
            call write_line(' ')
            call write_line(command,color='blue')
            read(*,'(I8)', iostat=ierr) selection


            select case(selection)
                case(1)
                    print_rbfs = .false.
                    set_rbf_type = .true.
                    do while (set_rbf_type)
                        if (print_rbfs) then
                            call list_radial_basis_functions()
                        end if

                        if (.not. found_bc_state) then
                            call write_line("The RBF specified wasn't found in the ChiDG library. &
                                             Maybe check that the string was input correctly?",color='red')
                        end if


                        command = "Enter RBF Base Type (? to list RBFs):"
                        call write_line(' ')
                        call write_line(command,color='blue')

                        read(*,'(A1024)', iostat=ierr) rbftype
                        ! Check for user exit
                        set_rbf_type = (trim(rbftype) /= "")


                        if (set_rbf_type) then
                            if ( trim(rbftype) == '?' ) then
                                print_rbfs = .true.
                            else
                                ! Call routine to set boundary condition in hdf file.
                                found_bc_state = check_rbf_registered(rbftype)
                                if (found_bc_state) then
                                    call set_rbf_type_hdf(sourcegroup_id, trim(rbftype))
                                    set_rbf_type = .false.
                                end if

                            end if
                        end if
                    end do



                case(2)
                    command = "Enter Radius-1:"
                    call write_line(' ')
                    call write_line(command,color='blue')
                    read_val = .true.
                    do while(read_val)
                        read(*,*,iostat=read_status) radius(1)
                        if (read_status == 0) then
                            read_val = .false.
                        end if
                    end do
        
                    command = "Enter Radius-2:"
                    call write_line(' ')
                    call write_line(command,color='blue')
                    read_val=.true.
                    do while(read_val)
                        read(*,*,iostat=read_status) radius(2)
                        if (read_status == 0) then
                            read_val = .false.
                        end if
                    end do
 
                    command = "Enter Radius-3:"
                    call write_line(' ')
                    call write_line(command,color='blue')
                    read_val=.true.
                    do while(read_val)
                        read(*,*,iostat=read_status) radius(3)
                        if (read_status == 0) then
                            read_val = .false.
                        end if
                    end do
 
                    call set_rbf_base_radius_hdf(sourcegroup_id, radius)

                case default
                    edit_state_group = .false.
            end select




        end do !edit_state_group



    end subroutine chidg_edit_rbf_base_parameters
    !******************************************************************************************

    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_rbf_explicit_parameters(fid,sourcegroup_id, sourcegroup_name)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: sourcegroup_id
        character(*),       intent(in)  :: sourcegroup_name
    
        integer(HID_T)              :: bcgroup_id, pmmfgroup_id
        integer(ik)                 :: ierr, selection, confirm
        character(1024)             :: state_name, bc_string, pname, rbftype
        character(:),   allocatable :: command
        logical                     :: group_exists, edit_state_group, print_bcs, set_bc, found_bc_state, &
                                       bc_state_exists, remove_bc, get_property, property_exists, edit_property, pmmf_exists

        type(svector_t)             :: pmmf_group_names
        type(string_t)              :: group_name

        integer(ik)                 :: read_status
        logical                     :: read_val, print_rbfs, set_rbf_type
        real(rk)                    :: base_fraction
        edit_state_group = .true.
        do while(edit_state_group)

            !
            ! Refresh display
            !
            call chidg_clear_screen()
            call print_overview(fid)
            call print_mm_groups(fid, active_topic='Groups', active_group=trim(sourcegroup_name))
            call print_rbf_parameters(sourcegroup_id, active_rbf='Explicit')
            !call print_rbf_overview(fid, sourcegroup_id)

            command = "1: Edit explicit RBF Type, 2: Edit base fraction, (0 to exit):"
            call write_line(' ')
            call write_line(command,color='blue')
            read(*,'(I8)', iostat=ierr) selection


            select case(selection)
                case(1)
                    print_rbfs = .false.
                    set_rbf_type = .true.
                    do while (set_rbf_type)
                        if (print_rbfs) then
                            call list_radial_basis_functions()
                        end if

                        if (.not. found_bc_state) then
                            call write_line("The RBF specified wasn't found in the ChiDG library. &
                                             Maybe check that the string was input correctly?",color='red')
                        end if


                        command = "Enter RBF Explicit Type (? to list RBFs):"
                        call write_line(' ')
                        call write_line(command,color='blue')

                        read(*,'(A1024)', iostat=ierr) rbftype
                        ! Check for user exit
                        set_rbf_type = (trim(rbftype) /= "")


                        if (set_rbf_type) then
                            if ( trim(rbftype) == '?' ) then
                                print_rbfs = .true.
                            else
                                ! Call routine to set boundary condition in hdf file.
                                found_bc_state = check_rbf_registered(rbftype)
                                if (found_bc_state) then
                                    call set_rbf_type_explicit_hdf(sourcegroup_id, trim(rbftype))
                                    set_rbf_type = .false.
                                end if

                            end if
                        end if
                    end do



                case(2)
                    command = "Enter base fraction:"
                    call write_line(' ')
                    call write_line(command,color='blue')
                    read_val = .true.
                    do while(read_val)
                        read(*,*,iostat=read_status) base_fraction
                        if (read_status == 0) then
                            read_val = .false.
                        end if
                    end do
        
                    call set_rbf_base_fraction_hdf(sourcegroup_id, base_fraction)

                case default
                    edit_state_group = .false.
            end select




        end do !edit_state_group



    end subroutine chidg_edit_rbf_explicit_parameters
    !******************************************************************************************

    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine chidg_edit_rbf_sources(fid, mm_group_id, mm_group_name)
        integer(HID_T), intent(in)  :: fid
        integer(HID_T), intent(in)  :: mm_group_id
        character(*),   intent(in)  :: mm_group_name

        character(:),   allocatable :: command
        character(1024)             :: group_name, group_family
        integer(ik)                 :: selection, ierr
        logical                     :: run_states, edit_group, open_group

        run_states = .true.
        do while(run_states)

            !
            ! Refresh display
            !
            call chidg_clear_screen()
            call print_overview(fid)
            call print_mm_groups(fid,active_topic='Groups', active_group=trim(mm_group_name))
            call print_rbf_sources(mm_group_id,active_topic = "Sources")


            command = "1: Create RBF source, 2: Edit RBF source, 3: Remove RBF source, (0 to exit):"
            call write_line(' ')
            call write_line(command,color='blue')

            read(*,'(I8)', iostat=ierr) selection


            select case(selection)
                !
                ! Create group
                !
                case(1)
                    command = "Enter RBF source name: "
                    call write_line(' ')
                    call write_line(command, color='blue')
                    read(*,'(A1024)') group_name

                    call create_rbf_src_group_hdf(mm_group_id, group_name)

                !
                ! Edit group
                !
                case(2)
                    edit_group = .true.
                    open_group = .true.
                    do while(edit_group)
                        ! Refresh display
                        call chidg_clear_screen()
                        call print_overview(fid)
                        call print_mm_groups(fid,active_topic='Groups', active_group = trim(mm_group_name))
                        call print_rbf_sources(mm_group_id,active_topic = "Sources")

                        command = "Enter RBF source name: "
                        call write_line(' ')
                        call write_line(command, color='blue')
                        read(*,'(A1024)') group_name

                        if (trim(group_name) == '') then
                            edit_group = .false.
                            open_group = .false.
                        else
                            open_group = .true.
                            edit_group = .false.
                        end if

                        if (open_group) call chidg_edit_rbf_source_group(fid, mm_group_id, mm_group_name, group_name)
                    end do
                
                !
                ! Remove group
                !
                case(3)
                    command = "Enter RBF source name: "
                    call write_line(' ')
                    call write_line(command, color='blue')
                    read(*,'(A1024)') group_name

                    call remove_rbf_src_group_hdf(mm_group_id,trim(group_name))

                case default
                    run_states = .false.
            end select

        end do ! run_states

    end subroutine chidg_edit_rbf_sources
    !******************************************************************************************


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_rbf_source_group(fid, mm_group_id, mm_group_name, group_name)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: mm_group_id
        character(*),       intent(in)  :: mm_group_name
        character(*),       intent(in)  :: group_name

    
        integer(HID_T)              :: bcgroup_id
        integer(ik)                 :: ierr, selection
        character(1024)             :: state_name, bc_string, pname
        character(:),   allocatable :: command
        logical                     :: group_exists, edit_state_group, print_bcs, set_bc, found_bc_state, &
                                       bc_state_exists, remove_bc, get_property, property_exists, edit_property


        group_exists = check_link_exists_hdf(mm_group_id,"RBF_SRC_"//trim(group_name))

        if (group_exists) then


            call h5gopen_f(mm_group_id,"RBF_SRC_"//trim(group_name),bcgroup_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_rbf_source_group: h5gopen_f")

            edit_state_group = .true.
            do while(edit_state_group)

                !
                ! Refresh display
                !
                call chidg_clear_screen()
                call print_overview(fid)
                call print_mm_groups(fid,active_topic='Groups',active_group=trim(mm_group_name))
                call print_rbf_sources(mm_group_id,active_topic = "Sources",active_source=trim(group_name))


                command = "1: Edit Patch Name, 2: Edit RBF MM Driver, (0 to exit):"
                call write_line(' ')
                call write_line(command,color='blue')
                read(*,'(I8)', iostat=ierr) selection


                select case(selection)
                    case(1)

                        call chidg_edit_rbf_source_patch_name(fid, mm_group_id, bcgroup_id, mm_group_name, group_name)
                    !
                    ! Remove bc_state
                    !
                    case(2)

                        call chidg_edit_rbf_mm_driver(fid, mm_group_id, bcgroup_id, mm_group_name, group_name)

                    case default
                        edit_state_group = .false.
                end select




            end do !edit_state_group


            ! Close group
            call h5gclose_f(bcgroup_id,ierr)

        end if !group_exists

    end subroutine chidg_edit_rbf_source_group
    !******************************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_rbf_source_patch_name(fid, mm_group_id, sourcegroup_id, mm_group_name, group_name)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: mm_group_id
        integer(HID_T),     intent(in)  :: sourcegroup_id
        character(*),       intent(in)  :: mm_group_name
        character(*),       intent(in)  :: group_name

        integer(HID_T)              :: patch_group, domain_id
        integer(ik)                 :: ierr, selection
        character(:),   allocatable :: dname, dname_trim, command
        character(1024)             :: group
        logical                     :: run_edit_domain, run, set_group



        !
        ! edit / boundary condition / domain
        !
        run_edit_domain = .true.
        do while ( run_edit_domain )

            !
            ! Refresh display
            !
            call chidg_clear_screen()
            call print_overview(fid)
            call print_bc_states(fid)
            call print_rbf_sources(mm_group_id,active_topic = "Sources", active_source=trim(group_name))


            !
            ! Print command options, accept user selection.
            !
            command = "1. Add or replace Patch Name, 2. Clear Patch Name, 0. Exit:"
            call write_line(' ')
            call write_line(command,color='blue')
            ierr = 1
            do while ( ierr /= 0 )
                read(*,'(I8)', iostat=ierr) selection
                if ( ierr /= 0 )  call write_line("Invalid input: expecting an integer index.")
            end do


            !
            ! Select option
            !
            select case(selection)
                case(1)
                    
                    call write_line('Enter the Patch Name to set: ')
                    read(*,'(A1024)', iostat=ierr) group
                    if ( ierr /= 0 )  call write_line("Invalid input: expecting Mesh Motion Group name.")
                    
                    ! Check for user exit
                    set_group = (trim(group) /= "")

                    if (set_group) then
                        call set_rbf_src_patchname_hdf(sourcegroup_id,trim(group))
                    end if

                case(2)

                    call set_rbf_src_patchname_hdf(sourcegroup_id,'empty')

                case default
                    run_edit_domain = .false.
            end select




        end do  ! run_bc_edit

    end subroutine chidg_edit_rbf_source_patch_name
    !******************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine chidg_edit_rbf_mm_driver(fid, mm_group_id, sourcegroup_id, mm_group_name, group_name)
        integer(HID_T), intent(in)  :: fid
        integer(HID_T), intent(in)  :: mm_group_id
        integer(HID_T), intent(in)  :: sourcegroup_id
        character(*),       intent(in)  :: mm_group_name
        character(*),       intent(in)  :: group_name

        integer(HID_T)              :: driver_id
        character(:),   allocatable :: command
        integer(ik)                 :: selection, ierr, confirm
        logical                     :: run_states, edit_group, open_group, mm_driver_exists, remove_bc

        run_states = .true.
        do while(run_states)

            !
            ! Refresh display
            !
            call chidg_clear_screen()
            call print_overview(fid)
            call print_mm_groups(fid, active_topic='Groups', active_group = trim(mm_group_name))
            call print_rbf_sources(mm_group_id,active_topic = "Sources", active_source=trim(group_name))


            command = "1: Add (or replace) RBF MM Driver, 2. Remove RBF MM Driver, 3. Edit RBF MM Driver, (0 to exit):"
            call write_line(' ')
            call write_line(command,color='blue')

            read(*,'(I8)', iostat=ierr) selection
            select case(selection)
                !
                ! Add bc_state
                !
                case(1)
                    mm_driver_exists = check_any_rbf_mm_driver_exists_hdf(sourcegroup_id)

                    mm_driver_exists = .false.
                    if (mm_driver_exists) then
                        call remove_rbf_mm_driver_hdf(sourcegroup_id)
                    end if
                    call chidg_create_rbf_mm_driver_group(fid, mm_group_id, sourcegroup_id, mm_group_name, group_name)

                !
                ! Remove bc_state
                !
                case(2)
                    remove_bc = .true.
                    do while (remove_bc)


                        ! Get boundary condition state from user to remove
                        command = "1: Confirm, 2: Cancel"
                        call write_line(' ')
                        call write_line(command,color='blue')
                        read(*,'(I8)', iostat=ierr) confirm 

                        select case(confirm)
                            case(1)
                                mm_driver_exists = check_rbf_mm_driver_exists_hdf(sourcegroup_id)

                                if (mm_driver_exists) then
                                   call remove_rbf_mm_driver_hdf(sourcegroup_id)
                                end if

                                remove_bc = .false.

                            case(2)

                                remove_bc = .false.

                            case default

                                remove_bc = .false.

                        end select
                    end do !remove_bc


                !
                ! Edit property
                !
                case(3)
                    mm_driver_exists = check_rbf_mm_driver_exists_hdf(sourcegroup_id)

                    if (mm_driver_exists) then
                        driver_id = open_rbf_mm_driver_hdf(sourcegroup_id)
                        call chidg_edit_rbf_mm_driver_group(fid, mm_group_id, sourcegroup_id, mm_group_name, group_name)
                        call close_rbf_mm_driver_hdf(driver_id)
                    end if


                case default
                    run_states = .false.
            end select

        end do ! run_states

    end subroutine chidg_edit_rbf_mm_driver
    !******************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine chidg_create_rbf_mm_driver_group(fid, mm_group_id, sourcegroup_id, mm_group_name, source_name)
        integer(HID_T), intent(in)  :: fid
        integer(HID_T), intent(in)  :: mm_group_id
        integer(HID_T), intent(in)  :: sourcegroup_id
        character(*),       intent(in)  :: mm_group_name
        character(*),       intent(in)  :: source_name

        character(:),   allocatable :: command
        integer(ik)                 :: selection, ierr
        logical                     :: run_states, edit_group, open_group

        run_states = .true.
        do while(run_states)

            !
            ! Refresh display
            !
            call chidg_clear_screen()
            call print_overview(fid)
            call print_mm_groups(fid, active_topic='Groups', active_group=trim(mm_group_name))
            call print_rbf_sources(mm_group_id,active_topic = "Sources", active_source=trim(source_name))

            command = "1: Create PMMF-RBF MM Driver, (0 to exit):"
            call write_line(' ')
            call write_line(command,color='blue')

            read(*,'(I8)', iostat=ierr) selection

            !
            ! NOTE: To implement further rbf_mm_driver_t extensions, add corresponding cases and specialized procedure calls below
            !
            select case(selection)
                !
                ! Create group
                !
                case(1)
                    call create_rbf_mm_driver_hdf(sourcegroup_id, "PMM")

                    run_states = .false.
                case default
                    run_states = .false.
            end select

        end do ! run_states

    end subroutine chidg_create_rbf_mm_driver_group
    !******************************************************************************************

!>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_rbf_mm_driver_group(fid, mm_group_id, sourcegroup_id, mm_group_name, source_name)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: mm_group_id
        integer(HID_T),     intent(in)  :: sourcegroup_id
        character(*),       intent(in)  :: mm_group_name
        character(*),       intent(in)  :: source_name

    
        integer(HID_T)              :: bcgroup_id
        integer(ik)                 :: ierr, selection
        character(1024)             :: state_name, bc_string, pname, fname
        character(:),   allocatable :: command
        logical                     :: group_exists, edit_state_group, print_bcs, set_bc, found_bc_state, &
                                       bc_state_exists, remove_bc, get_property, property_exists, edit_property
        class(bc_state_t), allocatable    :: bc_state


        group_exists = check_link_exists_hdf(sourcegroup_id,"RBF Driver")

        if (group_exists) then


            call h5gopen_f(sourcegroup_id,"RBF Driver",bcgroup_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_boundarycondition_state_group: h5gopen_f")


            !Get the MM Family and call the appropriate edit procedure 
            call h5ltget_attribute_string_f(bcgroup_id, ".", "Family", fname, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_mm_group: h5ltget_attribute_string_f")

            ! Close group
            call h5gclose_f(bcgroup_id,ierr)
            
            !
            ! NOTE: To implement further rbf_mm_driver_t extensions, add corresponding cases and specialized procedure calls below
            !

            select case(trim(fname))

            case("PMM")

                call chidg_edit_pmm_rbf_mm_driver_group(fid,mm_group_id,sourcegroup_id,mm_group_name,source_name)

            case default

                call chidg_signal(FATAL,"chidg_edit_mm_group: MM Family not recognized! Must be one of: 'RBF', or 'PMM'.")

            end select

        end if !group_exists

    end subroutine chidg_edit_rbf_mm_driver_group
    !******************************************************************************************

    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_pmm_rbf_mm_driver_group(fid,mm_group_id, sourcegroup_id,mm_group_name, source_name)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: mm_group_id
        integer(HID_T),     intent(in)  :: sourcegroup_id
        character(*),       intent(in)  :: mm_group_name
        character(*),       intent(in)  :: source_name

    
        integer(HID_T)              :: bcgroup_id
        integer(ik)                 :: ierr, selection
        character(1024)             :: state_name, bc_string, pname
        character(:),   allocatable :: command
        logical                     :: group_exists, edit_state_group, print_bcs, set_bc, found_bc_state, &
                                       bc_state_exists, remove_bc, get_property, property_exists, edit_property


        group_exists = check_link_exists_hdf(sourcegroup_id,"RBF Driver")

        if (group_exists) then


            call h5gopen_f(sourcegroup_id,"RBF Driver",bcgroup_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_pmm_group: h5gopen_f")


            call chidg_edit_pmmf_group(fid,bcgroup_id, mm_group_name)

            ! Close group
            call h5gclose_f(bcgroup_id,ierr)

        end if !group_exists

    end subroutine chidg_edit_pmm_rbf_mm_driver_group
    !******************************************************************************************












    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_pmm_group(fid,group_name)
        integer(HID_T),     intent(in)  :: fid
        character(*),       intent(in)  :: group_name

    
        integer(HID_T)              :: bcgroup_id
        integer(ik)                 :: ierr, selection
        character(1024)             :: state_name, bc_string, pname
        character(:),   allocatable :: command
        logical                     :: group_exists, edit_state_group, print_bcs, set_bc, found_bc_state, &
                                       bc_state_exists, remove_bc, get_property, property_exists, edit_property
        class(prescribed_mesh_motion_function_t), allocatable    :: bc_state


        group_exists = check_link_exists_hdf(fid,"MM_"//trim(group_name))

        if (group_exists) then


            call h5gopen_f(fid,"MM_"//trim(group_name),bcgroup_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_pmm_group: h5gopen_f")


            call chidg_edit_pmmf_group(fid,bcgroup_id, group_name)

            ! Close group
            call h5gclose_f(bcgroup_id,ierr)

        end if !group_exists

    end subroutine chidg_edit_pmm_group
    !******************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_pmmf_group(fid,sourcegroup_id, mm_group_name)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: sourcegroup_id
        character(*),       intent(in)  :: mm_group_name
    
        integer(HID_T)              :: bcgroup_id, pmmfgroup_id
        integer(ik)                 :: ierr, selection, confirm
        character(1024)             :: state_name, bc_string, pname
        character(:),   allocatable :: command
        logical                     :: group_exists, edit_state_group, print_bcs, set_bc, found_bc_state, &
                                       bc_state_exists, remove_bc, get_property, property_exists, edit_property, pmmf_exists
        class(prescribed_mesh_motion_function_t), allocatable    :: bc_state

        type(svector_t)             :: pmmf_group_names
        type(string_t)              :: group_name

        edit_state_group = .true.
        do while(edit_state_group)

            !
            ! Refresh display
            !
            call chidg_clear_screen()
            call print_overview(fid)
            call print_mm_groups(fid, active_topic='Groups', active_group=trim(mm_group_name))

            command = "1: Add (or replace) PMMF, 2: Remove PMMF, 3: Edit PMMF Option, (0 to exit):"
            call write_line(' ')
            call write_line(command,color='blue')
            read(*,'(I8)', iostat=ierr) selection


            select case(selection)
                !
                ! Add bc_state
                !
                case(1)
                    set_bc    = .true.
                    print_bcs = .false.
                    found_bc_state = .true.
                    do while (set_bc)


                        if (print_bcs) then
                            call list_prescribed_mesh_motion_functions()
                        end if

                        if (.not. found_bc_state) then
                            call write_line("The PMMF specified wasn't found in the ChiDG library. &
                                             Maybe check that the string was input correctly?",color='red')
                        end if

                        ! Get bc_state string from user
                        command = "Enter prescribed mesh motion function(? to list): "
                        call write_line(' ')
                        call write_line(command,color='blue')
                        read(*,"(A1024)") bc_string

                        ! Check for user exit
                        set_bc = (trim(bc_string) /= "")


                        if (set_bc) then
                            if ( trim(bc_string) == '?' ) then
                                print_bcs = .true.
                            else

                                ! Call routine to remove boundary condition state
                                bc_state_exists = check_any_pmmf_exists_hdf(sourcegroup_id)

                                if (bc_state_exists) then
                                    call remove_any_pmmf_hdf(sourcegroup_id)
                                end if

   
                                ! Call routine to set boundary condition in hdf file.
                                found_bc_state = check_pmmf_registered(bc_string)
                                if (found_bc_state) then
                                    call create_prescribed_mesh_motion_function(bc_state,bc_string)
                                    call add_pmmf_hdf(sourcegroup_id,bc_state)
                                    set_bc = .false.
                                end if

                            end if
                        end if


                    end do ! set_bc


                !
                ! Remove bc_state
                !
                case(2)
                    remove_bc = .true.
                    do while (remove_bc)


                        ! Get boundary condition state from user to remove
                        command = "1: Confirm, 2: Cancel"
                        call write_line(' ')
                        call write_line(command,color='blue')
                        read(*,'(I8)', iostat=ierr) confirm 

                        select case(confirm)
                            case(1)

                                ! Call routine to remove boundary condition state
                                bc_state_exists = check_any_pmmf_exists_hdf(sourcegroup_id)

                                if (bc_state_exists) then
                                    call remove_any_pmmf_hdf(sourcegroup_id)
                                end if

                                remove_bc = .false.

                            case(2)

                                remove_bc = .false.

                            case default

                                remove_bc = .false.

                        end select
                    end do !remove_bc


                !
                ! Edit property
                !
                case(3)

                    pmmf_exists = check_any_pmmf_exists_hdf(sourcegroup_id)
                    if (pmmf_exists) then
                        pmmf_group_names = get_pmmf_names_hdf(sourcegroup_id)
                        group_name = pmmf_group_names%at(1)

                        pmmfgroup_id = open_pmmf_hdf(sourcegroup_id, trim(group_name%get()))
                        edit_property = .true.
                    else
                        edit_property=.false.
                    end if

                    do while (edit_property)
                        get_property = .true.

                        do while ( get_property )

                            ! Refresh display
                            call chidg_clear_screen()
                            call write_line(trim(group_name%get()))
                            call print_overview(fid)
                            call print_mm_groups(fid, active_topic='Groups', active_group=trim(mm_group_name))
                            call print_pmmf_options(pmmfgroup_id)

                            ! Get property
                            call write_line(" ")
                            call write_line("Enter PMMF option: ",color='blue')
                            read(*,"(A1024)") pname

                            ! Check for user exit
                            if (trim(pname) == "") then
                                get_property = .false.
                                edit_property = .false.
                            end if

                            if (get_property) then
                                ! Check property exists
                                property_exists = check_pmmfo_exists_hdf(pmmfgroup_id,trim(pname))
                                if ( property_exists ) get_property = .false.

                                if (property_exists) then
                                    get_property = .false.
                                else
                                    call write_line("Invalid option", color='blue')
                                end if

                            else
                                property_exists = .false.
                            end if
                        end do 

                        if (property_exists) then
                            call chidg_edit_pmmf_option(fid,pmmfgroup_id,trim(pname))
                        else
                            call write_line("Invalid option", color='blue')
                        end if
                        
                    end do

                    if (pmmf_exists) call close_pmmf_hdf(pmmfgroup_id)


                case default
                    edit_state_group = .false.
            end select




        end do !edit_state_group



    end subroutine chidg_edit_pmmf_group
    !******************************************************************************************














    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_mm_domains(fid)
        integer(HID_T),     intent(in)  :: fid

        integer(ik)                         :: ierr, idom_hdf, ndom
        character(len=:),       allocatable :: command
        character(1024)                     :: domain_name
        logical                             :: run_domains, open_domain, domain_exists


        run_domains = .true.
        do while(run_domains)

            !
            ! Refresh display
            !
            call chidg_clear_screen()
            call print_overview(fid)
            call print_mm_overview(fid,active_topic='Domains')


            command = "Enter a domain to edit(Enter nothing to exit):"
            call write_line(' ')
            call write_line(command,color='blue')
            read(*,'(A1024)') domain_name

            
            domain_exists = check_domain_exists_hdf(fid,trim(domain_name))


            if (trim(domain_name) == '') then
                run_domains = .false.
                open_domain = .false.
            else if ( (trim(domain_name) /= '') .and. (.not. domain_exists) ) then
                run_domains = .true.
                open_domain = .false.
            else
                open_domain = .true.
            end if



            if (open_domain) call chidg_edit_mm_domain_group(fid,trim(domain_name))


        end do

    end subroutine chidg_edit_mm_domains
    !******************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_mm_domain_group(fid,domain_name)
        integer(HID_T),     intent(in)  :: fid
        character(*),       intent(in)  :: domain_name

        integer(HID_T)              :: patch_group, domain_id
        integer(ik)                 :: ierr, selection
        character(:),   allocatable :: dname, dname_trim, command
        character(1024)             :: group
        logical                     :: run_edit_domain, run, set_group

        ! Open domain
        domain_id = open_domain_hdf(fid,domain_name)


        !
        ! edit / boundary condition / domain
        !
        run_edit_domain = .true.
        do while ( run_edit_domain )

            !
            ! Refresh display
            !
            call chidg_clear_screen()
            call print_overview(fid,domain_id)
            call print_mm_overview(fid, active_topic = 'Domains', active_domain=domain_id)
            dname_trim = trim(adjustl(domain_name)) 


            !
            ! Print command options, accept user selection.
            !
            command = "1. Add or replace MM group, 2. Clear MM group, 0. Exit:"
            call write_line(' ')
            call write_line(command,color='blue')
            ierr = 1
            do while ( ierr /= 0 )
                read(*,'(I8)', iostat=ierr) selection
                if ( ierr /= 0 )  call write_line("Invalid input: expecting an integer index.")
            end do


            !
            ! Select option
            !
            select case(selection)
                case(1)
                    
                    call write_line('Enter the Mesh Motion Group to set: ')
                    read(*,'(A1024)', iostat=ierr) group
                    if ( ierr /= 0 )  call write_line("Invalid input: expecting Mesh Motion Group name.")
                    
                    ! Check for user exit
                    set_group = (trim(group) /= "")

                    if (set_group) then
                        call set_mm_domain_group_hdf(domain_id,trim(group))
                    end if

                case(2)

                    call set_mm_domain_group_hdf(domain_id,'empty')

                case default
                    run_edit_domain = .false.
            end select




        end do  ! run_bc_edit


        !
        ! Close groups
        !

        call close_domain_hdf(domain_id)

    end subroutine chidg_edit_mm_domain_group
    !******************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!  @note   Modified to include bc_states
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_pmmf_option(fid, pmmfgroup_id, oname)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: pmmfgroup_id
        character(*),       intent(in)  :: oname

        integer(HID_T)                  :: bcprop, bc_state
        integer(HSIZE_T)                :: adim
        character(len=:),   allocatable :: command
        character(len=1024)             :: option, function_string, gname
        logical                         :: option_exists, run, property_found, property_exists, set_fcn, list_fcns, read_val, run_function
        real(rk)                        :: val
        integer                         :: ierr, type, igrp, iop, nmembers, read_status
        type(svector_t)                 :: bc_state_strings
        type(string_t)                  :: string
        class(function_t),  allocatable :: func


        ! Refresh display
        call chidg_clear_screen()
        call print_overview(fid)
        !call print_mm_overview(fid,active_topic='Groups')
        call print_mm_groups(fid, active_topic='Groups')
        call print_pmmf_options(pmmfgroup_id)

       


        !
        ! Edit option loop
        !
        list_fcns = .false.
        run = .true.
        do while (run)

            ! Refresh display
            call chidg_clear_screen()
            call print_overview(fid)
            call print_mm_groups(fid, active_topic='Groups')
            !call print_bc_overview(fid,active_topic='States',active_group=trim(group_name))
            call print_pmmf_options(pmmfgroup_id, active_option=trim(oname))



            if (run) then

                ! Check option exists
                call h5aexists_f(pmmfgroup_id, trim(oname), option_exists, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_pmmf_option: error opening option attribute")


                !
                ! Handle option_exists
                !
                if ( option_exists ) then


                    command = "Set option value: "
                    call write_line(command, color='blue')
                    read_val = .true.
                    do while(read_val)
                        read(*,*,iostat=read_status) val
                        if (read_status == 0) then
                            read_val = .false.
                        end if
                    end do
        


                    !
                    ! Set option
                    !
                    adim = 1
                    call h5ltset_attribute_double_f(pmmfgroup_id, ".", trim(oname), [real(val,rdouble)], adim, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_boundarycondition_property: error setting option value")

                    run = .false.



                else

                    call write_line("Invalid option", color='blue')
                    run = .false.

                end if

            end if

        end do ! run


    end subroutine chidg_edit_pmmf_option
    !******************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_mm_overview(fid,active_topic,active_group,active_domain,active_face)
        integer(HID_T),     intent(in)              :: fid
        character(*),       intent(in), optional    :: active_topic
        character(*),       intent(in), optional    :: active_group
        integer(HID_T),     intent(in), optional    :: active_domain
        integer(ik),        intent(in), optional    :: active_face


        call write_line(' ')
        call write_line(' ')
        call print_mm_groups(fid,active_topic,active_group)

        call write_line(' ')
        call write_line(' ')
        call write_line(' ')

        call print_mm_domains(fid, active_topic,active_domain)
        call write_line(' ')

    end subroutine print_mm_overview
    !******************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_mm_groups(fid,active_topic,active_group)
        integer(HID_T), intent(in)              :: fid
        character(*),   intent(in), optional    :: active_topic
        character(*),   intent(in), optional    :: active_group

        integer(HID_T)              :: bcgroup_id
        integer(ik)                 :: igroup, istate, ierr
        type(svector_t)             :: bc_state_groups, bc_state_names, mm_group_names
        type(string_t)              :: group_name, bc_state_name
        character(:),   allocatable :: group_family, color
        logical                     :: bold



        color = 'none'
        if (present(active_topic)) then
            if (trim(active_topic) == 'Groups') then
                color = 'blue'
            end if
        end if


        !
        ! Write boundary state overview header
        !
        if (present(active_topic)) then 
            if (trim(active_topic) == 'Groups') then
                call write_line(":[Mesh Motion Groups]",color=color,bold=.true.)
            else
                call write_line(":Mesh Motion Groups")
            end if
        else
            call write_line(":Mesh Motion Groups")
        end if



        ! Write state group header
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
        call write_line("Group","Family", columns=.true.,column_width=25,delimiter=':',color=color)
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
            
        !
        ! Get boundary condition information
        !
        mm_group_names = get_mm_group_names_hdf(fid)


        do igroup = 1,mm_group_names%size()

            group_name = mm_group_names%at(igroup)
            bcgroup_id = open_mm_group_hdf(fid,group_name%get())

            group_family = get_mm_family_hdf(bcgroup_id)

            if (group_name%get() == trim(active_group)) then
                color = 'blue'
                bold  = .true.
            else
                color = 'none'
                bold = .false.
            end if


            call add_to_line(trim(group_name%get()), columns=.true., column_width=25,color=color,bold=bold)
            call add_to_line(trim(group_family), columns=.true., column_width=25,color=color,bold=bold)


            call close_mm_group_hdf(bcgroup_id)

            call send_line()

        end do !istate

        call write_line(" ")

    end subroutine print_mm_groups
    !******************************************************************************************




!>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_mm_domains(fid,active_topic,active_domain)
        integer(HID_T),     intent(in)              :: fid
        character(*),       intent(in), optional    :: active_topic
        integer(HID_T),     intent(in), optional    :: active_domain

        integer(HID_T)              :: bcgroup_id, dom_id
        integer(ik)                 :: igroup, istate, ierr, idom_hdf
        type(svector_t)             :: bc_state_groups, bc_state_names
        type(string_t)              :: group_name, bc_state_name
        character(:),   allocatable :: group_family, color
        logical                     :: bold

        character(len=1024), allocatable :: dnames(:)
        character(len=1024)         :: domain_name, active_domain_name, mm_group


        color = 'none'
        if (present(active_topic)) then
            if (trim(active_topic) == 'Domains') then
                color = 'blue'
            end if
        end if


        !
        ! Write boundary state overview header
        !
        if (present(active_topic)) then 
            if (trim(active_topic) == 'Domains') then
                call write_line(":[Mesh Motion Domains]",color=color,bold=.true.)
            else
                call write_line(":Mesh Motion Domains")
            end if
        else
            call write_line(":Mesh Motion Domains")
        end if



        ! Write state group header
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
        call write_line("Domain","MM Group", columns=.true.,column_width=25,delimiter=':',color=color)
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
            

        !
        ! Get boundary condition information
        !
        dnames = get_domain_names_hdf(fid)

        do idom_hdf = 1,size(dnames)

            !
            ! Get domain name and bc patch groups associated with idom_hdf
            !
            domain_name = dnames(idom_hdf)
            
            dom_id = open_domain_hdf(fid,trim(domain_name))

            mm_group = get_mm_domain_group_hdf(dom_id)
            call close_domain_hdf(dom_id)


            if (present(active_domain)) then

                active_domain_name = get_domain_name_hdf(active_domain)
            
                !
                ! Print active domain
                !
                if (trim(active_domain_name) == trim(domain_name)) then


                    ! Need to add information individually to selectively color the entries.
                    call add_to_line(domain_name, columns=.True., column_width=25, color='pink')

                    call add_to_line(trim(mm_group), columns=.True., column_width=25, color='pink')

                    call send_line()



                !
                ! Print non-active domains
                !
                else

                    call add_to_line( domain_name, columns=.True., column_width=25)
                    call add_to_line( trim(mm_group), columns=.True., column_width=25)
                    call send_line()

                end if

            else
                !
                ! Print domain info if no active domain is present
                !
                call add_to_line( domain_name, columns=.True., column_width=25)
                call add_to_line( trim(mm_group), columns=.True., column_width=25)
                call send_line()

            end if

        end do

        call write_line(" ")

    end subroutine print_mm_domains
    !******************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_rbf_parameters(sourcegroup_id, active_rbf)
        use iso_fortran_env
        integer(HID_T),     intent(in)              :: sourcegroup_id  
        character(*),       intent(in), optional    :: active_rbf


        integer(HID_T)                          :: rbfgroup_id  
        integer(ik)                             :: nattr, ierr, lines_printed
        integer(HSIZE_T)                        :: iattr, idx
        character(len=1024),    allocatable     :: option_keys(:)
        character(len=1024)                     :: fcn
        character(:),           allocatable     :: color
        real(rk),               allocatable     :: option_vals(:)
        type(h5o_info_t),       target          :: h5_info
        real(rdouble),   dimension(1)                :: buf

        logical :: rbf_exists, rbf_exp_exists
        character(1024) :: rbfname, rbfexpname
        real(rk)        :: radius(3), base_fraction

        rbf_exists = check_link_exists_hdf(sourcegroup_id, "RBF Function")
        rbf_exp_exists = check_link_exists_hdf(sourcegroup_id, "RBF Function Explicit")

        if (rbf_exists) then
            rbfgroup_id = open_rbf_function_hdf(sourcegroup_id)
 
            ! RBF Type (mandatory)
            ! e.g. "tps" (not compactly supported), "wc6", etc 
            call h5ltget_attribute_string_f(rbfgroup_id, ".", "RBF Type", rbfname, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_rbf_mm_hdf: error getting RBF Type.")

            ! Get the RBF radius
            ! To allow for ellipsoidal RBFs, we have three radii along each coordinate direction
            ! This might be useful for 1D/2D problems that are slender in 2/1 directions,
            ! where an isotropic RBF might lead to  an ill-conditioned linear system because 
            ! the isotropic RBF would be close to constant in  the slender directions.
            ! 
            ! To avoid this issue, take the radius in the slender direction to be on the order
            ! of the mesh length in the slender direction.
            
            call h5ltget_attribute_double_f(rbfgroup_id, ".", "RBF Radius-1", buf, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_rbf_mm_hdf: error getting RBF Radius-1.")
            radius(1) = real(buf(1), rk)

            call h5ltget_attribute_double_f(rbfgroup_id, ".", "RBF Radius-2", buf, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_rbf_mm_hdf: error getting RBF Radius-2.")
            radius(2) = real(buf(1), rk)

            call h5ltget_attribute_double_f(rbfgroup_id, ".", "RBF Radius-3", buf, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_rbf_mm_hdf: error getting RBF Radius-3.")
            radius(3) = real(buf(1), rk)

            call h5gclose_f(rbfgroup_id,ierr)
        else

            !Otherwise, use default values
            rbfname = "empty"
            radius = ZERO

        end if

        if (rbf_exp_exists) then
            rbfgroup_id = open_rbf_function_exp_hdf(sourcegroup_id)
            ! If so, read in parameters.
            !

            ! RBF Explicit type (optional, must be compactly supported!)
            ! e.g. "wc6"
            call h5ltget_attribute_string_f(rbfgroup_id, ".", "RBF Type Explicit", rbfexpname, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_rbf_mm_hdf: error getting RBF Base.")


            ! Get the fraction of source nodes to be used as base nodes (expect ~0.05-0.1)
            ! Optional
            call h5ltget_attribute_double_f(rbfgroup_id, ".", "RBF Base Fraction", buf, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"print_rbf_parameters: error getting RBF Base Fraction.")
            base_fraction = real(buf(1), rk)

            call h5gclose_f(rbfgroup_id,ierr)
            
            ! Set up the RBF

        else

            rbfexpname = 'empty'
            base_fraction = ONE
        end if

        color = 'none'
        if (present(active_rbf)) then
            if (trim(active_rbf) == 'Base') then
                color = 'blue'
            end if
        end if

        if (present(active_rbf)) then 
            if (trim(active_rbf) == 'Base') then
                call write_line(":[RBF Base Parameters]",color=color,bold=.true.)
            else
                call write_line(":RBF Base Parameters")
            end if
        else
            call write_line(":RBF Base Parameters")
        end if




        ! Write state group header
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
        call write_line("RBF Type","Radius-1","Radius-2","Radius-3",columns=.true.,column_width=25,delimiter=':',color=color)
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
         
        call write_line(trim(rbfname),radius(1),radius(2),radius(3),columns=.true.,column_width=25,delimiter=':',color=color)
        call write_line(" ")

        color = 'none'
        if (present(active_rbf)) then
            if (trim(active_rbf) == 'Explicit') then
                color = 'blue'
            end if
        end if

        if (present(active_rbf)) then 
            if (trim(active_rbf) == 'Explicit') then
                call write_line(":[RBF Explicit Parameters]",color=color,bold=.true.)
            else
                call write_line(":RBF Explicit Parameters")
            end if
        else
            call write_line(":RBF Explicit Parameters")
        end if



        ! Write state group header
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
        call write_line("RBF Type Explicit","Base Fraction",columns=.true.,column_width=25,delimiter=':',color=color)
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
        
        call write_line(trim(rbfexpname),base_fraction,columns=.true.,column_width=25,delimiter=':',color=color)
        call write_line(" ")


    end subroutine print_rbf_parameters
    !******************************************************************************************


!>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_rbf_sources(mm_group_id,active_topic,active_source)
        integer(HID_T), intent(in)              :: mm_group_id
        character(*),   intent(in), optional    :: active_topic
        character(*),   intent(in), optional    :: active_source

        integer(HID_T)              :: bcgroup_id
        integer(ik)                 :: igroup, istate, ierr
        type(svector_t)             :: bc_state_groups, bc_state_names, mm_group_names
        type(string_t)              :: group_name, bc_state_name
        character(:),   allocatable :: group_family, color, patch_name, source_comment
        logical                     :: bold


        !
        ! Write boundary state overview header
        !
        if (present(active_topic)) then 
            if (trim(active_topic) == 'Sources') then
                call write_line(":[RBF Sources]",color=color,bold=.true.)
            else
                call write_line(":RBF Sources")
            end if
        else
            call write_line(":RBF Sources")
        end if



        ! Write state group header
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
        call write_line("Source Name","Patch Group Name", "Driver Family", "Comment",columns=.true.,column_width=25,delimiter=':',color=color)
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
            
        !
        ! Get boundary condition information
        !
        mm_group_names = get_rbf_src_group_names_hdf(mm_group_id)


        do igroup = 1,mm_group_names%size()

            group_name = mm_group_names%at(igroup)
            bcgroup_id = open_rbf_src_hdf(mm_group_id,group_name%get())

            patch_name = get_rbf_src_patch_name_hdf(bcgroup_id)
            group_family = get_rbf_src_family_hdf(bcgroup_id)
            source_comment = get_rbf_src_comment_hdf(bcgroup_id)

            if (group_name%get() == trim(active_source)) then
                color = 'blue'
                bold  = .true.
            else
                color = 'none'
                bold = .false.
            end if


            call add_to_line(trim(group_name%get()), columns=.true., column_width=25,color=color,bold=bold)
            call add_to_line(trim(patch_name), columns=.true., column_width=25,color=color,bold=bold)
            call add_to_line(trim(group_family), columns=.true., column_width=25,color=color,bold=bold)
            call add_to_line(trim(source_comment), columns=.true., column_width=25,color=color,bold=bold)


            call close_rbf_src_hdf(bcgroup_id)

            call send_line()

        end do !istate

        call write_line(" ")

    end subroutine print_rbf_sources
    !******************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_pmmf_options(pmmfgroup_id, active_option)
        use iso_fortran_env
        integer(HID_T),     intent(in)              :: pmmfgroup_id  
        character(*),       intent(in), optional    :: active_option


        integer(ik)                             :: nattr, ierr, lines_printed
        integer(HSIZE_T)                        :: iattr, idx
        character(len=1024),    allocatable     :: option_keys(:)
        character(len=1024)                     :: fcn
        character(:),           allocatable     :: color
        real(rk),               allocatable     :: option_vals(:)
        type(h5o_info_t),       target          :: h5_info
        real(rdouble),   dimension(1)                :: buf

        call h5ltget_attribute_string_f(pmmfgroup_id, ".", "Function", fcn, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"print_pmmf_options: error getting PMMF name.")

        !
        ! Get number of attributes attached to the group id
        !
        call h5oget_info_f(pmmfgroup_id, h5_info, ierr)
        nattr = h5_info%num_attrs
        if (ierr /= 0) call chidg_signal(FATAL,"print_pmmf_options: error getting current number of attributes.")


        allocate(option_keys(nattr), option_vals(nattr), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Gather any existing attributes and their values
        !
        if ( nattr > 0 ) then
            idx = 0
            do iattr = 1,nattr
                idx = iattr - 1
                call h5aget_name_by_idx_f(pmmfgroup_id, ".", H5_INDEX_CRT_ORDER_F, H5_ITER_NATIVE_F, idx, option_keys(iattr), ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"print_property_options: error reading attribute name")

                if ( trim(option_keys(iattr)) /= 'Function' ) then  ! don't read function. not a floating point attribute.
                    call h5ltget_attribute_double_f(pmmfgroup_id, ".", option_keys(iattr), buf, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"print_property_options: error reading attribute value")
                    option_vals(iattr) = real(buf(1),rk)
                end if
            end do

        end if ! nattr


        color = 'none'
        ! Write state group header
        call write_line(" PMM Function Name: ", fcn)
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
        call write_line("PMMF Option","Value", columns=.true.,column_width=25,delimiter=':',color=color)
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
         

        !
        ! Print the gathered property options. Except the function attribute, which exists for all
        ! properties and was printed above.
        !
        lines_printed = 0
        do iattr = 1,nattr
            if (present(active_option)) then
                    if (trim(option_keys(iattr))==trim(active_option)) then
                    color = 'blue'
                else
                    color='none'
                end if
            end if

            if ( (iattr == 1) .and. (trim(option_keys(iattr)) /= "Function") ) then
                call write_line(option_keys(iattr), option_vals(iattr), columns=.true., column_width = 25,color=color)
                lines_printed = lines_printed + 1

            else if ( (iattr == 1) .and. (trim(option_keys(iattr)) == "Function") ) then
                ! Don't print this line. We don't want to print "Function"
            else
                if (lines_printed == 0) then
                    call write_line(option_keys(iattr), option_vals(iattr), columns=.true., column_width = 25,color=color)
                else
                    call write_line(option_keys(iattr), option_vals(iattr), columns=.true., column_width = 25,color=color)
                end if
                lines_printed = lines_printed + 1
            end if

        end do

        call write_line(" ")



    end subroutine print_pmmf_options
    !******************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/9/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_mesh_motion_group_help()
        call write_line("Boundary State Groups:", color='blue', width=100)
        call write_line("These are formed to impose a condition or set of &
                         conditions on a region. For example, one might do this for the RANS &
                         equations; creating a group called 'Inlet' and adding a Total Inlet &
                         boundary state in addition to an &
                         inlet state for the turbulence working variables. In order for these &
                         state functions to be applied on a boundary, they must be associated &
                         with a boundary condition patch, or set of boundary condition patches. &
                         The association can be made by editing a patch and choosing the association.", color='red', width=100)

    end subroutine print_mesh_motion_group_help
    !******************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/9/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_mesh_motion_domain_help()
        call write_line("Boundary Patches:", color='blue', width=100)
        call write_line("These specify regions over which Boundary State Groups &
                         are applied. For a given patch, one can associate a Boundary State Group &
                         that will be applied over the patch region. Currently, the patches are &
                         generated in the conversion process for block-structured grids. The &
                         original faces of the block-structured grids are the patches defining &
                         the boundary geometry for each domain. By editing a patch, one can choose &
                         which Boundary State Group the patch is associated with. It might be helpful &
                         to first, enter the Boundary State Group menu and create a group. Visually, &
                         it should make a bit more sense when a group exists first. Then a given patch &
                         can be associated with the group.", color='red', width=100)

    end subroutine print_mesh_motion_domain_help
    !******************************************************************************************









end module mod_chidg_edit_meshmotion
