module mod_chidg_edit_boundaryconditions
#include <messenger.h>
    use mod_kinds,          only: rk, ik, rdouble
    use mod_constants,      only: NFACES, XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use mod_function,       only: create_function, list_functions
    use mod_bc,             only: create_bc, list_bcs, check_bc_state_registered
    use mod_string,         only: string_t
    use type_bc_state,      only: bc_state_t
    use type_svector,       only: svector_t
    use type_function,      only: function_t
    use hdf5
    use h5lt

    use mod_chidg_edit_printoverview,   only: print_overview
    use mod_hdf_utilities,     only: get_ndomains_hdf, get_domain_names_hdf,                        &
                                     get_domain_name_hdf, get_bc_state_names_hdf,                   &
                                     delete_group_attributes_hdf,                                   &
                                     add_bc_state_hdf, set_bc_property_function_hdf,                &
                                     create_bc_group_hdf, remove_bc_state_group_hdf,                &
                                     get_bc_state_group_names_hdf, get_bc_state_group_family_hdf,   &
                                     get_bc_patch_group_hdf, set_bc_patch_group_hdf,                &
                                     check_bc_property_exists_hdf, remove_bc_state_hdf,             &
                                     check_bc_state_exists_hdf, check_link_exists_hdf,              &
                                     open_bc_group_hdf, close_bc_group_hdf, get_bc_state_group_family_hdf, &
                                     open_domain_hdf, close_domain_hdf, check_domain_exists_hdf
    implicit none



contains


    !-----------------------------------------------------------------------------------------
    !!
    !!  chidg_edit_boundaryconditions
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
    !!  print_bc_patches
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
    subroutine chidg_edit_boundaryconditions(fid)
        integer(HID_T),     intent(in)  :: fid

        integer(ik)                     :: ierr, idom_hdf, ndom
        logical                         :: run_bc_edit, print_help
        character(len=:),   allocatable :: command, help_message
        character(1024)                 :: selection

        run_bc_edit = .true.
        print_help  = .false.
        do while ( run_bc_edit )

            ! Refresh display
            call execute_command_line("clear")
            call print_overview(fid)
            call print_bc_overview(fid)


            !
            ! Print command options, accept user input.
            !
            command = "1: Edit Boundary State Groups, 2: Edit Boundary Patches, (0: Exit, ?: Help)"
            call write_line(' ')
            call write_line(command,color='blue')



            if (print_help) then
                call write_line(" ")
                call write_line(" ")
                call write_line("Select topic to edit: Boundary State Groups or Boundary Patches.", color='red', width=100)
                call write_line(" ")
                call print_boundary_state_group_help()
                call write_line(" ")
                call print_boundary_patch_help()
                print_help = .false.
            end if


            read(*,'(A1024)', iostat=ierr) selection


            !
            ! Direct user-selection
            !
            select case(trim(selection))
                case('1')
                    call chidg_edit_boundarycondition_states(fid)
                case('2')
                    call chidg_edit_boundarycondition_domains(fid)
                case('?')
                    print_help = .true.
                case default
                    run_bc_edit = .false.
            end select 


        end do  ! run_bc_edit

    end subroutine chidg_edit_boundaryconditions
    !******************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine chidg_edit_boundarycondition_states(fid)
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
            call execute_command_line("clear")
            call print_overview(fid)
            call print_bc_overview(fid,active_topic='States')


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
                    call write_line(command, color='blue')
                    read(*,'(A1024)') group_name

                    call create_bc_group_hdf(fid,group_name)
                    call chidg_edit_boundarycondition_state_group(fid,group_name)

                !
                ! Edit group
                !
                case(2)
                    edit_group = .true.
                    open_group = .true.
                    do while(edit_group)
                        ! Refresh display
                        call execute_command_line("clear")
                        call print_overview(fid)
                        call print_bc_overview(fid,active_topic='States')

                        command = "Enter group name: "
                        call write_line(command, color='blue')
                        read(*,'(A1024)') group_name

                        if (trim(group_name) == '') then
                            edit_group = .false.
                            open_group = .false.
                        else
                            open_group = .true.
                        end if

                        if (open_group) call chidg_edit_boundarycondition_state_group(fid,group_name)
                    end do
                
                !
                ! Remove group
                !
                case(3)
                    command = "Enter group name: "
                    call write_line(command, color='blue')
                    read(*,'(A1024)') group_name

                    call remove_bc_state_group_hdf(fid,trim(group_name))

                case default
                    run_states = .false.
            end select

        end do ! run_states

    end subroutine chidg_edit_boundarycondition_states
    !******************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_boundarycondition_state_group(fid,group_name)
        integer(HID_T),     intent(in)  :: fid
        character(*),       intent(in)  :: group_name

    
        integer(HID_T)              :: bcgroup_id
        integer(ik)                 :: ierr, selection
        character(1024)             :: state_name, bc_string, pname
        character(:),   allocatable :: command
        logical                     :: group_exists, edit_state_group, print_bcs, set_bc, found_bc_state, &
                                       bc_state_exists, remove_bc, get_property, property_exists, edit_property
        class(bc_state_t), allocatable    :: bc_state


        group_exists = check_link_exists_hdf(fid,"BCSG_"//trim(group_name))

        if (group_exists) then


            call h5gopen_f(fid,"BCSG_"//trim(group_name),bcgroup_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_boundarycondition_state_group: h5gopen_f")

            edit_state_group = .true.
            do while(edit_state_group)

                !
                ! Refresh display
                !
                call execute_command_line("clear")
                call print_overview(fid)
                call print_bc_overview(fid,active_topic='States',active_group=trim(group_name))

                command = "1: Add state function, 2: Remove state function, 3: Edit Property, (0 to exit):"
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
                                call list_bcs()
                            end if

                            if (.not. found_bc_state) then
                                call write_line("The boundary condition state specified wasn't found in the ChiDG library. &
                                                 Maybe check that the string was input correctly?",color='red')
                            end if

                            ! Get bc_state string from user
                            command = "Enter boundary condition state(? to list): "
                            call write_line(' ')
                            call write_line(command,color='blue')
                            read(*,"(A1024)") bc_string

                            ! Check for user exit
                            set_bc = (trim(bc_string) /= "")


                            if (set_bc) then
                                if ( trim(bc_string) == '?' ) then
                                    print_bcs = .true.
                                else
                                
                                    ! Call routine to set boundary condition in hdf file.
                                    found_bc_state = check_bc_state_registered(bc_string)
                                    if (found_bc_state) then
                                        call create_bc(bc_string,bc_state)
                                        call add_bc_state_hdf(bcgroup_id,bc_state)
                                        set_bc = .false.
                                    end if

                                end if
                            end if


                        end do ! set_bc


                    !
                    ! Remove bc_state
                    !
                    case(2)
                        bc_state_exists = .true.
                        remove_bc = .true.
                        do while (remove_bc)

                            if (.not. bc_state_exists) then
                                call write_line("The boundary condition state specified for removal wasn't found on the face")
                            end if


                            ! Get boundary condition state from user to remove
                            command = "Enter boundary condition state to remove: "
                            call write_line(' ')
                            call write_line(command,color='blue')
                            read(*,"(A1024)") bc_string

                            ! Call routine to remove boundary condition state
                            bc_state_exists = check_bc_state_exists_hdf(bcgroup_id,trim(bc_string))

                            if (bc_state_exists) then
                                call remove_bc_state_hdf(bcgroup_id,trim(bc_string))
                                remove_bc = .false.
                            end if

                            ! Exit if blank
                            if (trim(bc_string) == "") remove_bc = .false.

                        end do !remove_bc


                    !
                    ! Edit property
                    !
                    case(3)

                        edit_property = .true.
                        do while (edit_property)
                            get_property = .true.
                            do while ( get_property )

                                ! Refresh display
                                call execute_command_line("clear")
                                call print_overview(fid)
                                call print_bc_overview(fid,active_topic='States',active_group=trim(group_name))
                                call print_bc_state_properties(bcgroup_id)

                                ! Get property
                                call write_line(" ")
                                call write_line("Enter boundary condition property: ",color='blue')
                                read(*,"(A1024)") pname

                                ! Check for user exit
                                if (trim(pname) == "") then
                                    get_property = .false.
                                    edit_property = .false.
                                end if

                                ! Check property exists
                                property_exists = check_bc_property_exists_hdf(bcgroup_id,trim(pname))
                                if ( property_exists ) get_property = .false.


                            end do 

                            call chidg_edit_boundarycondition_property(fid,bcgroup_id,trim(group_name),trim(pname))
                        end do



                    case default
                        edit_state_group = .false.
                end select




            end do !edit_state_group


            ! Close group
            call h5gclose_f(bcgroup_id,ierr)

        end if !group_exists

    end subroutine chidg_edit_boundarycondition_state_group
    !******************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_boundarycondition_domains(fid)
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
            call execute_command_line("clear")
            call print_overview(fid)
            call print_bc_overview(fid,active_topic='Patches')


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



            if (open_domain) call chidg_edit_boundarycondition_domain_patches(fid,trim(domain_name))

!            ndom = get_ndomains_hdf(fid)
!            ierr = 1
!            do while ( ierr /= 0 )
!                read(*,'(I8)', iostat=ierr) idom_hdf
!                if ( ierr /= 0 )  call write_line("Invalid input: expecting an integer index.")
!
!                if ( idom_hdf > ndom ) then
!                    ierr = 1
!                    call write_line("Invalid domain range. Enter a number between 1 and ",ndom)
!                end if
!
!            end do
!
!
!            if (idom_hdf /= 0) then
!                call chidg_edit_boundarycondition_domain_patches(fid,idom_hdf)
!            else
!                run_domains = .false.
!            end if

        end do

    end subroutine chidg_edit_boundarycondition_domains
    !******************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    !subroutine chidg_edit_boundarycondition_domain_patches(fid,idom_hdf)
    subroutine chidg_edit_boundarycondition_domain_patches(fid,domain_name)
        integer(HID_T),     intent(in)  :: fid
        character(*),       intent(in)  :: domain_name

        integer(HID_T)              :: bcgroup, domain_id
        integer(ik)                 :: ierr, iface
        character(:),   allocatable :: dname, dname_trim, command
        logical                     :: run_edit_bc_domain

        ! Open domain
        domain_id = open_domain_hdf(fid,domain_name)


        !
        ! Open boundary condition group
        !
        call h5gopen_f(domain_id, "BoundaryConditions", bcgroup, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_boundarycondition_domain: h5gopen - BoundaryConditions")


        !
        ! edit / boundary condition / domain
        !
        run_edit_bc_domain = .true.
        do while ( run_edit_bc_domain )

            !
            ! Refresh display
            !
            call execute_command_line("clear")
            call print_overview(fid,domain_id)
            call print_bc_overview(fid,active_topic='Patches',active_domain=domain_id)
            dname_trim = trim(adjustl(domain_name)) 


            !
            ! Print command options, accept user selection.
            !
            command = "Select a patch for editing(1-6, 0 to exit):"
            call write_line(' ')
            call write_line(command,color='blue')
            ierr = 1
            do while ( ierr /= 0 )
                read(*,'(I8)', iostat=ierr) iface
                if ( ierr /= 0 )  call write_line("Invalid input: expecting an integer index.")

                if ( iface > NFACES ) then
                    ierr = 1
                    call write_line("Invalid domain range. Enter a number between 1 and ",NFACES)
                end if

            end do



            !
            ! Operate on particular domain
            !
            select case (iface)
                case (0)
                    run_edit_bc_domain = .false.

                case default
                    !call chidg_edit_boundarycondition_patch(fid,bcgroup,idom_hdf,iface)
                    call chidg_edit_boundarycondition_patch(fid,bcgroup,domain_id,iface)

            end select


        end do  ! run_bc_edit


        !
        ! Close groups
        !
        call h5gclose_f(bcgroup,ierr)

        call close_domain_hdf(domain_id)

    end subroutine chidg_edit_boundarycondition_domain_patches
    !******************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_edit_boundarycondition_patch(fid,bcgroup,domain_id,iface)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: bcgroup
        integer(HID_T),     intent(in)  :: domain_id
        integer(ik),        intent(in)  :: iface


        integer(HID_T)              :: patch_id
        character(10)               :: patches(NFACES)
        integer(ik)                 :: ierr, selection
        character(1024)             :: dname, group
        character(:),   allocatable :: command
        logical                     :: run, run_set, set_group


        
        !
        ! Open boundary condition face group
        !
        patches = ["  XI_MIN","  XI_MAX"," ETA_MIN"," ETA_MAX","ZETA_MIN","ZETA_MAX"]
        call h5gopen_f(bcgroup, trim(adjustl(patches(iface))), patch_id, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_boundarycondition_patch: error opening group for boundary condition face")


        !
        ! edit / boundary condition / domain / face
        !
        run = .true.
        do while ( run )

            ! Refresh display
            call execute_command_line("clear")
            call print_overview(fid,domain_id)
            call print_bc_overview(fid,active_topic='Patches',active_domain=domain_id,active_face=iface)


            ! Print command options, accept user selection.
            command = "1:Set boundary state group, 2:Clear boundary state group, 0:Exit"
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
                    
                    call write_line('Enter the Boundary State Group to set: ')
                    read(*,'(A1024)') group
                    
                    ! Check for user exit
                    set_group = (trim(group) /= "")

                    if (set_group) then
                        call set_bc_patch_group_hdf(patch_id,group)
                    end if

                case(2)

                    call set_bc_patch_group_hdf(patch_id,'empty')

                case default
                    run = .false.
            end select


        end do  ! run




        !
        ! Close boundary condition face group
        !
        call h5gclose_f(patch_id, ierr)

    end subroutine chidg_edit_boundarycondition_patch
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
    subroutine chidg_edit_boundarycondition_property(fid,bcgroup_id,group_name,pname)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: bcgroup_id
        character(*),       intent(in)  :: group_name
        character(*),       intent(in)  :: pname

        integer(HID_T)                  :: bcprop, bc_state
        integer(HSIZE_T)                :: adim
        character(len=:),   allocatable :: command
        character(len=1024)             :: option, function_string, gname
        logical                         :: option_exists, run, property_found, property_exists, set_fcn, list_fcns
        real(rk)                        :: val
        integer                         :: ierr, type, igrp, iop, nmembers
        type(svector_t)                 :: bc_state_strings
        type(string_t)                  :: string
        class(function_t),  allocatable :: func


        ! Refresh display
        call execute_command_line("clear")
        call print_overview(fid)
        call print_bc_overview(fid,active_topic='States',active_group=trim(group_name))
        call print_bc_state_properties(bcgroup_id,trim(pname))

        !
        ! Loop through the states to find the property name
        !
        call h5gn_members_f(bcgroup_id, ".", nmembers, ierr)


        !
        !  Loop through groups and delete properties
        !
        if ( nmembers > 0 ) then

            ! First get number of states. This could be different than number of groups.
            do igrp = 0,nmembers-1

                ! Get group name
                call h5gget_obj_info_idx_f(bcgroup_id, ".", igrp, gname, type, ierr)

                ! Test if group is a boundary condition function. 'BCS_'
                if (gname(1:4) == 'BCS_') then
                    call bc_state_strings%push_back(string_t(trim(gname)))
                end if

            end do  ! igrp

        end if



        !
        ! Find the state with the property
        !
        property_found = .false.
        do iop = 1,bc_state_strings%size()

            ! Open the state group
            string = bc_state_strings%at(iop)
            call h5gopen_f(bcgroup_id, string%get(), bc_state, ierr)

            ! Check if it contains a link to the property group
            call h5lexists_f(bc_state, "BCP_"//trim(pname), property_exists, ierr)

            if (property_exists) then
                call h5gopen_f(bc_state, "BCP_"//trim(pname), bcprop, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_boundarycondition_property: error opening bcproperty")
                property_found = .true.
                exit
            end if

            call h5gclose_f(bc_state,ierr)
        end do !iop



        !
        ! Edit option loop
        !
        list_fcns = .false.
        run = property_found
        do while (run)

            ! Refresh display
            call execute_command_line("clear")
            call print_overview(fid)
            call print_bc_overview(fid,active_topic='States',active_group=trim(group_name))
            call print_bc_state_properties(bcgroup_id,trim(pname))

            !
            ! Read option to edit
            !
            command = "Enter 'Function' or option name: "
            call write_line(" ")
            call write_line(command, color='blue')
            read(*,"(A1024)") option
            run = (trim(option) /= "")



            if (run) then

                ! Check option exists
                call h5aexists_f(bcprop, trim(option), option_exists, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_boundarycondition_property: error opening option attribute")


                !
                ! Handle option_exists
                !
                if ( option_exists ) then


                    if ( trim(option) == "Function" ) then

                        !
                        ! Set function
                        !
                        command = "Set function (? to list): "
                        call write_line(" ")
                        call write_line(command, color='blue')
                        read(*,"(A1024)") function_string




                        set_fcn = (trim(function_string) /= "")
                        list_fcns = (trim(function_string) == "?")

                        if (list_fcns) then
                            call list_functions()
                        end if

                        if (set_fcn .and. (.not. list_fcns)) then
                            call create_function(func,trim(function_string))
                            call set_bc_property_function_hdf(bcprop,func)
                            !run = .false.
                        end if


                    else

                        command = "Set option value: "
                        call write_line(command, color='blue')
                        read(*,*) val

                        !
                        ! Set option
                        !
                        adim = 1
                        call h5ltset_attribute_double_f(bc_state, "BCP_"//trim(pname), trim(option), [real(val,rdouble)], adim, ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_boundarycondition_property: error setting option value")

                        run = .false.

                    end if


                else

                    call write_line("Invalid option", color='blue')

                end if

            end if

        end do ! run


        ! Close bcproperty
        if (property_found) then
            call h5gclose_f(bcprop,ierr)
        end if

    end subroutine chidg_edit_boundarycondition_property
    !******************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_bc_overview(fid,active_topic,active_group,active_domain,active_face)
        integer(HID_T),     intent(in)              :: fid
        character(*),       intent(in), optional    :: active_topic
        character(*),       intent(in), optional    :: active_group
        integer(HID_T),     intent(in), optional    :: active_domain
        integer(ik),        intent(in), optional    :: active_face


        call write_line(' ')
        call write_line(' ')
        call print_bc_states(fid,active_topic,active_group)

        call write_line(' ')
        call write_line(' ')
        call write_line(' ')
        call write_line(' ')

        call print_bc_patches(fid,active_topic,active_domain,active_face)
        call write_line(' ')
        call write_line(' ')

    end subroutine print_bc_overview
    !******************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_bc_states(fid,active_topic,active_group)
        integer(HID_T), intent(in)              :: fid
        character(*),   intent(in), optional    :: active_topic
        character(*),   intent(in), optional    :: active_group

        integer(HID_T)              :: bcgroup_id
        integer(ik)                 :: igroup, istate, ierr
        type(svector_t)             :: bc_state_groups, bc_state_names
        type(string_t)              :: group_name, bc_state_name
        character(:),   allocatable :: group_family, color
        logical                     :: bold


        color = 'none'
        if (present(active_topic)) then
            if (trim(active_topic) == 'States') then
                color = 'blue'
            end if
        end if


        !
        ! Write boundary state overview header
        !
        if (present(active_topic)) then 
            if (trim(active_topic) == 'States') then
                call write_line(":[Boundary State Groups]",color=color,bold=.true.)
            else
                call write_line(":Boundary State Groups")
            end if
        else
            call write_line(":Boundary State Groups")
        end if



        ! Write state group header
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
        call write_line("Group","Family","States", columns=.true.,column_width=25,delimiter=':',color=color)
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
            
        !
        ! Get boundary condition information
        !
        bc_state_groups = get_bc_state_group_names_hdf(fid)


        do igroup = 1,bc_state_groups%size()

            group_name = bc_state_groups%at(igroup)
            bcgroup_id = open_bc_group_hdf(fid,group_name%get())

            group_family = get_bc_state_group_family_hdf(bcgroup_id)

            if (group_name%get() == trim(active_group)) then
                color = 'blue'
                bold  = .true.
            else
                color = 'none'
                bold = .false.
            end if


            call add_to_line(trim(group_name%get()), columns=.true., column_width=25,color=color,bold=bold)
            call add_to_line(trim(group_family), columns=.true., column_width=25,color=color,bold=bold)



            bc_state_names = get_bc_state_names_hdf(bcgroup_id)
            do istate = 1,bc_state_names%size()
                bc_state_name = bc_state_names%at(istate)
                call add_to_line(trim(bc_state_name%get()), columns=.true., column_width=25,color=color,bold=bold)
            end do

            call close_bc_group_hdf(bcgroup_id)

            call send_line()

        end do !istate


    end subroutine print_bc_states
    !******************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_bc_patches(fid,active_topic,active_domain,active_face)
        integer(HID_T),     intent(in)              :: fid
        character(*),       intent(in), optional    :: active_topic
        integer(HID_T),     intent(in), optional    :: active_domain
        integer(ik),        intent(in), optional    :: active_face


        character(len=10)                   :: faces(NFACES)
        integer(ik)                         :: idom_hdf, iface, ierr
        integer(HID_T)                      :: patch_id
        type(svector_t),    allocatable     :: bcs(:)
        character(:),       allocatable     :: bc_patch_group
        type(string_t)                      :: bc_patch_groups(NFACES)
        character(len=1024), allocatable    :: dnames(:)
        character(:),       allocatable     :: color, domain_name, active_domain_name


        color = 'none'
        if (present(active_topic)) then
            if (trim(active_topic) == 'Patches') color = 'blue'
        end if

        !
        ! Write boundary patch overview header
        !
        if (present(active_topic)) then 
            if (trim(active_topic) == 'Patches') then
                call write_line(":[Boundary Patches]",color=color,bold=.true.)
            else
                call write_line(":Boundary Patches")
            end if
        else
            call write_line(":Boundary Patches")
        end if


        !
        ! Get boundary condition information
        !
        dnames = get_domain_names_hdf(fid)


        !
        ! Write domain and boundary conditions. 
        ! TODO: Assumes NFACES=6. Could generalize.
        !
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
        call write_line("Domain name",  &
                        "1 - XI_MIN",   &
                        "2 - XI_MAX",   &
                        "3 - ETA_MIN",  &
                        "4 - ETA_MAX",  &
                        "5 - ZETA_MIN", &
                        "6 - ZETA_MAX", columns=.True., column_width=16, delimiter=':',color=color)
        call write_line("----------------------------------------------------------------------------------------------------------------------",color=color)
        do idom_hdf = 1,size(dnames)

            !
            ! Get domain name and bc patch groups associated with idom_hdf
            !
            domain_name = dnames(idom_hdf)
            

            faces = ["  XI_MIN","  XI_MAX"," ETA_MIN"," ETA_MAX","ZETA_MIN","ZETA_MAX"]
            do iface = 1,NFACES

                call h5gopen_f(fid, "/D_"//trim(domain_name)//"/BoundaryConditions/"//trim(adjustl(faces(iface))), patch_id, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"print_bc_patches: h5gopen_f - patch (ex. XI_MIN)")
                
                bc_patch_group = get_bc_patch_group_hdf(patch_id)
                call bc_patch_groups(iface)%set(bc_patch_group)

                call h5gclose_f(patch_id,ierr)

            end do !iface



            if (present(active_domain)) then

                active_domain_name = get_domain_name_hdf(active_domain)
            
                !
                ! Print active domain
                !
                !if (idom_hdf == active_domain) then
                if (trim(active_domain_name) == trim(domain_name)) then


                    ! Need to add information individually to selectively color the entries.
                    call add_to_line(domain_name, columns=.True., column_width=15, color='pink')
                    do iface = 1,6

                        if ( present(active_face) ) then
                            if ( iface == active_face ) then
                                call add_to_line( "["//bc_patch_groups(iface)%get()//"]", columns=.True., column_width=15, color='blue',bold=.true.)
                            else
                                call add_to_line( bc_patch_groups(iface)%get(), columns=.True., column_width=15, color='pink')
                            end if
                        else
                            call add_to_line( bc_patch_groups(iface)%get(), columns=.True., column_width=15, color='pink')
                        end if

                    end do
                    call send_line()



                !
                ! Print non-active domains
                !
                else

                    call add_to_line( domain_name, columns=.True., column_width=15)
                    do iface = 1,6
                        call add_to_line( bc_patch_groups(iface)%get(), columns=.True., column_width=15)
                    end do
                    call send_line()

                end if

            else
                !
                ! Print domain info if no active domain is present
                !
                call add_to_line( domain_name, columns=.True., column_width=15)
                do iface = 1,6
                    call add_to_line( bc_patch_groups(iface)%get(), columns=.True., column_width=15)
                end do
                call send_line()

            end if

        end do



    end subroutine print_bc_patches
    !******************************************************************************************












    !>  Print the properties associated with a given boundary condition face group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!  @note   Modified to include bc_states
    !!
    !!  @param[in]  bcface  HDF5 group identifier. Should be associated with a boundary condition.
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_bc_state_properties(bcgroup_id,active_property)
        integer(HID_T),     intent(in)              :: bcgroup_id
        character(*),       intent(in), optional    :: active_property

        integer(HID_T)          :: bcprop, bc_state
        integer                 :: nmembers, igrp, type, ierr, iop
        character(len=1024)     :: gname, bcname
        type(svector_t)         :: bc_state_strings
        type(string_t)          :: bcstring


        !
        ! Get state groups
        !
        call h5gn_members_f(bcgroup_id, ".", nmembers, ierr)

        if (nmembers > 0) then
            do igrp = 0,nmembers-1
                
                    ! Get group name
                    call h5gget_obj_info_idx_f(bcgroup_id, ".", igrp, gname, type, ierr)

                    if (gname(1:4) == 'BCS_') then
                        call bc_state_strings%push_back(string_t(trim(gname(5:))))
                    end if
            end do
        end if

        call write_line("======================================================================================================================")

        !
        ! Loop through and print states + properties
        !
        do iop = 1,bc_state_strings%size()

            ! Write state name
            bcstring = bc_state_strings%at(iop)
            call write_line(" ")
            call write_line('['//trim(adjustl(bcstring%get()))//']', color='blue')

            call write_line("----------------------------------------------------------------------------------------------------------------------")
            call write_line("Property","Function","Function Options","Values", columns=.true.,column_width=25,delimiter=':')
            call write_line("----------------------------------------------------------------------------------------------------------------------")
            

            ! Open state group
            call h5gopen_f(bcgroup_id, "BCS_"//bcstring%get(), bc_state, ierr)
            


            !
            !  Get number of groups linked to the current bcgroup_id
            !
            call h5gn_members_f(bc_state, ".", nmembers, ierr)


            !
            !  Loop through groups and print property options
            !
            if ( nmembers > 0 ) then
                do igrp = 0,nmembers-1


                    ! Get group name
                    call h5gget_obj_info_idx_f(bc_state, ".", igrp, gname, type, ierr)

                    ! Test if group is a boundary condition function. 'BCP_'
                    if (gname(1:4) == 'BCP_') then

                        ! Open bcproperty_t group
                        call h5gopen_f(bc_state, gname, bcprop, ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,"print_bc_state_properties: error opening bcproperty")


                        ! Print property options
                        call print_bc_state_property_options(bcprop,gname(5:),active_property)


                        ! Close bcproperty
                        call h5gclose_f(bcprop,ierr)

                    end if


                end do  ! igrp
            end if ! nmembers

            ! Close state group
            call h5gclose_f(bc_state,ierr)

        end do !iop

        call write_line("======================================================================================================================")

    end subroutine print_bc_state_properties
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
    subroutine print_bc_state_property_options(bcprop,property_name,active_property)
        use iso_fortran_env
        integer(HID_T),     intent(in)              :: bcprop
        character(*),       intent(in)              :: property_name
        character(*),       intent(in), optional    :: active_property


        integer(ik)                             :: nattr, ierr, lines_printed
        integer(HSIZE_T)                        :: iattr, idx
        character(len=1024),    allocatable     :: option_keys(:)
        character(len=1024)                     :: fcn
        character(:),           allocatable     :: color
        real(rk),               allocatable     :: option_vals(:)
        type(h5o_info_t),       target          :: h5_info
        real(rdouble),   dimension(1)                :: buf



        !
        ! Get property function
        !
        call h5ltget_attribute_string_f(bcprop, ".", 'Function', fcn, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"print_property_options: error retrieving property function name")


        !
        ! Get number of attributes attached to the group id
        !
        call h5oget_info_f(bcprop, h5_info, ierr)
        nattr = h5_info%num_attrs
        if (ierr /= 0) call chidg_signal(FATAL,"print_property_options: error getting current number of attributes.")


        allocate(option_keys(nattr), option_vals(nattr), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Gather any existing attributes and their values
        !
        if ( nattr > 0 ) then
            idx = 0
            do iattr = 1,nattr
                idx = iattr - 1
                call h5aget_name_by_idx_f(bcprop, ".", H5_INDEX_CRT_ORDER_F, H5_ITER_NATIVE_F, idx, option_keys(iattr), ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"print_property_options: error reading attribute name")

                if ( trim(option_keys(iattr)) /= 'Function' ) then  ! don't read function. not a floating point attribute.
                    call h5ltget_attribute_double_f(bcprop, ".", option_keys(iattr), buf, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"print_property_options: error reading attribute value")
                    option_vals(iattr) = real(buf(1),rk)
                end if
            end do

        end if ! nattr


        if (present(active_property)) then
            if (trim(active_property) == trim(property_name)) then
                color = 'blue'
            else
                color = 'none'
            end if
        else
            color = 'none'
        end if


        !
        ! Print the gathered property options. Except the function attribute, which exists for all
        ! properties and was printed above.
        !
        lines_printed = 0
        do iattr = 1,nattr

            if ( (iattr == 1) .and. (trim(option_keys(iattr)) /= "Function") ) then
                call write_line(trim(adjustl(property_name)), trim(adjustl(fcn)), option_keys(iattr), option_vals(iattr), columns=.true., column_width = 25,color=color)
                lines_printed = lines_printed + 1

            else if ( (iattr == 1) .and. (trim(option_keys(iattr)) == "Function") ) then
                ! Don't print this line. We don't want to print "Function"
            else
                if (lines_printed == 0) then
                    call write_line(trim(adjustl(property_name)), trim(adjustl(fcn)), option_keys(iattr), option_vals(iattr), columns=.true., column_width = 25,color=color)
                else
                    call write_line(' ', ' ', option_keys(iattr), option_vals(iattr), columns=.true., column_width = 25,color=color)
                end if
                lines_printed = lines_printed + 1
            end if

        end do




    end subroutine print_bc_state_property_options
    !******************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/9/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_boundary_state_group_help()
        call write_line("Boundary State Groups:", color='blue', width=100)
        call write_line("These are formed to impose a condition or set of &
                         conditions on a region. For example, one might do this for the RANS &
                         equations; creating a group called 'Inlet' and adding a Total Inlet &
                         boundary state in addition to an &
                         inlet state for the turbulence working variables. In order for these &
                         state functions to be applied on a boundary, they must be associated &
                         with a boundary condition patch, or set of boundary condition patches. &
                         The association can be made by editing a patch and choosing the association.", color='red', width=100)

    end subroutine print_boundary_state_group_help
    !******************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/9/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_boundary_patch_help()
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

    end subroutine print_boundary_patch_help
    !******************************************************************************************









end module mod_chidg_edit_boundaryconditions
