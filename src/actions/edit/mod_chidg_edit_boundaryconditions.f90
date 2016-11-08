module mod_chidg_edit_boundaryconditions
#include <messenger.h>
    use mod_kinds,          only: rk, ik, rdouble
    use mod_constants,      only: NFACES, XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use type_bc_state,   only: bc_state_t
    use mod_bc,             only: create_bc, list_bcs, check_bc_state_registered
    use type_svector,       only: svector_t
    use mod_string,         only: string_t
    use type_function,      only: function_t
    use mod_function,       only: create_function, list_functions
    use hdf5
    use h5lt

    use mod_hdf_utilities,              only: get_ndomains_hdf, get_domain_names_hdf, get_bcnames_hdf, &
                                              get_domain_name_hdf, delete_group_attributes_hdf, &
                                              add_bc_state_hdf, set_bc_property_function_hdf, &
                                              get_bc_state_group_names_hdf,                 &
                                              check_bc_property_exists_hdf, remove_bc_state_hdf, &
                                              check_bc_state_exists_hdf
    use mod_chidg_edit_printoverview,   only: print_overview
    implicit none



contains


    !-----------------------------------------------------------------------------------------
    !!
    !!  chidg_edit_boundaryconditions
    !!  chidg_edit_boundarycondition_states
    !!  chidg_edit_boundarycondition_domains
    !!  chidg_edit_boundarycondition_domain_patches
    !!  chidg_edit_boundarycondition_patch
    !!  chidg_edit_boundarycondition_property
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

        integer(ik)                     :: ierr, idom_hdf, ndom, selection
        logical                         :: run_bc_edit, run
        character(len=:),   allocatable :: command

        run_bc_edit = .true.
        do while ( run_bc_edit )

            ! Refresh display
            call execute_command_line("clear")
            call print_overview(fid)
            call print_bc_overview(fid)


            ! Print command options, accept user input.
            command = "1: Edit Boundary State Groups, 2: Edit Boundary Patches, 0: Exit"
            call write_line(' ')
            call write_line(command,color='blue')
            read(*,'(I8)', iostat=ierr) selection

            select case(selection)
                case(1)
                    call chidg_edit_boundarycondition_states(fid)
                case(2)
                    call chidg_edit_boundarycondition_domains(fid)
                case default
                    run_bc_edit = .false.
            end select 


        end do  ! run_bc_edit

    end subroutine chidg_edit_boundaryconditions
    !************************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine chidg_edit_boundarycondition_states(fid)
        integer(HID_T), intent(in)  :: fid

        character(:),   allocatable :: command
        integer(ik)                 :: selection, ierr
        logical                     :: run_states

        run_states = .true.
        do while(run_states)

            !
            ! Refresh display
            !
            call execute_command_line("clear")
            call print_overview(fid)
            call print_bc_overview(fid,active_topic='States')


            command = "1: Create a group for boundary condition states, 2: Select a group, 3: Remove a group, (0 to exit):"
            call write_line(' ')
            call write_line(command,color='blue')

            read(*,'(I8)', iostat=ierr) selection

            select case(selection)
                case(1)

                case(2)

                case(3)

                case default
                    run_states = .false.
            end select

        end do ! run_states

    end subroutine chidg_edit_boundarycondition_states
    !************************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine chidg_edit_boundarycondition_domains(fid)
        integer(HID_T),     intent(in)  :: fid

        integer(ik)                         :: ierr, idom_hdf, ndom
        character(len=:),       allocatable :: command
        logical                             :: run_domains


        run_domains = .true.
        do while(run_domains)

            !
            ! Refresh display
            !
            call execute_command_line("clear")
            call print_overview(fid)
            call print_bc_overview(fid,active_topic='Patches')


            command = "Select a domain for editing(0 to exit):"
            call write_line(' ')
            call write_line(command,color='blue')


            ndom = get_ndomains_hdf(fid)
            ierr = 1
            do while ( ierr /= 0 )
                read(*,'(I8)', iostat=ierr) idom_hdf
                if ( ierr /= 0 )  call write_line("Invalid input: expecting an integer index.")

                if ( idom_hdf > ndom ) then
                    ierr = 1
                    call write_line("Invalid domain range. Enter a number between 1 and ",ndom)
                end if

            end do


            if (idom_hdf /= 0) then
                call chidg_edit_boundarycondition_domain_patches(fid,idom_hdf)
            else
                run_domains = .false.
            end if

        end do

    end subroutine chidg_edit_boundarycondition_domains
    !*************************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine chidg_edit_boundarycondition_domain_patches(fid,idom_hdf)
        integer(HID_T),     intent(in)  :: fid
        integer(ik),        intent(in)  :: idom_hdf

        integer(HID_T)              :: bcgroup
        integer(ik)                 :: ierr, iface
        character(:),   allocatable :: dname, dname_trim, command
        logical                     :: run_edit_bc_domain

        !
        ! Get domain name associated with idom_hdf, the domain index as represented in the hdf file.
        !
        dname = get_domain_name_hdf(fid,idom_hdf)


        !
        ! Open boundary condition group
        !
        call h5gopen_f(fid, "/"//trim(adjustl(dname))//"/BoundaryConditions", bcgroup, ierr)
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
            call print_overview(fid,idom_hdf)
            call print_bc_overview(fid,active_topic='Patches',active_domain=idom_hdf)
            dname_trim = trim(adjustl(dname)) 


            !
            ! Print command options, accept user selection.
            !
            command = "Select a patch for editing(0 to exit):"
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
                    call chidg_edit_boundarycondition_patch(fid,bcgroup,idom_hdf,iface)

            end select


        end do  ! run_bc_edit


        !
        ! Close BoundaryCondition group
        !
        call h5gclose_f(bcgroup,ierr)


    end subroutine chidg_edit_boundarycondition_domain_patches
    !**************************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine chidg_edit_boundarycondition_patch(fid,bcgroup,idom_hdf,iface)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: bcgroup
        integer(ik),        intent(in)  :: idom_hdf
        integer(ik),        intent(in)  :: iface


        integer(HID_T)                      :: bcface
        character(len=10)                   :: faces(NFACES)
        integer(ik)                         :: ierr, int_action
        type(svector_t),        allocatable :: bcnames(:)
        type(string_t)                      :: bcstring
        character(len=1024)                 :: dname
        character(len=:),       allocatable :: command, dname_trim
        character(len=1024)                 :: bc_string, pname
        logical                             :: run_edit_bc_face, get_property, property_exists, set_bc, &
                                               print_bcs, remove_bc, bc_state_exists, found_bc_state
        class(bc_state_t),      allocatable :: bc_state


        faces = ["  XI_MIN","  XI_MAX"," ETA_MIN"," ETA_MAX","ZETA_MIN","ZETA_MAX"]

        
        dname = get_domain_name_hdf(fid,idom_hdf)


        !
        ! Open boundary condition face group
        !
        call h5gopen_f(bcgroup, trim(adjustl(faces(iface))), bcface, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_boundarycondition_patch: error opening group for boundary condition face")


        !
        ! edit / boundary condition / domain / face
        !
        run_edit_bc_face = .true.
        do while ( run_edit_bc_face )

            !
            ! Get updated boundary conditions
            !
            bcnames = get_bcnames_hdf(fid,dname)


            !
            ! Refresh display
            !
            call execute_command_line("clear")
            call print_overview(fid,idom_hdf)
            call print_bc_overview(fid,active_topic='Patches',active_domain=idom_hdf,active_face=iface)
            dname_trim = trim(adjustl(dname)) 


            !
            ! Check current face boundary condition
            !
            bcstring = bcnames(iface)%at(1)
            if ( bcstring%get() == 'empty' ) then

            else
                call print_bc_state_properties(bcface)
            end if



            !
            ! Print command options, accept user selection.
            !
            command = "1:Add boundary state, 2:Remove boundary state, 3: Select property, 0:Exit"
            call write_line(' ')
            call write_line(command,color='blue')
            ierr = 1
            do while ( ierr /= 0 )
                read(*,'(I8)', iostat=ierr) int_action
                if ( ierr /= 0 )  call write_line("Invalid input: expecting an integer index.")
            end do



            !
            ! Operate on user action
            !
            select case (int_action)
                case (0)
                    run_edit_bc_face = .false.


                !
                ! Add bc_state case
                !
                case (1)
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
                                    call add_bc_state_hdf(bcface,bc_state)
                                    set_bc = .false.
                                end if

                            end if
                        end if


                    end do ! set_bc



                !
                ! Remove bc_state case
                !
                case (2)
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
                        bc_state_exists = check_bc_state_exists_hdf(bcface,trim(bc_string))

                        if (bc_state_exists) then
                            call remove_bc_state_hdf(bcface,trim(bc_string))
                            remove_bc = .false.
                        end if

                        ! Exit if blank
                        if (trim(bc_string) == "") remove_bc = .false.

                    end do !remove_bc



                !
                ! Edit property case
                !
                case (3)

                    get_property = .true.
                    do while ( get_property )

                        ! Refresh display
                        call execute_command_line("clear")
                        call print_overview(fid,idom_hdf)
                        call print_bc_overview(fid,active_topic='Patches',active_domain=idom_hdf,active_face=iface)
                        call print_bc_state_properties(bcface)

                        ! Call edit option
                        call write_line(" ")
                        call write_line("Enter boundary condition property: ",color='blue')
                        read(*,"(A1024)") pname

                        ! Check for user exit
                        if (trim(pname) == "") then
                            get_property = .false.
                        end if

                        ! Check property exists
                        property_exists = check_bc_property_exists_hdf(bcface,trim(pname))
                        if ( property_exists ) get_property = .false.

                    end do 

                    call chidg_edit_boundarycondition_property(fid,idom_hdf,iface,bcface,trim(pname))


                case default

            end select



        end do  ! run_edit_bc_face




        !
        ! Close boundary condition face group
        !
        call h5gclose_f(bcface, ierr)

    end subroutine chidg_edit_boundarycondition_patch
    !**************************************************************************************************













    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!  @note   Modified to include bc_states
    !!
    !--------------------------------------------------------------------------------------------
    subroutine chidg_edit_boundarycondition_property(fid,idom_hdf,iface,bcface,pname)
        integer(HID_T),     intent(in)  :: fid
        integer(ik),        intent(in)  :: idom_hdf
        integer(ik),        intent(in)  :: iface
        integer(HID_T),     intent(in)  :: bcface
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
        call print_overview(fid,idom_hdf)
        call print_bc_overview(fid,active_topic='Patches',active_domain=idom_hdf,active_face=iface)
        call print_bc_state_properties(bcface,pname)

        !
        ! Loop through the states to find the property name
        !
        call h5gn_members_f(bcface, ".", nmembers, ierr)


        !
        !  Loop through groups and delete properties
        !
        if ( nmembers > 0 ) then

            ! First get number of states. This could be different than number of groups.
            do igrp = 0,nmembers-1

                ! Get group name
                call h5gget_obj_info_idx_f(bcface, ".", igrp, gname, type, ierr)

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
            call h5gopen_f(bcface, string%get(), bc_state, ierr)

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
                            run = .false.
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
    !********************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine print_bc_overview(fid,active_topic,active_state,active_domain,active_face)
        integer(HID_T),     intent(in)              :: fid
        character(*),       intent(in), optional    :: active_topic
        integer(ik),        intent(in), optional    :: active_state
        integer(ik),        intent(in), optional    :: active_domain
        integer(ik),        intent(in), optional    :: active_face


        call write_line(' ')
        call write_line(' ')
        call print_bc_states(fid,active_topic,active_state)

        call write_line(' ')
        call write_line(' ')
        call write_line(' ')
        call write_line(' ')

        call print_bc_patches(fid,active_topic,active_domain,active_face)
        call write_line(' ')

    end subroutine print_bc_overview
    !************************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------------
    subroutine print_bc_states(fid,active_topic,active_state)
        integer(HID_T), intent(in)              :: fid
        character(*),   intent(in), optional    :: active_topic
        integer(ik),    intent(in), optional    :: active_state

        type(svector_t) :: bc_state_groups

        !
        ! Write boundary state overview header
        !
        if (present(active_topic)) then 
            if (trim(active_topic) == 'States') then
                call write_line(":[Boundary State Groups]",color='blue')
            else
                call write_line(":Boundary State Groups")
            end if
        else
            call write_line(":Boundary State Groups")
        end if



        ! Write state group header
        call write_line("----------------------------------------------------------------------------------------------------------------------")
        call write_line("Group","Family","States", columns=.true.,column_width=25,delimiter=':')
        call write_line("----------------------------------------------------------------------------------------------------------------------")
            
        !
        ! Get boundary condition information
        !
        bc_state_groups = get_bc_state_group_names_hdf(fid)

    end subroutine print_bc_states
    !************************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------------
    subroutine print_bc_patches(fid,active_topic,active_domain,active_face)
        integer(HID_T),     intent(in)              :: fid
        character(*),       intent(in), optional    :: active_topic
        integer(ik),        intent(in), optional    :: active_domain
        integer(ik),        intent(in), optional    :: active_face


        integer(ik)                         :: idom_hdf, ndom, iface
        type(svector_t),    allocatable     :: bcs(:)
        type(string_t)                      :: bcstring
        character(len=1024)                 :: dname


        !
        ! Write boundary patch overview header
        !
        if (present(active_topic)) then 
            if (trim(active_topic) == 'Patches') then
                call write_line(":[Boundary Patches]",color='blue')
            else
                call write_line(":Boundary Patches")
            end if
        else
            call write_line(":Boundary Patches")
        end if


        !
        ! Get boundary condition information
        !
        ndom   = get_ndomains_hdf(fid)


        !
        ! Write domain and boundary conditions. 
        ! TODO: Assumes NFACES=6. Could generalize.
        !
        call write_line("----------------------------------------------------------------------------------------------------------------------")
        call write_line("Domain name",  &
                        "1 - XI_MIN",   &
                        "2 - XI_MAX",   &
                        "3 - ETA_MIN",  &
                        "4 - ETA_MAX",  &
                        "5 - ZETA_MIN", &
                        "6 - ZETA_MAX", columns=.True., column_width=16, delimiter=':')
        call write_line("----------------------------------------------------------------------------------------------------------------------")
        do idom_hdf = 1,ndom

            !
            ! Get domain name and bcs associated with idom_hdf
            !
            dname = get_domain_name_hdf(fid,idom_hdf)
            bcs   = get_bcnames_hdf(fid,dname)
            

            if (present(active_domain)) then
            
                !
                ! Print active domain
                !
                if (idom_hdf == active_domain) then

                    ! Need to add information individually to selectively color the entries.
                    call add_to_line(dname(3:), columns=.True., column_width=15, color='pink')
                    do iface = 1,6
                        if ( present(active_face) ) then
                            if ( iface == active_face ) then
                                bcstring = bcs(iface)%at(1)
                                call add_to_line( "["//trim(bcstring%get())//"]", columns=.True., column_width=15, color='blue')
                            else
                                bcstring = bcs(iface)%at(1)
                                call add_to_line( trim(bcstring%get()), columns=.True., column_width=15, color='pink')
                            end if
                        else
                            bcstring = bcs(iface)%at(1)
                            call add_to_line( trim(bcstring%get()), columns=.True., column_width=15, color='pink')
                        end if
                    end do
                    call send_line()



                !
                ! Print non-active domains
                !
                else

                    call add_to_line( dname(3:), columns=.True., column_width=15)
                    do iface = 1,6
                        bcstring = bcs(iface)%at(1)
                        call add_to_line( bcstring%get(), columns=.True., column_width=15)
                    end do
                    call send_line()

                end if

            else
                !
                ! Print domain info if no active domain is present
                !
                call add_to_line( dname(3:), columns=.True., column_width=15)
                do iface = 1,6
                    bcstring = bcs(iface)%at(1)
                    call add_to_line( bcstring%get(), columns=.True., column_width=15)
                end do
                call send_line()

            end if

        end do



    end subroutine print_bc_patches
    !************************************************************************************************












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
    !-------------------------------------------------------------------------------------------------
    subroutine print_bc_state_properties(bcface,active_property)
        integer(HID_T),     intent(in)              :: bcface
        character(*),       intent(in), optional    :: active_property

        integer(HID_T)          :: bcprop, bc_state
        integer                 :: nmembers, igrp, type, ierr, iop
        character(len=1024)     :: gname, bcname
        type(svector_t)         :: bc_state_strings
        type(string_t)          :: bcstring


        !
        ! Get state groups
        !
        call h5gn_members_f(bcface, ".", nmembers, ierr)

        if (nmembers > 0) then
            do igrp = 0,nmembers-1
                
                    ! Get group name
                    call h5gget_obj_info_idx_f(bcface, ".", igrp, gname, type, ierr)

                    if (gname(1:4) == 'BCS_') then
                        call bc_state_strings%push_back(string_t(trim(gname(5:))))
                    end if
            end do
        end if


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
            call h5gopen_f(bcface, "BCS_"//bcstring%get(), bc_state, ierr)
            


            !
            !  Get number of groups linked to the current bcface
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


    end subroutine print_bc_state_properties
    !**********************************************************************************************












    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------
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
    !***************************************************************************************************






















end module mod_chidg_edit_boundaryconditions
