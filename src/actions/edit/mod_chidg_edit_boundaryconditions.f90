module mod_chidg_edit_boundaryconditions
#include <messenger.h>
    use mod_kinds,          only: rk, ik, rdouble
    use mod_constants,      only: NFACES, XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use type_bc_state,   only: bc_state_t
    use mod_bc,             only: create_bc, list_bcs, check_bc_state_registered
    use type_svector,       only: svector_t
    use mod_string,         only: string_t
    use type_function,      only: function_t
    use mod_function,       only: create_function
    use hdf5
    use h5lt

    use mod_hdf_utilities,              only: get_ndomains_hdf, get_domain_names_hdf, get_bcnames_hdf, get_domain_name_hdf, &
                                              delete_group_attributes
    use mod_chidg_edit_printoverview,   only: print_overview
    implicit none



contains





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------
    subroutine chidg_edit_boundaryconditions(fid)
        integer(HID_T),     intent(in)  :: fid


        integer(ik)                     :: ierr, idom_hdf, ndom
        logical                         :: run_bc_edit
        character(len=:),   allocatable :: command



        ndom = get_ndomains_hdf(fid)


        run_bc_edit = .true.
        do while ( run_bc_edit )

            !
            ! Refresh display
            !
            call execute_command_line("clear")
            call print_overview(fid)
            call print_bc_overview(fid)



            !
            ! Print command options, accept user selection.
            !
            command = "Select a domain for editing(0 to exit):"
            call write_line(' ')
            call write_line(command,color='blue')

            ierr = 1
            do while ( ierr /= 0 )
                read(*,'(I8)', iostat=ierr) idom_hdf
                if ( ierr /= 0 )  call write_line("Invalid input: expecting an integer index.")

                if ( idom_hdf > ndom ) then
                    ierr = 1
                    call write_line("Invalid domain range. Enter a number between 1 and ",ndom)
                end if

            end do



            !
            ! Operate on particular domain
            !
            select case (idom_hdf)
                case (0)
                    run_bc_edit = .false.

                case default
                    call chidg_edit_boundarycondition_domain(fid,idom_hdf)

            end select



        end do  ! run_bc_edit



    end subroutine chidg_edit_boundaryconditions
    !************************************************************************************************












    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine print_bc_overview(fid,active_domain,active_face)
        integer(HID_T),     intent(in)              :: fid
        integer(ik),        intent(in), optional    :: active_domain
        integer(ik),        intent(in), optional    :: active_face


        integer(ik)                         :: idom_hdf, ndom, iface
        type(svector_t),    allocatable     :: bcs(:)
        type(string_t)                      :: bcstring
        character(len=1024)                 :: dname


        !
        ! Write boundary condition overview header
        !
        call write_line(" ")
        call write_line(":Boundary Conditions")
        call write_line("______________________________________________________________________________________________________")
        call write_line(" ")
        call write_line(" ")


        !
        ! Get boundary condition information
        !
        ndom   = get_ndomains_hdf(fid)


        !
        ! Write domain and boundary conditions. 
        ! TODO: Assumes NFACES=6. Could generalize.
        !
        call write_line("Domain name","1 - XI_MIN","2 - XI_MAX","3 - ETA_MIN","4 - ETA_MAX","5 - ZETA_MIN","6 - ZETA_MAX", columns=.True., column_width=22)
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
                    call add_to_line(dname(3:), columns=.True., column_width=22, color='pink')
                    do iface = 1,6
                        if ( present(active_face) ) then
                            if ( iface == active_face ) then
                                bcstring = bcs(iface)%at(1)
                                call add_to_line( "["//trim(bcstring%get())//"]", columns=.True., column_width=22, color='blue')
                            else
                                bcstring = bcs(iface)%at(1)
                                call add_to_line( trim(bcstring%get()), columns=.True., column_width=22, color='pink')
                            end if
                        else
                            bcstring = bcs(iface)%at(1)
                            call add_to_line( trim(bcstring%get()), columns=.True., column_width=22, color='pink')
                        end if
                    end do
                    call send_line()



                !
                ! Print non-active domains
                !
                else

                    call add_to_line( dname(3:), columns=.True., column_width=22)
                    do iface = 1,6
                        bcstring = bcs(iface)%at(1)
                        call add_to_line( bcstring%get(), columns=.True., column_width=22)
                    end do
                    call send_line()

                end if

            else
                !
                ! Print domain info if no active domain is present
                !
                call add_to_line( dname(3:), columns=.True., column_width=22)
                do iface = 1,6
                    bcstring = bcs(iface)%at(1)
                    call add_to_line( bcstring%get(), columns=.True., column_width=22)
                end do
                call send_line()

            end if

        end do


    end subroutine print_bc_overview
    !************************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine print_bc_domain(dname)
        character(*),   intent(in)      :: dname

        call write_line(" ")
        call write_line("::Domain",dname)
        call write_line("______________________________________________________________________________________________________")
        call write_line(" ")
        call write_line(" ")

    end subroutine print_bc_domain
    !*************************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine print_bc_domain_face(dname,fname)
        character(*),   intent(in)      :: dname
        character(*),   intent(in)      :: fname


        call write_line(" ")
        call write_line("::Domain",dname,"     :::Face", fname)
        call write_line("______________________________________________________________________________________________________")
        call write_line(" ")
        call write_line(" ")


    end subroutine print_bc_domain_face
    !*************************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine chidg_edit_boundarycondition_domain(fid,idom_hdf)
        integer(HID_T),     intent(in)  :: fid
        integer(ik),        intent(in)  :: idom_hdf

        integer(HID_T)                      :: bcgroup
        integer(ik)                         :: ierr, iface
        character(len=1024)                 :: dname
        character(len=:),       allocatable :: command, dname_trim
        logical                             :: run_edit_bc_domain


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
            call print_bc_overview(fid,idom_hdf)

            dname_trim = trim(adjustl(dname)) 
            call print_bc_domain(dname_trim(3:))


            !
            ! Print command options, accept user selection.
            !
            command = "Select a boundary for editing(0 to exit):"
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
                    call chidg_edit_boundarycondition_face(fid,bcgroup,idom_hdf,iface)

            end select


        end do  ! run_bc_edit


        !
        ! Close BoundaryCondition group
        !
        call h5gclose_f(bcgroup,ierr)


    end subroutine chidg_edit_boundarycondition_domain
    !*************************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine chidg_edit_boundarycondition_face(fid,bcgroup,idom_hdf,iface)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: bcgroup
        integer(ik),        intent(in)  :: idom_hdf
        integer(ik),        intent(in)  :: iface


        integer(HID_T)                      :: bcface
        character(len=10)                   :: faces(NFACES)
        integer(ik)                         :: ierr, int_action
        !character(len=1024),    allocatable :: bcnames(:)
        type(svector_t),        allocatable :: bcnames(:)
        type(string_t)                      :: bcstring
        character(len=1024)                 :: dname
        character(len=:),       allocatable :: command, dname_trim
        character(len=1024)                 :: bc_string, pname
        logical                             :: run_edit_bc_face, get_property, property_exists, set_bc, &
                                               print_bcs, remove_bc, bc_state_exists


        faces = ["  XI_MIN","  XI_MAX"," ETA_MIN"," ETA_MAX","ZETA_MIN","ZETA_MAX"]

        
        dname = get_domain_name_hdf(fid,idom_hdf)


        !
        ! Open boundary condition face group
        !
        call h5gopen_f(bcgroup, trim(adjustl(faces(iface))), bcface, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_boundarycondition_face: error opening group for boundary condition face")


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
            call print_bc_overview(fid,idom_hdf,iface)
            dname_trim = trim(adjustl(dname)) 
            call print_bc_domain_face(dname_trim(3:), trim(adjustl(faces(iface))))


            !
            ! Check current face boundary condition
            !
            bcstring = bcnames(iface)%at(1)
            if ( bcstring%get() == 'empty' ) then

            else
                call print_bc_properties(bcface)
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
                    do while (set_bc)

                        ! Refresh display
                        call execute_command_line("clear")
                        call print_overview(fid,idom_hdf)
                        call print_bc_overview(fid,idom_hdf,iface)
                        dname_trim = trim(adjustl(dname)) 
                        call print_bc_domain_face(dname_trim(3:), trim(adjustl(faces(iface))))

                        if (print_bcs) then
                            call list_bcs()
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
                                call add_bc_state_hdf(bcface,bc_string)
                                set_bc = .false.

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

                        ! Refresh display
                        call execute_command_line("clear")
                        call print_overview(fid,idom_hdf)
                        call print_bc_overview(fid,idom_hdf,iface)
                        dname_trim = trim(adjustl(dname))
                        call print_bc_domain_face(dname_trim(3:), trim(adjustl(faces(iface))))

                        if (.not. bc_state_exists) then
                            call write_line("The boundary condition state specified for removal wasn't found on the face")
                        end if


                        ! Get boundary condition state from user to remove
                        command = "Enter boundary condition state to remove: "
                        call write_line(' ')
                        call write_line(command,color='blue')
                        read(*,"(A1024)") bc_string

                        ! Call routine to remove boundary condition state
                        bc_state_exists = check_bc_state_exists(bcface,trim(bc_string))

                        if (bc_state_exists) then
                            call delete_bc_state_hdf(bcface,trim(bc_string))
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
                        call print_bc_overview(fid,idom_hdf)
                        dname_trim = trim(adjustl(dname)) 
                        call print_bc_domain_face(dname_trim(3:), trim(adjustl(faces(iface))))
                        call print_bc_properties(bcface)

                        ! Call edit option
                        call write_line("Enter boundary condition property: ",color='blue')
                        read(*,"(A1024)") pname

                        ! Check for user exit
                        if (trim(pname) == "") then
                            get_property = .false.
                        end if

                        ! Check property exists
                        property_exists = check_property_exists(bcface,trim(pname))
                        if ( property_exists ) get_property = .false.

                    end do 



                    call edit_property(bcface,trim(pname))



                case default

            end select



        end do  ! run_edit_bc_face




        !
        ! Close boundary condition face group
        !
        call h5gclose_f(bcface, ierr)

    end subroutine chidg_edit_boundarycondition_face
    !**************************************************************************************************











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
    !---------------------------------------------------------------------------------------------------
    subroutine add_bc_state_hdf(bcface,state_string)
        integer(HID_T),     intent(in)      :: bcface
        character(*),       intent(in)      :: state_string

        class(bc_state_t),   allocatable :: bc_state
        integer(ik)                         :: nfcns, ifcn, ierr
        integer(HID_T)                      :: op_id
        logical                             :: link_exists, state_found


        if ( trim(state_string) == 'empty' ) then
            !
            ! If 'empty' do not allocate new bc
            !
            
        else

            ! Check to make sure the bc_state wasn't previously added
            call h5lexists_f(bcface, "BCS_"//trim(adjustl(state_string)), link_exists, ierr)


            if (.not. link_exists) then

                ! Check bc_state exists in the register. 
                ! If not, user probably entered the wrong string, so do nothing
                state_found = check_bc_state_registered(trim(adjustl(state_string)))

                if (state_found) then
                    ! Create a new group for the bc_state_t
                    call h5gcreate_f(bcface, "BCS_"//trim(adjustl(state_string)), op_id, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"add_boundarycondition_state_hdf: error creating new group for bcfunction")

                    ! Create an instance of the specified boundary condition to query its options
                    call create_bc(state_string,bc_state)

                    ! Add bc_state properties to the group that was created
                    call add_bc_properties_hdf(op_id,bc_state)

                    ! Close function group
                    call h5gclose_f(op_id,ierr)
                end if

            end if

        end if



    end subroutine add_bc_state_hdf
    !***************************************************************************************************









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
    !----------------------------------------------------------------------------------------------------
    subroutine add_bc_properties_hdf(op_id,bc_state)
        integer(HID_T),         intent(in)      :: op_id
        class(bc_state_t),   intent(inout)   :: bc_state

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
            call h5gcreate_f(op_id, "BCP_"//trim(adjustl(pstring)), prop_id, ierr)
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
    !----------------------------------------------------------------------------------------------------
    subroutine delete_bc_state_hdf(bcface,state_string)
        integer(HID_T),     intent(in)  :: bcface
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
        call h5gopen_f(bcface, "BCS_"//trim(state_string), bc_state, ierr)


        !
        ! Delete overall boundary condition face attributes
        !
        call delete_group_attributes(bc_state)


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
                call delete_bc_property_hdf(bc_state,pnames(iprop))
            end do


        end if ! nmembers


        !
        ! Close the bc_state group
        !
        call h5gclose_f(bc_state,ierr)

        !
        ! Unlink the bc_state group
        !
        call h5gunlink_f(bcface,"BCS_"//trim(state_string),ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"delete_bc_state_hdf: error unlinking bc_state group")


    end subroutine delete_bc_state_hdf
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
    subroutine delete_bc_property_hdf(bc_state,pname)
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
        call delete_group_attributes(bcprop)


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


    end subroutine delete_bc_property_hdf
    !******************************************************************************************************

















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
    !--------------------------------------------------------------------------------------------------------
    subroutine print_bc_properties(bcface)
        integer(HID_T),     intent(in)  :: bcface

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
            call write_line('['//trim(adjustl(bcstring%get()))//']', color='pink')

            
            ! Open state group
            call h5gopen_f(bcface, "BCS_"//bcstring%get(), bc_state, ierr)
            


            !
            !  Get number of groups linked to the current bcface
            !
            call h5gn_members_f(bc_state, ".", nmembers, ierr)


            !
            !  Loop through groups and delete properties
            !
            if ( nmembers > 0 ) then
                do igrp = 0,nmembers-1


                    !
                    ! Get group name
                    !
                    call h5gget_obj_info_idx_f(bc_state, ".", igrp, gname, type, ierr)

                    !
                    ! Test if group is a boundary condition function. 'BCP_'
                    !
                    if (gname(1:4) == 'BCP_') then
                        !
                        ! Open bcproperty_t group
                        !
                        call h5gopen_f(bc_state, gname, bcprop, ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,"print_bc_properties: error opening bcproperty")


                        !
                        ! Print property name
                        !
                        call write_line(trim(gname(5:)))


                        !
                        ! Print property options
                        !
                        call print_property_options(bcprop)


                        !
                        ! Close bcproperty
                        !
                        call h5gclose_f(bcprop,ierr)
                    end if


                end do  ! igrp
            end if ! nmembers

            ! Close state group
            call h5gclose_f(bc_state,ierr)

        end do !iop


    end subroutine print_bc_properties
    !*****************************************************************************************












    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine print_property_options(bcprop)
        use iso_fortran_env
        integer(HID_T),     intent(in)      :: bcprop


        integer(ik)                             :: nattr, ierr
        integer(HSIZE_T)                        :: iattr, idx
        character(len=1024),    allocatable     :: option_keys(:)
        character(len=1024)                     :: fcn
        real(rk),               allocatable     :: option_vals(:)
        type(h5o_info_t),       target          :: h5_info
        real(rdouble),   dimension(1)                :: buf



        !
        ! Get property function
        !
        call h5ltget_attribute_string_f(bcprop, ".", 'function', fcn, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"print_property_options: error retrieving property function name")
        call write_line(" ")
        call write_line(' f(t,x,y,z) = '//trim(adjustl(fcn)))
        call write_line("____________________________________")
        call write_line(" ")
        



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

                if ( trim(option_keys(iattr)) /= 'function' ) then  ! don't read function. not a floating point attribute.
                    call h5ltget_attribute_double_f(bcprop, ".", option_keys(iattr), buf, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"print_property_options: error reading attribute value")
                    option_vals(iattr) = real(buf(1),rk)
                end if
            end do

        end if ! nattr




        !
        ! Print the gathered property options. Except the function attribute, which exists for all
        ! properties and was printed above.
        !
        do iattr = 1,nattr
            if ( trim(option_keys(iattr)) == 'function' ) then

            else
                call write_line(option_keys(iattr), option_vals(iattr), columns=.true., column_width = 15)
            end if
        end do




    end subroutine print_property_options
    !*******************************************************************************************








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
    subroutine edit_property(bcface,pname)
        integer(HID_T),     intent(in)  :: bcface
        character(*),       intent(in)  :: pname

        integer(HID_T)                  :: bcprop, bc_state
        integer(HSIZE_T)                :: adim
        character(len=:),   allocatable :: command
        character(len=1024)             :: option, new_function, gname
        logical                         :: option_exists, run, property_found, property_exists
        real(rk)                        :: val
        integer                         :: ierr, type, igrp, iop, nmembers
        type(svector_t)                 :: bc_state_strings
        type(string_t)                  :: string


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
                if (ierr /= 0) call chidg_signal(FATAL,"edit_property: error opening bcproperty")
                property_found = .true.
                exit
            end if

            call h5gclose_f(bc_state,ierr)
        end do !iop



        !
        ! Edit option loop
        !
        run = property_found
        do while (run)

            !
            ! Read option to edit
            !
            command = "Enter 'function' or option name: "
            call write_line(command, color='blue')
            read(*,*) option


            !
            ! Check option exists
            !
            call h5aexists_f(bcprop, trim(option), option_exists, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"edit_property: error opening option attribute")


            !
            ! Handle option_exists
            !
            if ( option_exists ) then


                if ( trim(option) == "function" ) then

                    !
                    ! Set function
                    !
                    command = "Set function: "
                    call write_line(command, color='blue')
                    read(*,*) new_function

                    call set_function(bcprop, trim(new_function))


                else

                    command = "Set option value: "
                    call write_line(command, color='blue')
                    read(*,*) val

                    !
                    ! Set option
                    !
                    adim = 1
                    call h5ltset_attribute_double_f(bc_state, "BCP_"//trim(pname), trim(option), [real(val,rdouble)], adim, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"edit_property: error setting option value")

                end if

                ! Close routine
                run = .false.

            else

                call write_line("Invalid option", color='blue')

            end if


        end do ! run


        ! Close bcproperty
        if (property_found) then
            call h5gclose_f(bcprop,ierr)
        end if

    end subroutine edit_property
    !********************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/15/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    function check_bc_state_exists(bcface,bc_state) result(exist_status)
        integer(HID_T),     intent(in)  :: bcface
        character(len=*),   intent(in)  :: bc_state

        integer(ik) :: ierr
        logical     :: exist_status

        ! Check if face contains the bc_state
        call h5lexists_f(bcface, "BCS_"//trim(bc_state), exist_status, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"check_bc_state_exists: Error in call to h5lexists_f")


    end function check_bc_state_exists
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
    function check_property_exists(bcface,pname) result(exist_status)
        integer(HID_T),     intent(in)  :: bcface
        character(*),       intent(in)  :: pname

        integer(HID_T)          :: bc_state
        integer                 :: ierr, nmembers, igrp, type, iop
        character(len=1024)     :: gname
        logical                 :: exist_status
        type(svector_t)         :: bc_state_strings
        type(string_t)          :: string

        
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
        exist_status = .false.
        do iop = 1,bc_state_strings%size()

            ! Open the state group
            string = bc_state_strings%at(iop)
            call h5gopen_f(bcface, string%get(), bc_state, ierr)

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



    end function check_property_exists
    !********************************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine set_function(bcprop, fname)
        integer(HID_T),     intent(in)  :: bcprop
        character(*),       intent(in)  :: fname

        integer(HSIZE_T)                :: adim
        class(function_t),  allocatable :: fcn
        character(len=:),   allocatable :: option
        real(rk)                        :: val
        integer(ik)                     :: nopt, iopt
        integer                         :: ierr
        

        !
        ! Delete bcproperty attributes
        !
        call delete_group_attributes(bcprop)


        !
        ! Set 'function' attribute
        !
        call h5ltset_attribute_string_f(bcprop, ".", "function", trim(fname), ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"set_function: error setting function name")


        !
        ! Create specified function
        !
        call create_function(fcn,fname)


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
            call h5ltset_attribute_double_f(bcprop, ".", trim(option), [real(val,rdouble)], adim, ierr)

        end do ! iopt

    end subroutine set_function
    !********************************************************************************************





end module mod_chidg_edit_boundaryconditions
