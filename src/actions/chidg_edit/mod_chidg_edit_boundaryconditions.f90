module mod_chidg_edit_boundaryconditions
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: NFACES, XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use type_bc,        only: bc_t
    use mod_bc,         only: create_bc
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
            call print_overview_boundaryconditions(fid)



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
    subroutine print_overview_boundaryconditions(fid,active_domain)
        integer(HID_T),     intent(in)              :: fid
        integer(ik),        intent(in), optional    :: active_domain


        integer(ik)                         :: idom_hdf, ndom
        character(len=1024),    allocatable :: bcs(:)
        character(len=1024)                 :: dname


        !
        ! Write boundary condition overview header
        !
        call write_line(" ")
        call write_line("---------------------------------------------------------------------------------")
        call write_line(":Boundary Conditions")
        call write_line("---------------------------------------------------------------------------------")




        !
        ! Get boundary condition information
        !
        ndom   = get_ndomains_hdf(fid)




        !
        ! Write domain and boundary conditions. 
        ! TODO: Assumes NFACES=6. Could generalize.
        !
        call write_line("Domain name","1 - XI_MIN","2 - XI_MAX","3 - ETA_MIN","4 - ETA_MAX","5 - ZETA_MIN","6 - ZETA_MAX", columns=.True., column_width=20)
        do idom_hdf = 1,ndom

            !
            ! Get domain name and bcs associated with idom_hdf
            !
            dname = get_domain_name_hdf(fid,idom_hdf)
            bcs   = get_bcnames_hdf(fid,dname)
            

            call write_line(' ')


            if (present(active_domain)) then
            
                if (idom_hdf == active_domain) then
                    call write_line(dname(3:),  bcs(XI_MIN),   bcs(XI_MAX),  &
                                                bcs(ETA_MIN),  bcs(ETA_MAX), & 
                                                bcs(ZETA_MIN), bcs(ZETA_MAX), columns=.True., column_width=20, color='pink')
                else
                    call write_line(dname(3:),  bcs(XI_MIN),   bcs(XI_MAX),  &
                                                bcs(ETA_MIN),  bcs(ETA_MAX), & 
                                                bcs(ZETA_MIN), bcs(ZETA_MAX), columns=.True., column_width=20)
                end if

            else
                call write_line(dname(3:),  bcs(XI_MIN),   bcs(XI_MAX),  &
                                            bcs(ETA_MIN),  bcs(ETA_MAX), & 
                                            bcs(ZETA_MIN), bcs(ZETA_MAX), columns=.True., column_width=20)
            end if

        end do




    end subroutine print_overview_boundaryconditions
    !************************************************************************************************














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
        if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5gopen - BoundaryConditions")




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
            call print_overview_boundaryconditions(fid,idom_hdf)


            !
            ! Write boundary condition domain overview header
            !
            dname_trim = trim(adjustl(dname)) 
            call write_line(" ")
            call write_line("---------------------------------------------------------------------------------")
            call write_line("::Domain",dname_trim(3:))
            call write_line("---------------------------------------------------------------------------------")
            


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
        character(len=1024),    allocatable :: bcnames(:)
        character(len=1024)                 :: dname
        character(len=:),       allocatable :: command, dname_trim
        character(len=1024)                 :: bc_string, pname
        logical                             :: run_edit_bc_face


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
            call print_overview_boundaryconditions(fid,idom_hdf)


            !
            ! Write boundary condition domain overview header
            !
            dname_trim = trim(adjustl(dname)) 
            call write_line(" ")
            call write_line("---------------------------------------------------------------------------------")
            call write_line("::Domain",dname_trim(3:),"     :::Face", iface)
            call write_line("---------------------------------------------------------------------------------")
            



            !
            ! Check current face boundary condition
            !
            if ( bcnames(iface) == 'empty' ) then

            else
                !call print_boundarycondition_options(bcface)
            end if




            !
            ! Print command options, accept user selection.
            !
            command = "1:Set boundary condition, 2: Select property, 0:Exit"
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


                case (1)
                    !
                    ! Get boundary condition string from user
                    !
                    command = "Enter boundary condition(0 to exit): "
                    call write_line(' ')
                    call write_line(command,color='blue')
                    read(*,*) bc_string

                
                    !
                    ! Call routine to set boundary condition in hdf file.
                    !
                    call set_boundarycondition_face(bcgroup,bcface,bc_string)



                case (2)
                    !
                    ! Call edit option
                    !
                    command = "Enter boundary condition property: "
                    call write_line(command,color='blue')
                    read(*,*) pname

                    !call select_boundarycondition_property(bcface,pname)



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
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine set_boundarycondition_face(bcgroup,bcface,bcstring)
        integer(HID_T),     intent(in)  :: bcgroup
        integer(HID_T),     intent(in)  :: bcface
        character(*),       intent(in)  :: bcstring


        !
        ! Delete current boundary condition and settings
        !
        print*, 'deleting current boundary condition information'
        call delete_boundarycondition_hdf(bcgroup,bcface)



        !
        ! Add new boundary condition and settings
        !
        print*, 'adding new boundary condition'
        call add_boundarycondition_hdf(bcgroup,bcface,bcstring)


    end subroutine set_boundarycondition_face
    !***************************************************************************************************












    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine add_boundarycondition_hdf(bcgroup,bcface,bcstring)
        integer(HID_T),     intent(in)      :: bcgroup
        integer(HID_T),     intent(in)      :: bcface
        character(*),       intent(in)      :: bcstring

        class(bc_t),    allocatable :: bc
        integer(ik)                 :: nfcns, ifcn, ierr

        print*, 'adding bc_name'
        !
        ! Set boundary condition name attribute: bc_name
        !
        call h5ltset_attribute_string_f(bcface, ".", "bc_name", bcstring, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5ltset_attribute_string")


        print*, 'creating boundary condition'
        !
        ! Create an instance of the specified boundary condition to query its options
        !
        call create_bc(bcstring,bc)
    

        print*, 'adding bcproperties'
        !
        ! Add boundary condition functions
        !
        call add_bcproperties_hdf(bcface,bc)


    end subroutine add_boundarycondition_hdf
    !***************************************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine add_bcproperties_hdf(bcface,bc)
        integer(HID_T),     intent(in)      :: bcface
        class(bc_t),        intent(inout)      :: bc

        integer(HID_T)                  :: prop_id
        integer(HSIZE_T)                :: adim
        integer(ik)                     :: iprop, nprop, iopt, nopt, ierr
        character(len=1024)             :: pstring
        character(len=:),   allocatable :: option_key
        real(rk)                        :: option_value


        print*, 'getting number of properties'
        !
        ! Get number of functions in the boundary condition
        !
        nprop = bc%get_nproperties()


        !
        ! Loop through and add properties
        !
        do iprop = 1,nprop

            print*, 'getting property name'
            !
            ! Get string the property is associated with
            !
            pstring = bc%get_property_name(iprop)



            print*, 'creating property group'
            !
            ! Create a new group for the property
            !
            call h5gcreate_f(bcface, "BCP_"//trim(adjustl(pstring)), prop_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"add_bcfunction_hdf: error creating new group for bcfunction")



            print*, 'getting number of options'
            !
            ! Get number of options available for the current property
            !
            nopt = bc%get_noptions(iprop)


            if (nopt > 0 ) then
                do iopt = 1,nopt

                    print*, 'getting option key/value'
                    !
                    ! Get the current option and default value.
                    !
                    option_key   = bc%get_option_key(iprop,iopt)
                    option_value = bc%get_option_value(iprop,option_key)

                    print*, 'setting the option key/value as attribute'
                    !
                    ! Set the option as a real attribute
                    !
                    adim = 1
                    call h5ltset_attribute_double_f(prop_id, ".", option_key, [option_value], adim, ierr)

                end do
            end if


            !
            ! Close function group
            !
            call h5gclose_f(prop_id,ierr)
            


        end do !ifcn

        print*, 'finished adding options'


    end subroutine add_bcproperties_hdf
    !****************************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine delete_boundarycondition_hdf(bcgroup, bcface)
        integer(HID_T),     intent(in)  :: bcgroup
        integer(HID_T),     intent(in)  :: bcface

        integer(HSIZE_T)                        :: iattr
        integer(ik)                             :: nattr
        integer                                 :: nmembers, igrp, type, ierr
        character(len=10)                       :: faces(NFACES)
        character(len=1024),    allocatable     :: anames(:)
        character(len=1024)                     :: gname



        !
        ! Delete overall boundary condition face attributes
        !
        print*, "deleting boundary condition face attributes"
        call delete_group_attributes(bcface)



        !
        !  Get number of groups linked to the current bcface
        !
        call h5gn_members_f(bcface, "/", nmembers, ierr)


        !
        !  Loop through groups and delete properties
        !
        if ( nmembers > 0 ) then

            do igrp = 0,nmembers-1
                !
                ! Get group name
                !
                call h5gget_obj_info_idx_f(bcface,"/", igrp, gname, type, ierr)

                !
                ! Test if group is a boundary condition function. 'BCP_'
                !
                if (gname(1:4) == 'BCP_') then

                    call delete_bcproperty_hdf(bcface,gname)

                end if

            end do  ! igrp

        end if ! nmembers


    end subroutine delete_boundarycondition_hdf
    !****************************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------------
    subroutine delete_bcproperty_hdf(bcface,pname)
        integer(HID_T),     intent(in)      :: bcface
        character(*),       intent(in)      :: pname

        integer(HID_T)  :: bcprop
        integer(ik)     :: ierr

        !
        ! Open bcproperty group
        !
        call h5gopen_f(bcface, trim(adjustl(pname)), bcprop, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"delete_bcproperty_hdf: error opening bcproperty group")


        !
        ! Delete bcproperty attributes
        !
        print*, "deleting boundary condition property attributes"
        call delete_group_attributes(bcprop)


        !
        ! Close bcproperty group
        !
        call h5gclose_f(bcprop, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"delete_bcproperty_hdf: error closing bcproperty group")


        !
        ! Now that the data in bcproperty has been removed, unlink the bcproperty group.
        !
        call h5gunlink_f(bcface,trim(pname),ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"delete_bcfunction_hdf: error unlinking bcproperty group")


    end subroutine delete_bcproperty_hdf
    !******************************************************************************************************

















    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------------
    subroutine print_boundarycondition_options(bcface)
        integer(HID_T),     intent(in)  :: bcface







    end subroutine print_boundarycondition_options















end module mod_chidg_edit_boundaryconditions
