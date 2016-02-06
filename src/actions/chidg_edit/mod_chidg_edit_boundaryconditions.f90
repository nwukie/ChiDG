module mod_chidg_edit_boundaryconditions
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: NFACES, XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use type_bc,        only: bc_t
    use mod_bc,         only: create_bc
    use type_function,  only: function_t
    use mod_function,   only: create_function
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
            call print_boundaryconditions_overview(fid)



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
    subroutine print_boundaryconditions_overview(fid,active_domain)
        integer(HID_T),     intent(in)              :: fid
        integer(ik),        intent(in), optional    :: active_domain


        integer(ik)                         :: idom_hdf, ndom
        character(len=1024),    allocatable :: bcs(:)
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
        call write_line("Domain name","1 - XI_MIN","2 - XI_MAX","3 - ETA_MIN","4 - ETA_MAX","5 - ZETA_MIN","6 - ZETA_MAX", columns=.True., column_width=20)
        do idom_hdf = 1,ndom

            !
            ! Get domain name and bcs associated with idom_hdf
            !
            dname = get_domain_name_hdf(fid,idom_hdf)
            bcs   = get_bcnames_hdf(fid,dname)
            

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


    end subroutine print_boundaryconditions_overview
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
    subroutine print_boundaryconditions_domain(dname)
        character(*),   intent(in)      :: dname

        call write_line(" ")
        call write_line("::Domain",dname)
        call write_line("______________________________________________________________________________________________________")
        call write_line(" ")
        call write_line(" ")

    end subroutine print_boundaryconditions_domain
    !*************************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine print_boundaryconditions_domain_face(dname,fname)
        character(*),   intent(in)      :: dname
        character(*),   intent(in)      :: fname


        call write_line(" ")
        call write_line("::Domain",dname,"     :::Face", fname)
        call write_line("______________________________________________________________________________________________________")
        call write_line(" ")
        call write_line(" ")


    end subroutine print_boundaryconditions_domain_face
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
            call print_boundaryconditions_overview(fid,idom_hdf)

            dname_trim = trim(adjustl(dname)) 
            call print_boundaryconditions_domain(dname_trim(3:))


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
        logical                             :: run_edit_bc_face, get_property, property_exists


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
            call print_boundaryconditions_overview(fid,idom_hdf)
            dname_trim = trim(adjustl(dname)) 
            call print_boundaryconditions_domain_face(dname_trim(3:), trim(adjustl(faces(iface))))


            !
            ! Check current face boundary condition
            !
            if ( bcnames(iface) == 'empty' ) then

            else
                call print_boundarycondition_properties(bcface)
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


                !
                ! Set boundary condition case
                !
                case (1)
                    !
                    ! Refresh display
                    !
                    call execute_command_line("clear")
                    call print_overview(fid,idom_hdf)
                    call print_boundaryconditions_overview(fid,idom_hdf)
                    dname_trim = trim(adjustl(dname)) 
                    call print_boundaryconditions_domain_face(dname_trim(3:), trim(adjustl(faces(iface))))


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



                !
                ! Edit property case
                !
                case (2)

                    get_property = .true.
                    do while ( get_property )
                        !
                        ! Refresh display
                        !
                        call execute_command_line("clear")
                        call print_overview(fid,idom_hdf)
                        call print_boundaryconditions_overview(fid,idom_hdf)
                        dname_trim = trim(adjustl(dname)) 
                        call print_boundaryconditions_domain_face(dname_trim(3:), trim(adjustl(faces(iface))))
                        call print_boundarycondition_properties(bcface)

                        !
                        ! Call edit option
                        !
                        command = "Enter boundary condition property: "
                        call write_line(command,color='blue')
                        read(*,*) pname

                        !
                        ! Check property exists
                        !
                        property_exists = check_property_exists(bcface,pname)

                        if ( property_exists ) get_property = .false.

                    end do 



                    call edit_property(bcface,pname)



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
        call delete_boundarycondition_hdf(bcgroup,bcface)


        !
        ! Add new boundary condition and settings
        !
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

        !
        ! Set boundary condition name attribute: bc_name
        !
        call h5ltset_attribute_string_f(bcface, ".", "bc_name", bcstring, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5ltset_attribute_string")


        if ( trim(bcstring) == 'empty' ) then
            !
            ! If 'empty' do not allocate new bc
            !
            
        else
            !
            ! Create an instance of the specified boundary condition to query its options
            !
            call create_bc(bcstring,bc)
        

            !
            ! Add boundary condition functions
            !
            call add_bcproperties_hdf(bcface,bc)
        end if


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
        character(len=:),   allocatable :: option_key, fcn_name
        real(rk)                        :: option_value


        !
        ! Get number of functions in the boundary condition
        !
        nprop = bc%get_nproperties()


        !
        ! Loop through and add properties
        !
        do iprop = 1,nprop

            !
            ! Get string the property is associated with
            !
            pstring = bc%get_property_name(iprop)


            !
            ! Create a new group for the property
            !
            call h5gcreate_f(bcface, "BCP_"//trim(adjustl(pstring)), prop_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"add_bcproperties_hdf: error creating new group for bcfunction")


            !
            ! Print property function attribute
            !
            fcn_name = bc%bcproperties%bcprop(iprop)%fcn%get_name()
            call h5ltset_attribute_string_f(prop_id, ".", "function", fcn_name, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"add_bcproperties_hdf: error setting function attribute")


            !
            ! Get number of options available for the current property
            !
            nopt = bc%get_noptions(iprop)


            if (nopt > 0 ) then
                do iopt = 1,nopt

                    !
                    ! Get the current option and default value.
                    !
                    option_key   = bc%get_option_key(iprop,iopt)
                    option_value = bc%get_option_value(iprop,option_key)

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

        integer(HSIZE_T)                        :: iattr, idx
        integer(ik)                             :: nattr
        integer                                 :: nmembers, igrp, type, ierr, iter, iprop, nprop
        character(len=10)                       :: faces(NFACES)
        character(len=1024),    allocatable     :: anames(:), pnames(:)
        character(len=1024)                     :: gname
        type(h5o_info_t), target                :: h5_info


        !
        ! Delete overall boundary condition face attributes
        !
        call delete_group_attributes(bcface)


        !
        !  Get number of groups linked to the current bcface
        !
        call h5gn_members_f(bcface, ".", nmembers, ierr)


        !
        !  Loop through groups and delete properties
        !
        if ( nmembers > 0 ) then

            !
            ! First get number of properties. This could be different than number of groups.
            !
            nprop = 0
            do igrp = 0,nmembers-1
                !
                ! Get group name
                !
                call h5gget_obj_info_idx_f(bcface, ".", igrp, gname, type, ierr)

                !
                ! Test if group is a boundary condition function. 'BCP_'
                !
                if (gname(1:4) == 'BCP_') then

                    ! increment nprop
                    nprop = nprop + 1

                end if

            end do  ! igrp


            !
            ! Second, get all property names
            !
            allocate(pnames(nprop), stat=ierr)
            if (ierr /= 0) call AllocationError
            iprop = 1
            do igrp = 0,nmembers-1
                call h5gget_obj_info_idx_f(bcface, ".", igrp, gname, type, ierr)

                !
                ! Test if group is a boundary condition function. 'BCP_'
                !
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

                call delete_bcproperty_hdf(bcface,pnames(iprop))
            
            end do


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

















    !>  Print the properties associated with a given boundary condition face group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @param[in]  bcface  HDF5 group identifier. Should be associated with a boundary condition.
    !!
    !--------------------------------------------------------------------------------------------------------
    subroutine print_boundarycondition_properties(bcface)
        integer(HID_T),     intent(in)  :: bcface

        integer(HID_T)                          :: bcprop
        integer                                 :: nmembers, igrp, type, ierr
        character(len=1024)                     :: gname, bcname



        !
        !  Print bc name
        !
        call h5ltget_attribute_string_f(bcface, ".", "bc_name", bcname, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5ltset_attribute_string")
        call write_line('['//trim(adjustl(bcname))//']')



        !
        !  Get number of groups linked to the current bcface
        !
        call h5gn_members_f(bcface, ".", nmembers, ierr)


        !
        !  Loop through groups and delete properties
        !
        if ( nmembers > 0 ) then
            do igrp = 0,nmembers-1


                !
                ! Get group name
                !
                call h5gget_obj_info_idx_f(bcface, ".", igrp, gname, type, ierr)

                !
                ! Test if group is a boundary condition function. 'BCP_'
                !
                if (gname(1:4) == 'BCP_') then
                    !
                    ! Open bcproperty_t group
                    !
                    call h5gopen_f(bcface, gname, bcprop, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"print_boundarycondition_properties: error opening bcproperty")


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



    end subroutine print_boundarycondition_properties
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
        real(rk),   dimension(1)                :: buf



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
                    option_vals(iattr) = buf(1)
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
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine edit_property(bcface,pname)
        integer(HID_T),     intent(in)  :: bcface
        character(*),       intent(in)  :: pname

        integer(HID_T)                  :: bcprop
        integer(HSIZE_T)                :: adim
        character(len=:),   allocatable :: command
        character(len=1024)             :: option, new_function
        logical                         :: option_exists, run
        real(rk)                        :: val
        integer                         :: ierr

        !
        ! Open bcproperty_t group
        !
        call h5gopen_f(bcface, "BCP_"//trim(pname), bcprop, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"edit_property: error opening bcproperty")



        !
        ! Edit option loop
        !
        run = .true.
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
                    call h5ltset_attribute_double_f(bcface, "BCP_"//trim(pname), trim(option), [val], adim, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"edit_property: error setting option value")

                end if

                !
                ! Close routine
                !
                run = .false.
            else

                call write_line("Invalid option", color='blue')

            end if


        end do ! run


        !
        ! Close bcproperty
        !
        call h5gclose_f(bcprop,ierr)

    end subroutine edit_property
    !********************************************************************************************












    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function check_property_exists(bcface,pname) result(exist_status)
        integer(HID_T),     intent(in)  :: bcface
        character(*),       intent(in)  :: pname

        integer :: ierr
        logical :: exist_status

        
        !
        ! Call object exist check
        !
        call h5lexists_f(bcface, "BCP_"//trim(pname), exist_status, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"check_property_exists: error check if property exists.")


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
            call h5ltset_attribute_double_f(bcprop, ".", trim(option), [val], adim, ierr)

        end do ! iopt

    end subroutine set_function
    !********************************************************************************************





end module mod_chidg_edit_boundaryconditions
