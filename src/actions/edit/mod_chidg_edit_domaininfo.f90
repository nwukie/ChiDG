module mod_chidg_edit_domaininfo
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: NFACES, XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use mod_equations,  only: equation_builder_factory
    use hdf5
    use h5lt

    use mod_hdf_utilities,              only: get_ndomains_hdf, get_domain_names_hdf, get_domain_name_hdf,  &
                                              open_domain_hdf, close_domain_hdf, set_domain_equation_set_hdf, &
                                              check_domain_exists_hdf, set_domain_name_hdf
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
    subroutine chidg_edit_domaininfo(fid)
        integer(HID_T),     intent(in)  :: fid


        integer(ik)                     :: ierr
        logical                         :: run, open_domain, domain_exists
        character(1024)                 :: domain_name
        character(len=:),   allocatable :: bc_commands



        run = .true.
        do while ( run )

            !
            ! Refresh display
            !
            call execute_command_line("clear")
            call print_overview(fid)



            !
            ! Print command options, accept user selection.
            !
            bc_commands = "Enter a domain to edit(Enter blank to exit):"
            call write_line(' ')
            call write_line(bc_commands,color='blue')



            read(*,'(A1024)') domain_name

            domain_exists = check_domain_exists_hdf(fid,trim(domain_name))

            if (trim(domain_name) == '') then
                run         = .false.
                open_domain = .false.
            else if ( (trim(domain_name) /= '' ) .and. &
                      (domain_exists .eqv. .false.) ) then
                run         = .true.
                open_domain = .false.
            else
                run         = .true.
                open_domain = .true.
            end if



            if (open_domain) call chidg_edit_domaininfo_domain(fid,trim(domain_name))




        end do  ! run_bc_edit





    end subroutine chidg_edit_domaininfo
    !************************************************************************************************
















    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------
    subroutine chidg_edit_domaininfo_domain(fid,domain_name)
        integer(HID_T),     intent(in)  :: fid
        character(*),       intent(in)  :: domain_name

        integer(HID_T)                  :: domain_id
        integer(ik)                     :: ierr, iedit
        logical                         :: run_domain_edit
        character(len=:),   allocatable :: domain_commands




        domain_id = open_domain_hdf(fid,trim(domain_name))



        run_domain_edit = .true.
        do while ( run_domain_edit )

            !
            ! Refresh display
            !
            call execute_command_line("clear")
            call print_overview(fid,domain_id)



            !
            ! Print command options, accept user selection.
            !
            call write_line(' ')
            call write_line("1: Edit name", "2: Edit equation set", "0: exit",color='blue')
            ierr = 1
            do while ( ierr /= 0 )
                read(*,'(I8)', iostat=ierr) iedit
                if ( ierr /= 0 )  call write_line("Invalid input: expecting an integer index.")

                if ( (iedit /= 0) .and. (iedit /= 1) .and. (iedit /= 2 ) ) then
                    ierr = 1
                    call write_line("Invalid input. Enter 0, 1, or 2.")
                end if
            end do



            !
            ! Operate on particular domain
            !
            select case (iedit)
                case (0)
                    run_domain_edit = .false.

                case (1) ! edit name
                    call chidg_edit_domaininfo_domain_name(fid,domain_id)

                case (2) ! edit equation set
                    call chidg_edit_domaininfo_domain_equation_set(fid,domain_id)
                
                case default
                    run_domain_edit = .false.

            end select


        end do  ! run_domain_edit



        call close_domain_hdf(domain_id)


    end subroutine chidg_edit_domaininfo_domain
    !************************************************************************************************

















    !>  Reset the name of a specified domain in a ChiDG HDF file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF5 file identifier
    !!  @param[in]  idom    Domain index to be renamed
    !!
    !------------------------------------------------------------------------------------------------
    subroutine chidg_edit_domaininfo_domain_name(fid,domain_id)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: domain_id


        integer(ik)                         :: ierr, idom, iind
        character(len=1024)                 :: dname_current, dname_new



        dname_current = get_domain_name_hdf(domain_id)


        !
        ! Refresh display
        !
        call execute_command_line("clear")
        call print_overview(fid,domain_id)


        !
        ! Print command options, accept user selection.
        !
        call write_line(' ')
        call write_line("Enter new domain name: ",color='blue')
        read(*,'(A1024)') dname_new


        !
        ! Move link to rename
        !
        if ( (trim(dname_new) /= '') .and. (trim(dname_new) /= trim(dname_current)) ) then
            !call h5gmove_f(fid, trim(adjustl(dname_current)), "D_"//trim(adjustl(dname_new)), ierr)
            call h5lmove_f(fid, "D_"//trim(adjustl(dname_current)), fid, "D_"//trim(adjustl(dname_new)), ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_domaininfo_domain_name: error renaming domain")

            call set_domain_name_hdf(domain_id,trim(dname_new))
        end if

    end subroutine chidg_edit_domaininfo_domain_name
    !************************************************************************************************














    !>  Reset the name of a specified domain in a ChiDG HDF file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF5 file identifier
    !!  @param[in]  idom    Domain index to be renamed
    !!
    !------------------------------------------------------------------------------------------------
    subroutine chidg_edit_domaininfo_domain_equation_set(fid,domain_id)
        integer(HID_T),     intent(in)  :: fid
        integer(HID_T),     intent(in)  :: domain_id


        logical                             :: list_equations
        logical                             :: equation_not_found
        character(:),           allocatable :: not_found_string
        integer(ik)                         :: ierr, idom, iind
        integer(HID_T)                      :: did
        character(len=1024), allocatable    :: dnames(:)
        character(len=1024)                 :: equation_set
        character(len=1024)                 :: dname
        integer(ik),        allocatable     :: dindices(:)


        list_equations = .false.
        equation_not_found = .false.
        do

            !
            ! Refresh display
            !
            call execute_command_line("clear")
            call print_overview(fid,domain_id)



            !
            ! Print command options, accept user selection.
            !
            call write_line(' ')
            call write_line("Enter equation set(? to list): ",color='blue')


            if (list_equations) call equation_builder_factory%list()

            not_found_string = "We didn't find '"//trim(equation_set)//"' registered &
                                in ChiDG :/. Try another equation and remember you can &
                                enter '?' to list the registered equation sets."
            if (equation_not_found) call write_line(not_found_string)



            read(*,'(A1024)', iostat=ierr) equation_set
            if ( (ierr/=0)  ) print*, "Invalid input. Try again :)"




            ! List equations
            if (trim(equation_set) == '?') then
                list_equations = .true.
                equation_not_found = .false.

            ! Exit without changes
            else if (trim(equation_set) == '') then
                exit

            ! Check registered, exit if valid
            else
                if ( equation_builder_factory%has(trim(equation_set)) ) exit
                list_equations     = .false.
                equation_not_found = .true.
            end if

        end do


        !
        ! Set equation set
        !
        if (trim(equation_set) /= '') then
            call set_domain_equation_set_hdf(domain_id,trim(adjustl(equation_set)))
        end if

    end subroutine chidg_edit_domaininfo_domain_equation_set
    !************************************************************************************************







































end module mod_chidg_edit_domaininfo
