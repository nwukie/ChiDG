module mod_chidg_edit_domaininfo
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: NFACES, XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use hdf5
    use h5lt

    use mod_hdf_utilities,              only: get_ndomains_hdf, get_domain_names_hdf, get_bcnames_hdf, get_domain_indices_hdf
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


        integer(ik)                     :: ierr, idom_hdf, ndom
        logical                         :: run_bc_edit
        character(len=:),   allocatable :: bc_commands



        ndom = get_ndomains_hdf(fid)



        bc_commands = "Select a domain for editing(0 to exit):"


        run_bc_edit = .true.
        do while ( run_bc_edit )

            !
            ! Refresh display
            !
            call execute_command_line("clear")
            call print_overview(fid)



            !
            ! Print command options, accept user selection.
            !
            call write_line(' ')
            call write_line(bc_commands,color='blue')
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
                    call chidg_edit_domaininfo_domain(fid,idom_hdf)

            end select








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
    subroutine chidg_edit_domaininfo_domain(fid,idom_hdf)
        integer(HID_T),     intent(in)  :: fid
        integer(ik),        intent(in)  :: idom_hdf


        integer(ik)                     :: ierr, ndom, iedit
        logical                         :: run_domain_edit
        character(len=:),   allocatable :: domain_commands



        ndom = get_ndomains_hdf(fid)





        run_domain_edit = .true.
        do while ( run_domain_edit )

            !
            ! Refresh display
            !
            call execute_command_line("clear")
            call print_overview(fid,idom_hdf)



            !
            ! Print command options, accept user selection.
            !
            call write_line(' ')
            call write_line("1: edit domain name", "2: edit coordinate expansion", "0: exit",color='blue')
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
                    call chidg_edit_domaininfo_domain_name(fid,idom_hdf)

                case (2) ! edit coordinate expansion order. 'mapping'
                    call chidg_edit_domaininfo_domain_gridorder(fid,idom_hdf)
                
                case default
                    run_domain_edit = .false.

            end select


        end do  ! run_domain_edit





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
    subroutine chidg_edit_domaininfo_domain_name(fid,idom_hdf)
        integer(HID_T),     intent(in)  :: fid
        integer(ik),        intent(in)  :: idom_hdf


        integer(ik)                         :: ierr, idom, iind
        character(len=1024), allocatable    :: dnames(:)
        character(len=1024)                 :: dname_current, dname_new
        integer(ik),        allocatable     :: dindices(:)


        dindices      = get_domain_indices_hdf(fid)
        dnames        = get_domain_names_hdf(fid)


        !
        ! Find idom_hdf in the list
        !
        do iind = 1,size(dindices)
            if ( dindices(iind) == idom_hdf ) then

                dname_current = dnames(iind)

            end if
        end do


        !
        ! Refresh display
        !
        call execute_command_line("clear")
        call print_overview(fid,idom_hdf)


        !
        ! Print command options, accept user selection.
        !
        call write_line(' ')
        call write_line("Enter new domain name: ",color='blue')
        read(*,*) dname_new


        !
        ! Modify domain name
        !
        call h5gmove_f(fid, trim(adjustl(dname_current)), "D_"//trim(adjustl(dname_new)), ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_domaininfo_domain_name: error renaming domain")


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
    subroutine chidg_edit_domaininfo_domain_gridorder(fid,idom_hdf)
        integer(HID_T),     intent(in)  :: fid
        integer(ik),        intent(in)  :: idom_hdf


        integer(HSIZE_T)                    :: adim
        integer(ik)                         :: ierr, idom, iind
        character(len=1024), allocatable    :: dnames(:)
        integer(ik)                         :: grid_order_new
        character(len=1024)                 :: dname
        integer(ik),        allocatable     :: dindices(:)


        !
        ! Refresh display
        !
        call execute_command_line("clear")
        call print_overview(fid,idom_hdf)


        !
        ! Find idom_hdf in the list
        !
        dindices = get_domain_indices_hdf(fid)
        dnames   = get_domain_names_hdf(fid)

        do iind = 1,size(dindices)
            if ( dindices(iind) == idom_hdf ) then

                dname = dnames(iind)

            end if
        end do



        !
        ! Print command options, accept user selection.
        !
        call write_line(' ')
        call write_line("Enter order of coordinate expansion: ",color='blue')
        read(*,*) grid_order_new



        !
        ! Modify grid mapping attribute
        !
        adim = 1
        call h5ltset_attribute_int_f(fid, trim(dname), 'mapping', [grid_order_new], adim, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit_domaininfo_domain_gridorder: error setting 'mapping' attribute.")



    end subroutine chidg_edit_domaininfo_domain_gridorder
    !************************************************************************************************


















end module mod_chidg_edit_domaininfo
