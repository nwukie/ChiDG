module mod_chidg_edit_boundaryconditions
#include <messenger.h>
    use mod_kinds,  only: rk, ik
    use hdf5
    use h5lt

    use mod_hdf_utilities,              only: get_ndomains_hdf
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


        integer(ik)                     :: ierr, int_input, ndom
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
            !print*, bc_commands
            ierr = 1
            do while ( ierr /= 0 )
                read(*,'(I8)', iostat=ierr) int_input
                if ( ierr /= 0 )  print*, "Invalid input: expecting an integer index."

                if ( int_input > ndom ) then
                    ierr = 1
                    print*, "Invalid domain range. Enter a number between 1 and ",ndom
                end if

            end do



            !
            ! Operate on particular domain
            !
            select case (int_input)
                case (0)
                    run_bc_edit = .false.

                case default
                    !call print_overview_boundaryconditions(fid,int_input)

            end select








        end do  ! run_bc_edit





    end subroutine chidg_edit_boundaryconditions
    !************************************************************************************************















end module mod_chidg_edit_boundaryconditions
