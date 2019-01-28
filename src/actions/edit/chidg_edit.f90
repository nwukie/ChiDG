!>
!!
!!  @author Nathan A. Wukie
!!  @date   2/3/2016
!!
!!
!!
!!
!----------------------------------------------------------------------------------------
module mod_chidg_edit
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_hdf_utilities,  only: open_file_hdf, close_file_hdf
    use hdf5
    use h5lt

    use type_chidg,                         only: chidg_t
    use mod_chidg_edit_domaininfo,          only: chidg_edit_domaininfo
    use mod_chidg_edit_boundaryconditions,  only: chidg_edit_boundaryconditions
    use mod_chidg_edit_meshmotion,          only: chidg_edit_meshmotion
    use mod_chidg_edit_matrixsolver,        only: chidg_edit_matrixsolver
    use mod_chidg_edit_timescheme,          only: chidg_edit_timescheme
    use mod_chidg_edit_printoverview,       only: print_overview, chidg_clear_screen
    implicit none




contains



    !>  ChiDG action: edit utility
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @parma[in]  file    Character string specifying an .h5 file in ChiDG format for editing.
    !!
    !--------------------------------------------------------------------------------------------
    subroutine chidg_edit(filename)
        character(*),   intent(in)  :: filename

        type(chidg_t)                   :: chidg

        logical                         :: run, fileexists
        integer(ik)                     :: ierr

        character(len=:),   allocatable :: char_input
        integer(ik)                     :: int_input
        real(rk)                        :: real_input

        character(len=:),   allocatable :: command_options

        integer(HID_T)     :: fid

        ! Start Alternate Screen so terminal state can be restored afterwards
        print*, achar(27)//"[?1049h"

        !
        ! Initialize chidg environment
        !
        call chidg%start_up('mpi')
        call chidg%start_up('core')

        !
        ! Send ChiDG output to screen
        !
        IO_DESTINATION = 'screen'



        !
        ! Clear screen for editing
        ! TODO: portability check
        !
        call chidg_clear_screen()



        !
        ! Edit loop
        !
        fid = open_file_hdf(filename)
        run = .true.
        do while ( run )

            !
            ! Refresh display with overview of file contents
            !
            call chidg_clear_screen()
            call print_overview(fid)


            !
            ! Print command options, accept user selection.
            !
            call write_line(' ')
            call write_line("Select command: ")
            call write_line("1:domain info","2:boundary conditions","0:exit", columns=.True., column_width=20, color='blue')

            ierr = 1
            do while ( ierr /= 0 )
                read(*,'(I8)', iostat=ierr) int_input
                if ( (ierr/=0) .or. (abs(int_input)>3) ) print*, "Invalid input: expecting 0, 1, or 2."
            end do

            
            !
            ! Sub-menu functions. Enter based on user-selection.
            !
            select case (int_input)
                case (0)
                    exit
                case (1)
                    call chidg_edit_domaininfo(fid)
                case (2)
                    call chidg_edit_boundaryconditions(fid)
                case (3)
                    call chidg_edit_meshmotion(fid)
                !case (3)
                !    call chidg_edit_timescheme(fid)
                !case (4)
                !    call chidg_edit_matrixsolver(fid)
                case default

            end select


        end do



        !
        ! Close HDF5 file and Fortran interface
        !
        call close_file_hdf(fid)



        !
        ! Clear terminal upon exit
        !
        call chidg_clear_screen()


        ! Restore terminal state from Alternate Screen
        print*, achar(27)//"[?1049l"


    end subroutine chidg_edit
    !********************************************************************************************













end module mod_chidg_edit
