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
    use mod_chidg_edit_matrixsolver,        only: chidg_edit_matrixsolver
    use mod_chidg_edit_timescheme,          only: chidg_edit_timescheme
    use mod_chidg_edit_printoverview,       only: print_overview
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


        !
        ! Initialize chidg environment
        !
        call chidg%start_up('mpi')
        call chidg%start_up('core')

        !
        ! Send ChiDG output to screen
        !
        IO_DESTINATION = 'screen'


        !call check_extension(file,'.h5')

        !
        ! Clear screen for editing
        ! TODO: portability check
        !
        call execute_command_line("clear")




        fid = open_file_hdf(filename)
!        !
!        ! Check that file can be found
!        !
!        inquire(file=filename, exist=fileexists)
!        if ( .not. fileexists ) then
!            call chidg_signal_one(FATAL,"chidg_edit: file not found for editing.",filename)
!        end if
!
!
!        !
!        ! Initialize Fortran interface
!        !
!        call h5open_f(ierr)
!        if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit: HDF5 Fortran interface had an error during initialization.")
!
!
!        !
!        ! Get HDF file identifier
!        !
!        call h5fopen_f(filename, H5F_ACC_RDWR_F, fid, ierr)
!        if (ierr /= 0) call chidg_signal_one(FATAL,"chidg_edit: Error opening HDF5 file for editing.",filename)



        !
        ! Edit loop
        !
        run = .true.
        do while ( run )

            !
            ! Refresh display with overview of file contents
            !
            call execute_command_line("clear")
            call print_overview(fid)


            !
            ! Print command options, accept user selection.
            !
            call write_line(' ')
            call write_line("Select command: ")
            call write_line("1:domain info","2:boundary conditions","3:time scheme","4:matrix solver","0:exit", columns=.True., column_width=20, color='blue')

            ierr = 1
            do while ( ierr /= 0 )
                read(*,'(I8)', iostat=ierr) int_input
                if ( (ierr/=0) .or. (abs(int_input)>3) ) print*, "Invalid input: expecting 0, 1, 2, or 3."
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
                    call chidg_edit_timescheme(fid)
                case (4)
                    call chidg_edit_matrixsolver(fid)
                case default

            end select


        end do







        !
        ! Close HDF5 file and Fortran interface
        !
        !call h5fclose_f(fid, ierr)  ! Close HDF5 File
        !if (ierr /= 0) call chidg_signal(FATAL,"chidg_edit: error closing file.")
        !call h5close_f(ierr)        ! Close HDF5 Fortran interface
        call close_file_hdf(fid)



        !
        ! Clear terminal upon exit
        !
        call execute_command_line("clear")

    end subroutine chidg_edit
    !********************************************************************************************



















end module mod_chidg_edit
