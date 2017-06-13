!>  Clone a ChiDG file configuration from a source file to a target file.
!!
!!  @author Nathan A. Wukie (AFRL)
!!  @date   6/13/2017
!!
!! Usage:   chidg clone source_file.h5 target_file.h5
!!
!!
!---------------------------------------------------------------------------------------------
module mod_chidg_clone
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_hdf_utilities,  only: copy_bc_state_groups_hdf, copy_patches_attributes_hdf, &
                                  open_file_hdf, close_file_hdf
    use hdf5
    implicit none







contains




    !>  Call copy configuration routine in high-order ChiDG library.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/13/2017
    !!
    !-----------------------------------------------------------------------------------
    subroutine chidg_clone(source_file,target_file)
        character(*),   intent(in)  :: source_file
        character(*),   intent(in)  :: target_file

        integer(HID_T)  :: fid_a, fid_b
        integer(ik)     :: user_input


        call execute_command_line("clear")


        !
        ! Write header and instructions...
        !
        call write_line("-------------------------------------------------------------------------------------------",color='blue')
        call write_line(" ", ltrim=.false.,color='blue')
        call write_line("                Clone configuration data from one ChiDG file to another.                   ",ltrim=.false.,color='blue')
        call write_line(" ", ltrim=.false.)
        call write_line("Copying from:", trim(source_file), color='blue', ltrim=.false.)
        call write_line("Copying to:  ", trim(target_file), color='blue', ltrim=.false.)
        call write_line("-------------------------------------------------------------------------------------------",color='blue')
        call write_line(" ", ltrim=.false.)
        call write_line("Mode", "Description", columns=.true.,column_width=42,color='blue')
        call write_line("-------------------------------------------------------------------------------------------")
        call write_line("1", "Copy all.",columns=.true.,column_width=42)
        call write_line("2", "Copy only boundary condition state groups.",column_width=42,columns=.true.)
        call write_line("3", "Copy only patch attributes.",column_width=42,columns=.true.)
        call write_line("0", "Exit", column_width=42, columns=.true.)
        call write_line("-------------------------------------------------------------------------------------------")
        call write_line("Enter mode: ")


        !
        ! Accept user input...
        ! 
        read(*,*) user_input



        !
        ! Open files
        !
        fid_a = open_file_hdf(source_file)
        fid_b = open_file_hdf(target_file)



        select case(user_input)
            case(1)
                call copy_bc_state_groups_hdf(fid_a,fid_b)
                call copy_patches_attributes_hdf(fid_a,fid_b)

            case(2)
                call copy_bc_state_groups_hdf(fid_a,fid_b)

            case(3)
                call copy_patches_attributes_hdf(fid_a,fid_b)

        end select
                

        call close_file_hdf(fid_a)
        call close_file_hdf(fid_b)


        call execute_command_line("clear")

    end subroutine chidg_clone
    !***********************************************************************************




end module mod_chidg_clone
