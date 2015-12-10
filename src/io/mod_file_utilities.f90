module mod_file_utilities
#include <messenger.h>
    use mod_hdf_utilities,      only: get_properties_hdf
    use mod_string_utilities,   only: get_file_extension
    use type_file_properties,   only: file_properties_t
    implicit none




contains



    !>
    !!
    !!
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------
    function get_file_properties(filename) result(file_props)
        character(*),   intent(in)  :: filename

        type(file_properties_t)         :: file_props
        character(len=:), allocatable   :: extension
        character(len=5), allocatable   :: extensions(:)


        extensions = ['.h5']

        !
        ! Get file extension
        !
        extension = get_file_extension(filename, extensions)


        print*, filename, extension

        !
        ! Call specialized routine for returning file properties
        !
        if ( extension == '.h5' ) then
            file_props = get_properties_hdf(filename)
        else
            call chidg_signal(FATAL, "Error: get_fileproperties -- file extension not recognized")
        end if



    end function get_file_properties
    !#####################################################################################















    !> Copy source file to target file.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine copy_file(sourcefile, targetfile)
        character(*),   intent(in)  :: sourcefile
        character(*),   intent(in)  :: targetfile

        integer     :: unit_src, unit_tar, ierr, irec
        character   :: char

        
        !
        ! Open files for source and target.
        !
        open(newunit=unit_src, file=sourcefile, access='direct', status='old', action='read', iostat=ierr, recl=1)
        open(newunit=unit_tar, file=targetfile, access='direct', status='replace', action='write', iostat=ierr, recl=1)

        
        !
        ! Copy file by reading from source, character-by-character and writing to target.
        !
        irec=1
        do
            read(unit=unit_src, rec=irec, iostat=ierr) char
            if (ierr /= 0) exit

            write(unit=unit_tar, rec=irec) char
            irec = irec + 1
        end do


        !
        ! Close files.
        !
        close(unit_src)
        close(unit_tar)

    end subroutine copy_file
    !#######################################################################################














    !> Delete file.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine delete_file(sourcefile)
        character(*),   intent(in)  :: sourcefile

        integer     :: unit_src, ierr, irec
        
        !
        ! Open source file for deletion.
        !
        open(newunit=unit_src, file=sourcefile, status='old', iostat=ierr)

        !
        ! Close and delete.
        !
        if (ierr == 0) close(unit_src, status='delete')
        

    end subroutine delete_file
    !#######################################################################################


























end module mod_file_utilities
