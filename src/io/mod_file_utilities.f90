module mod_file_utilities
#include <messenger.h>
    use mod_string,             only: get_file_extension
    use type_file_properties,   only: file_properties_t
    implicit none


contains



    !>  Copy source file to target file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !---------------------------------------------------------------------------------------------------------
    subroutine copy_file(sourcefile, targetfile)
        character(*),   intent(in)  :: sourcefile
        character(*),   intent(in)  :: targetfile

        integer             :: unit_src, unit_tar, ierr, irec
        character(len=1)    :: c

        ! Open files for source and target.
        open(newunit=unit_tar, file=targetfile, access='stream', action='write', iostat=ierr)
        if ( ierr /= 0 ) call chidg_signal(FATAL,"copy_file: error opening target file.")

        open(newunit=unit_src, file=sourcefile, access='stream', action='read',  iostat=ierr)
        if ( ierr /= 0 ) call chidg_signal(FATAL,"copy_file: error opening source file.")

        ! While there is a character coming from the source file, get character write to target file
        ierr = 0
        do while ( ierr == 0 )
            ! Get from source file
            read(unit=unit_src, iostat=ierr) c

            ! Write to target file
            if ( ierr == 0 ) then
                write(unit_tar) c
            end if
        end do

        ! Close files
        close(unit_src)
        close(unit_tar)

    end subroutine copy_file
    !*************************************************************************************************************



    !>  Delete file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine delete_file(sourcefile)
        character(*),   intent(in)  :: sourcefile

        integer     :: unit_src, ierr, irec
        
        ! Open source file for deletion.
        open(newunit=unit_src, file=sourcefile, status='old', iostat=ierr)

        ! Close and delete.
        if (ierr == 0) close(unit_src, status='delete')
        
    end subroutine delete_file
    !*************************************************************************************************************




    !>  Text file header
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2017
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine write_file_header(sourcefile,file_exists,func_name,ref_geom,aux_geom)
        character(*),   intent(in)  :: sourcefile
        logical,        intent(in)  :: file_exists
        character(*),   intent(in)  :: func_name
        character(*),   intent(in)  :: ref_geom
        character(*),   intent(in)  :: aux_geom

        integer     :: unit_src, ierr, irec
        
        ! Open file: create a new file if it does not exist, otherwise replace it
        if (file_exists) then
            open(newunit=unit_src, file=sourcefile, status='replace', iostat=ierr)
        else
            open(newunit=unit_src, file=sourcefile, status='new', iostat=ierr)
        end if

        if ( ierr /= 0 ) call chidg_signal(FATAL,"write_file_header: error opening target file.")

        write(unit_src, '(2A)') 'Functional name:    ', trim(func_name)
        write(unit_src, '(2A)') 'Reference geometry: ', trim(ref_geom)
        write(unit_src, '(2A)') 'Auxiliary geometry: ', trim(aux_geom)
        write(unit_src, '(A)' ) 'Computed data [value,step,time]:'

       close(unit_src) 

    end subroutine write_file_header
    !*************************************************************************************************************




end module mod_file_utilities
