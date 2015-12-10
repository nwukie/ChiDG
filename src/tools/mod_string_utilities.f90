module mod_string_utilities






contains





    !> Given a file name and a list of possible extensions, return an extension if 
    !! it is contained in the file string.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  file        Character string containing a file name
    !!  @param[in]  extensions  Array of character strings. Each string is an accepted extension
    !!
    !---------------------------------------------------------------------------------------------
    function get_file_extension(file,extensions) result(extension)
        character(*),       intent(in)  :: file
        character(len=5),   intent(in)  :: extensions(:)

        character(len=:), allocatable   :: extension
        integer                         :: iext, extloc


        do iext = 1,size(extensions)

            !
            ! Check for extension in grid file
            !
            extloc = index(file, trim(extensions(iext)))

            !
            ! If extloc is nonzero, then the extension was found in the filename
            !
            if ( extloc /= 0 ) then
                extension = trim(extensions(iext))
                exit
            end if

        end do


    end function get_file_extension
    !##############################################################################################












end module mod_string_utilities
