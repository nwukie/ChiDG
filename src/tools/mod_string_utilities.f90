module mod_string_utilities






contains





    !> Given a file name and a list of possible extensions, return an extension if 
    !! it is contained in the file string.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
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
    !**********************************************************************************************







    !>  Converts a character string to all lower-case by comparing the index of each character in the
    !!  ASCII set and shifting the index to the appropriate case of the character.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    function string_to_lower(string) result(string_lower)
        character(*),   intent(in)  :: string

        character(len=len(string))  :: string_lower
        integer                         :: i, k

        do i = 1,len(trim(string))
            
            k = iachar(string(i:i))


            if ( k >= iachar('A') .and. k <= iachar('Z') ) then
                k = k + (iachar('a') - iachar('A'))
                string_lower(i:i) = achar(k)
            else
                string_lower(i:i) = string(i:i)
            end if

        end do ! i
    end function string_to_lower
    !************************************************************************************************








    !>  Converts a character string to all upper-case by comparing the index of each character in the
    !!  ASCII set and shifting the index to the appropriate case of the character.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    function string_to_upper(string) result(string_upper)
        character(*),   intent(in)  :: string

        character(len=len(string))  :: string_upper
        integer                         :: i, k

        do i = 1,len(trim(string))
            
            k = iachar(string(i:i))


            if ( k >= iachar('a') .and. k <= iachar('z') ) then
                k = k - (iachar('a') - iachar('A'))
                string_upper(i:i) = achar(k)
            else
                string_upper(i:i) = string(i:i)
            end if

        end do ! i
    end function string_to_upper
    !************************************************************************************************










end module mod_string_utilities
