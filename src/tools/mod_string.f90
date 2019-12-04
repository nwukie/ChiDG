module mod_string
#include <messenger.h>


    type, public :: string_t
        character(len=:), allocatable   :: str
    contains
        procedure   :: lower
        procedure   :: upper
        procedure   :: set
        procedure   :: get
        procedure   :: strip
        procedure   :: replace
!        procedure, private :: write_formatted
!        generic :: write(formatted) => write_formatted
    end type string_t


!    interface string
!        module procedure new_string
!    end interface
    
    public str
    interface str
        module procedure to_string
    end interface


contains




    !>  Convert data to string representation and return string.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   11/20/2019
    !!
    !-------------------------------------------------------------------------------------------------
    function to_string(data_in) result(string_out)
        class(*),    intent(in)  :: data_in

        character(100)              :: write_internal
        character(:),   allocatable :: string_out

        string_out = ''

        select type(data_in)

            type is(character(len=*))
                string_out = data_in

            type is(integer)

                if ( data_in == 0 ) then
                    write_internal = '0'

                ! Positive integers
                else if ( data_in < 1e1 .and. data_in > 0) then
                    write(write_internal, '(I1.0)') data_in
                else if ( data_in < 1e2 .and. data_in > 0) then
                    write(write_internal, '(I2.0)') data_in
                else if ( data_in < 1e3 .and. data_in > 0) then
                    write(write_internal, '(I3.0)') data_in
                else if ( data_in < 1e4 .and. data_in > 0) then
                    write(write_internal, '(I4.0)') data_in
                else if ( data_in < 1e5 .and. data_in > 0) then
                    write(write_internal, '(I5.0)') data_in
                else if ( data_in < 1e6 .and. data_in > 0) then
                    write(write_internal, '(I6.0)') data_in
                else if ( data_in < 1e7 .and. data_in > 0) then
                    write(write_internal, '(I7.0)') data_in
                else if ( data_in < 1e8 .and. data_in > 0) then
                    write(write_internal, '(I8.0)') data_in
                else if ( data_in < 1e9 .and. data_in > 0) then
                    write(write_internal, '(I9.0)') data_in
                else if ( data_in < 1e10 .and. data_in > 0) then
                    write(write_internal, '(I10.0)') data_in
                else if ( data_in < 1e11 .and. data_in > 0) then
                    write(write_internal, '(I11.0)') data_in

                ! Negative integers: extra char width for neg sign
                else if ( abs(data_in) < 1e1 .and. data_in < 0) then
                    write(write_internal, '(I2.0)') data_in
                else if ( abs(data_in) < 1e2 .and. data_in < 0) then
                    write(write_internal, '(I3.0)') data_in
                else if ( abs(data_in) < 1e3 .and. data_in < 0) then
                    write(write_internal, '(I4.0)') data_in
                else if ( abs(data_in) < 1e4 .and. data_in < 0) then
                    write(write_internal, '(I5.0)') data_in
                else if ( abs(data_in) < 1e5 .and. data_in < 0) then
                    write(write_internal, '(I6.0)') data_in
                else if ( abs(data_in) < 1e6 .and. data_in < 0) then
                    write(write_internal, '(I7.0)') data_in
                else if ( abs(data_in) < 1e7 .and. data_in < 0) then
                    write(write_internal, '(I8.0)') data_in
                else if ( abs(data_in) < 1e8 .and. data_in < 0) then
                    write(write_internal, '(I9.0)') data_in
                else if ( abs(data_in) < 1e9 .and. data_in < 0) then
                    write(write_internal, '(I10.0)') data_in
                else if ( abs(data_in) < 1e10 .and. data_in < 0) then
                    write(write_internal, '(I11.0)') data_in
                else if ( abs(data_in) < 1e11 .and. data_in < 0) then
                    write(write_internal, '(I12.0)') data_in
                else
                    write(write_internal, '(I10.0)') data_in
                end if
                string_out = write_internal

            type is(integer(8))

                if ( data_in == 0 ) then
                    write_internal = '0'

                ! Positive integers
                else if ( data_in < 1e1 .and. data_in > 0) then
                    write(write_internal, '(I1.0)') data_in
                else if ( data_in < 1e2 .and. data_in > 0) then
                    write(write_internal, '(I2.0)') data_in
                else if ( data_in < 1e3 .and. data_in > 0) then
                    write(write_internal, '(I3.0)') data_in
                else if ( data_in < 1e4 .and. data_in > 0) then
                    write(write_internal, '(I4.0)') data_in
                else if ( data_in < 1e5 .and. data_in > 0) then
                    write(write_internal, '(I5.0)') data_in
                else if ( data_in < 1e6 .and. data_in > 0) then
                    write(write_internal, '(I6.0)') data_in
                else if ( data_in < 1e7 .and. data_in > 0) then
                    write(write_internal, '(I7.0)') data_in
                else if ( data_in < 1e8 .and. data_in > 0) then
                    write(write_internal, '(I8.0)') data_in
                else if ( data_in < 1e9 .and. data_in > 0) then
                    write(write_internal, '(I9.0)') data_in
                else if ( data_in < 1e10 .and. data_in > 0) then
                    write(write_internal, '(I10.0)') data_in
                else if ( data_in < 1e11 .and. data_in > 0) then
                    write(write_internal, '(I11.0)') data_in

                ! Negative integers: extra char width for neg sign
                else if ( abs(data_in) < 1e1 .and. data_in < 0) then
                    write(write_internal, '(I2.0)') data_in
                else if ( abs(data_in) < 1e2 .and. data_in < 0) then
                    write(write_internal, '(I3.0)') data_in
                else if ( abs(data_in) < 1e3 .and. data_in < 0) then
                    write(write_internal, '(I4.0)') data_in
                else if ( abs(data_in) < 1e4 .and. data_in < 0) then
                    write(write_internal, '(I5.0)') data_in
                else if ( abs(data_in) < 1e5 .and. data_in < 0) then
                    write(write_internal, '(I6.0)') data_in
                else if ( abs(data_in) < 1e6 .and. data_in < 0) then
                    write(write_internal, '(I7.0)') data_in
                else if ( abs(data_in) < 1e7 .and. data_in < 0) then
                    write(write_internal, '(I8.0)') data_in
                else if ( abs(data_in) < 1e8 .and. data_in < 0) then
                    write(write_internal, '(I9.0)') data_in
                else if ( abs(data_in) < 1e9 .and. data_in < 0) then
                    write(write_internal, '(I10.0)') data_in
                else if ( abs(data_in) < 1e10 .and. data_in < 0) then
                    write(write_internal, '(I11.0)') data_in
                else if ( abs(data_in) < 1e11 .and. data_in < 0) then
                    write(write_internal, '(I12.0)') data_in
                else
                    write(write_internal, '(I10.0)') data_in
                end if
                string_out = write_internal

            type is(real)
                if (abs(data_in) > 0.1) then
                    write(write_internal, '(F24.14)') data_in
                else
                    write(write_internal, '(E24.14)') data_in
                end if
                string_out = write_internal

            type is(real(8))
                if (abs(data_in) < 0.1) then
                    write(write_internal, '(E24.14)') data_in
                else if ( (abs(data_in) > 0.1) .and. (abs(data_in) < 1.e9) ) then
                    write(write_internal, '(F24.14)') data_in
                else
                    write(write_internal, '(E24.14)') data_in
                end if
                string_out = write_internal

            type is(logical(1))
                if (data_in) then
                    string_out = 'True'
                else
                    string_out = 'False'
                end if

            type is (logical(2))
                if (data_in) then
                    string_out = 'True'
                else
                    string_out = 'False'
                end if

            type is (logical(4))
                if (data_in) then
                    string_out = 'True'
                else
                    string_out = 'False'
                end if

            class default
                print*, 'Error: no IO rule for provided data in add_to_line'
                stop
        end select


    end function to_string
    !***********************************************************************************************
















    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine set(self,str)
        class(string_t),    intent(inout)   :: self
        character(len=*),   intent(in)      :: str

        self%str = str

    end subroutine set
    !*********************************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    function get(self) result(str)
        class(string_t),    intent(in)  :: self

        character(len=:),   allocatable :: str

        str = self%str

    end function get
    !*********************************************************************************************




    !>
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function lower(self) result(str_lower)
        class(string_t),    intent(in)  :: self

        character(len=:), allocatable   :: str_lower

        str_lower = string_to_lower(self%str)

    end function lower
    !*********************************************************************************************





    !>
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function upper(self) result(str_upper)
        class(string_t),    intent(in)  :: self

        character(len=:), allocatable   :: str_upper

        str_upper = string_to_upper(self%str)

    end function upper
    !*********************************************************************************************




    !> It reads in a file name and returns only the prefix, based on the expected extension 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   2/27/2017
    !!
    !!  @param[in]  file        Character string containing a file name
    !!  @param[in]  extensions  Character string of the expected extension
    !!
    !---------------------------------------------------------------------------------------------
    function get_file_prefix(file,extension) result(prefix)
        character(*),       intent(in)  :: file
        character(*),       intent(in)  :: extension

        character(len=:), allocatable   :: prefix
        integer                         :: iext, extloc
        
        !
        ! Check if the file name has already the extension
        !
        extloc = index(file, trim(extension))

        if ( extloc == 0 ) then
            prefix = file
        else
            prefix = trim(file(1:extloc-1))
        end if

    end function get_file_prefix
    !**********************************************************************************************










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






    !>  Split the string at the given separator and return an array of strings 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   2/12/2017
    !!
    !-----------------------------------------------------------------------------------------------
    function strip(self,separator) result(string_list)
        class(string_t),    intent(in)  :: self
        character(*),       intent(in)  :: separator
        
        type(string_t), allocatable     :: string_list(:)
        character(len=100)              :: temp_list(100)
        integer                         :: ierr, i, index_l, pos
        logical                         :: searching
        character(len=:), allocatable   :: usr_msg_1, usr_msg_2, word, string

        searching = .true.
        index_l   = 0
        string = self%get()

        ! Split the string in words and add each word to temp_list
        do while (searching)
            
            index_l = index_l + 1
            pos = index(string,separator)  

            if (pos == 0 ) then
                word = string
                searching = .false.
            else
                word = string(1:pos-1)
            end if
            
            ! Check if the number of words is greater than the size of temp_list
            usr_msg_1 = "string%strip: The number of words found in the stripped string "
            usr_msg_2 = "are too many."
            if (index_l > size(temp_list)) call chidg_signal_two(FATAL,usr_msg_1,self%str,usr_msg_2)
            
            ! Add word to temp_list
            temp_list(index_l) = word
            string = string(pos+1:)

        end do
        
        ! Allocate the vector of strings based on the word found
        allocate(string_list(index_l), stat=ierr)
        if (ierr /= 0) call chidg_signal(FATAL, "strip: error in allocating the string_list")

        ! Add words to the output vector
        do i = 1,index_l
            ! Adjust and trim the word and add it to string_list
            string_list(i)%str = trim(adjustl(temp_list(i)))
        end do

    end function strip
    !************************************************************************************************





    !>  Replace a charcter (before) with a new character (after) 
    !!  
    !!  @author Matteo Ugolotti
    !!  @date   2/12/2017
    !!
    !-----------------------------------------------------------------------------------------------
    function replace(self,before,after) result(string_out)
        class(string_t),    intent(in)  :: self
        character(*),       intent(in)  :: before
        character(*),       intent(in)  :: after
        
        type(string_t)                  :: string_out
        integer                         :: ierr, i, index_l, pos, after_len, before_len
        logical                         :: searching
        character(len=:), allocatable   :: word, string

        searching = .true.
        string = self%get()

        after_len  = len(after)
        before_len = len(before)
        i = 1

        ! Split the string in words and add each word to temp_list
        do while (searching)
            
            pos = index(string,before)  

            if (pos == 0 ) then
                word = string
                searching = .false.
            else
                word = string(:pos-1)// after // string(pos+before_len:)
            end if
            
            string = word

            if (i == 100) then
                searching = .false.
            end if
            
            i = i + 1

        end do

        ! Adjust and trim the word and add it to string_list
        call string_out%set(trim(adjustl(word)))

    end function replace
    !************************************************************************************************






end module mod_string
