module messenger
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: IO_DESTINATION
    implicit none


    character(len=:), allocatable   :: line                         ! Line that gets assembled and written
    character(len=2), parameter     :: default_delimiter = '  '     ! Delimiter of line parameters
    character(len=:), allocatable   :: current_delimiter            ! Delimiter of line parameters
    integer                         :: unit                         ! Unit of log file
    integer                         :: max_msg_length = 300         ! Maximum width of message line
    logical                         :: log_initialized = .false.    ! Status of log file

    character(len=:), allocatable   :: color_begin
    character(len=:), allocatable   :: color_end

contains



    !> Log initialization
    !!
    !!  Gets new available file unit and opens log file. 'unit' is a module 
    !!  variable that can be used throughout the module to access the log file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!------------------------------------------------------------------------------------------------------------
    subroutine log_init()

        logical :: file_opened = .false.

        !
        ! Open file
        !
        inquire(file='chidg.log', opened=file_opened)

        if ( .not. file_opened ) then
            open(newunit=unit, file='chidg.log')
        end if

        log_initialized = .true.

    end subroutine log_init
    !*************************************************************************************************************






    !> Log finalization
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!------------------------------------------------------------------------------------------------------------
    subroutine log_finalize()

        logical :: file_opened = .false.

        !
        ! Close file
        !
        inquire(file='chidg.log', opened=file_opened)

        if ( file_opened ) then
            close(unit)
        end if

    end subroutine log_finalize
    !*************************************************************************************************************







    !> Message routine for handling warnings and errors. Reports file name, line number,
    !! and warn/error message. This would usually not be called directly. Rather, use
    !! the macro defined in message.h that automatically inserts filename and linenumber.
    !!
    !! 'level' controls the action.
    !! - Warn            :: 1
    !! - Non-fatal error :: 2
    !! - Fatal error     :: 3
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!  @param[in]  pathname    Path and name of the file that the message is coming from.
    !!  @param[in]  linenum     Line number in the file that 'message' was called from.
    !!  @param[in]  sig         Signal level of the message. Defined above.
    !!  @param[in]  msg         Accompanying message to print.
    !!  @param[in]  info_one    Optional auxiliary information to be reported.
    !!  @param[in]  info_two    Optional auxiliary information to be reported.
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine message(pathname, linenum, sig, msg, info_one, info_two, info_three)
        character(*), intent(in)                        :: pathname
        integer(ik),  intent(in)                        :: linenum
        integer(ik),  intent(in)                        :: sig
        character(*), intent(in)                        :: msg
        class(*),     intent(in), target,   optional    :: info_one
        class(*),     intent(in), target,   optional    :: info_two
        class(*),     intent(in), target,   optional    :: info_three

        integer                         :: iaux, pathstart
        integer(ik)                     :: ierr
        character(len=:), allocatable   :: subpath, temppath
        class(*), pointer               :: auxdata => null()
        character(100)                  :: warnstr, errstr, killstr, genstr, starstr, linechar, dashstr, blankstr
        logical                         :: print_info_one   = .false.
        logical                         :: print_info_two   = .false.
        logical                         :: print_info_three = .false.


        warnstr =  '***************************************  Warning  ***************************************'
        errstr  =  '***********************************  Non-fatal error  ***********************************'
        killstr =  '*************************************  Fatal error  *************************************'
        starstr =  '*****************************************************************************************'
        dashstr =  '-----------------------------------------------------------------------------------------'
        blankstr = '               '



        !
        ! Chop off unimportant beginning of file path
        !
        temppath = pathname
        pathstart = index(temppath, '/ChiDG/')
        if (pathstart == 0) then
            subpath = temppath      ! The intel compiler provides just the file name, without the path. So here, we just take the file name and don't chop anything
        else
            subpath = temppath(pathstart:len(pathname))
        end if


        !
        ! Assemble string including file name and line number
        !
        write(linechar, '(i10)') linenum
        genstr = trim(subpath) // ' at ' // adjustl(trim(linechar))


        !
        ! Print message header
        !
        call write_line(blankstr)
        call write_line(trim(dashstr))




        !
        ! Select signal message
        !
        select case (sig)
            case (0)    ! Normal message -- Code continues
                call write_line(trim(msg))

            case (1)    ! Warning message -- Code continues
                call write_line(trim(warnstr))
                call write_line(trim(genstr))
                call write_line(blankstr)
                call write_line('Message:')
                call write_line(trim(msg))

            case (2)    ! Non-Fatal Error -- Code continues
                call write_line(trim(errstr))
                call write_line(trim(genstr))
                call write_line(blankstr)
                call write_line('Message:')
                call write_line(trim(msg))

            case (3)    ! Fatal Error -- Code terminates
                
                call write_line(trim(killstr))
                call write_line(trim(genstr))
                call write_line(blankstr)
                call write_line('Message:')
                call write_line(trim(msg))

            case default
                print*, "Messenger:message -- message code not recognized"
                stop
        end select
        call write_line(blankstr)



        !
        ! Loop through auxiliary variables. If present, try to print.
        !
        do iaux = 1,3

            print_info_one   = ( present(info_one)   .and. (iaux == 1) )
            print_info_two   = ( present(info_two)   .and. (iaux == 2) )
            print_info_three = ( present(info_three) .and. (iaux == 3) )


            if ( print_info_one )   auxdata => info_one
            if ( print_info_two )   auxdata => info_two
            if ( print_info_three ) auxdata => info_three


            !
            ! auxdata pointer is used to point to current auxiliary data variable and then go through the available IO types
            !
            if ( associated(auxdata) ) then
                call write_line('Case specific info:')

                select type(auxdata)
                    type is(integer)
                        call write_line(auxdata)

                    type is(integer(8))
                        call write_line(auxdata)

                    type is(real)
                        call write_line(auxdata)

                    type is(real(8))
                        call write_line(auxdata)

                    type is(character(*))
                        call write_line(auxdata)

                    class default
                        print*, '', "Data type not implemented for I/O in messege.f90"
                end select

            end if ! present(info_one)



            !
            ! Disassociate pointer so it doesn't try to print the same thing twice in some cases.
            !
            auxdata => null()

        end do ! iaux



!        !
!        ! If msg procedure is present, call that.
!        !
!        if ( present(msg_proc) ) then
!            !call msg_proc()
!            ierr = msg_proc()
!        end if
!



        !
        ! Print message footer
        !
        call write_line(blankstr)
        call write_line(dashstr)



        !
        ! Select signal action
        !
        select case (sig)
            case (3)    ! Fatal Error -- Code terminates
                stop
                !call MPI_Abort(MPI_COMM_WORLD,sig,ierr)

            case default

        end select



    end subroutine message
    !**********************************************************************************************************







    !> This subroutine writes a line to IO that is composed of 8 optional incoming variables.
    !! This is accomplished by first passing each component to the 'add_to_line' subroutine, which
    !! assembles the data into the 'line' module-global variable. Then, the 'send_line' subroutine is called
    !! to handle the destination of the line to either the screen, a file, or both.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!  @param[in]  a-h             Optional polymorphic variables that can be used to compose a line sent to IO.
    !!  @param[in]  delimiter       Optional delimiter that is used to separate line components. Default is set to ' '
    !!  @param[in]  columns         Logical optional to indicate if incoming arguments should be aligned in columns.
    !!  @param[in]  column_width    Optional integer indicating the column width if columns was indicated.
    !!
    !----------------------------------------------------------------------------------------------------------
    subroutine write_line(a,b,c,d,e,f,g,h,delimiter,columns,column_width,color,ltrim)
        class(*),           intent(in), target, optional        :: a
        class(*),           intent(in), target, optional        :: b
        class(*),           intent(in), target, optional        :: c
        class(*),           intent(in), target, optional        :: d
        class(*),           intent(in), target, optional        :: e
        class(*),           intent(in), target, optional        :: f
        class(*),           intent(in), target, optional        :: g
        class(*),           intent(in), target, optional        :: h
        character(len=*),   intent(in),         optional        :: delimiter
        logical,            intent(in),         optional        :: columns
        integer(ik),        intent(in),         optional        :: column_width
        character(*),       intent(in),         optional        :: color
        logical,            intent(in),         optional        :: ltrim
        

        class(*), pointer               :: auxdata => null()

        integer :: iaux
        logical :: print_info_one, print_info_two, print_info_three, print_info_four, print_info_five
        logical :: print_info_six, print_info_seven, print_info_eight

        
        !
        ! Loop through variables and compose line to write
        !
        do iaux = 1,8

            print_info_one   = ( present(a) .and. (iaux == 1) )
            print_info_two   = ( present(b) .and. (iaux == 2) )
            print_info_three = ( present(c) .and. (iaux == 3) )
            print_info_four  = ( present(d) .and. (iaux == 4) )
            print_info_five  = ( present(e) .and. (iaux == 5) )
            print_info_six   = ( present(f) .and. (iaux == 6) )
            print_info_seven = ( present(g) .and. (iaux == 7) )
            print_info_eight = ( present(h) .and. (iaux == 8) )

            if ( print_info_one   )  auxdata => a
            if ( print_info_two   )  auxdata => b
            if ( print_info_three )  auxdata => c
            if ( print_info_four  )  auxdata => d
            if ( print_info_five  )  auxdata => e
            if ( print_info_six   )  auxdata => f
            if ( print_info_seven )  auxdata => g
            if ( print_info_eight )  auxdata => h



            if ( associated(auxdata) ) then

                    !
                    ! Add data to line
                    !
                    call add_to_line(auxdata,delimiter,columns,column_width,color,ltrim)

            end if



            !
            ! Unassociate pointer
            !
            auxdata => null()


        end do




        !
        ! Send line to output
        !
        call send_line()



    end subroutine write_line
    !************************************************************************************************************








    !> Adds data to the module-global 'line' character string
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!  @param[in]  linedata        Polymorphic data component that gets converted to a string and added to the IO line
    !!  @param[in]  delimiter       Optional delimiter that is used to separate line components. Default is set to ' '
    !!  @param[in]  columns         Logical optional to indicate if incoming arguments should be aligned in columns.
    !!  @param[in]  column_width    Optional integer indicating the column width if columns was indicated.
    !!
    !--------------------------------------------------------------------------------------------------------------
    subroutine add_to_line(linedata,delimiter,columns,column_width,color,ltrim)
        class(*),       intent(in)              :: linedata
        character(*),   intent(in), optional    :: delimiter
        logical,        intent(in), optional    :: columns
        integer(ik),    intent(in), optional    :: column_width
        character(*),   intent(in), optional    :: color
        logical,        intent(in), optional    :: ltrim

        character(100)                  :: write_internal
        character(len=:),   allocatable :: temp, temp_a, temp_b
        integer(ik)                     :: current_width, ispace, extra_space, test_blank
        logical                         :: blank_line


        !
        ! Initialize temp
        !
        temp = ''


        !
        ! Select delimiter
        !
        if (  present(delimiter) ) then
            current_delimiter = delimiter
        else
            current_delimiter = default_delimiter
        end if


        !
        ! Set color
        !
        if ( present(color) ) then
            call set_color(color)
        else
            call set_color('none')
        end if




        !
        ! Add to line. Since variable is polymorphic, we have to test for each type and handle
        ! appropriately. Numeric data gets first written to a string variable and then concatatenated to 
        ! the module-global 'line' variable.
        !
        select type(linedata)

            type is(character(len=*))
                temp = linedata

            type is(integer)
                if ( linedata == 0 ) then
                    write_internal = '0'
                else
                    write(write_internal, '(I10.0)') linedata
                end if
                temp = write_internal
            type is(integer(8))
                if ( linedata == 0 ) then
                    write_internal = '0'
                else
                    write(write_internal, '(I10.0)') linedata
                end if
                temp = write_internal

            type is(real)
                if (linedata > 0.1) then
                    write(write_internal, '(F24.14)') linedata
                else
                    write(write_internal, '(E24.14)') linedata
                end if
                temp = write_internal

            type is(real(8))
                if (linedata > 0.1) then
                    write(write_internal, '(F24.14)') linedata
                else
                    write(write_internal, '(E24.14)') linedata
                end if
                temp = write_internal

            class default
                print*, 'Error: no IO rule for provided data in add_to_line'
                stop
        end select


        !
        ! Test for blank line
        !
        test_blank = verify(temp, set=' ')
        if ( test_blank == 0 ) then
            blank_line = .true.
        else
            blank_line = .false.
        end if


        !
        ! Rules for neatening up the string. Check blank string. Check ltrim.
        !
        if ( blank_line ) then
            temp_a = temp   ! blank line, done do any modification.

        else if ( present(ltrim) ) then
            if ( ltrim ) then
                temp_a = trim(adjustl(temp))    ! trim line if explicitly requested.
            else 
                temp_a = temp                   ! explicitly requested to not trim line.
            end if

        else
            temp_a = trim(adjustl(temp))        ! default, trim the line.
        end if





        !
        ! Append spaces for column alignment
        !
        if ( present(columns) ) then
            if (columns) then
                current_width = len_trim(temp_a)  ! length without trailing blanks
                extra_space = column_width - current_width

                !
                ! Add spaces on front and back. Could cause misalignment with dividing integer by 2.
                !
                do while ( len(temp_a) < column_width )
                    temp_a = ' '//temp_a//' '
                end do

                !
                ! Make sure we are at exactly column_width. Could have gone over the the step above.
                !
                temp_b = temp_a(1:column_width)
                 
            end if
        else
            
            temp_b = temp_a

        end if



        ! Append new text to line
        !line = line//temp_b
        line = line//color_begin//temp_b//color_end

        ! Append delimiter
        line = line//current_delimiter


    end subroutine add_to_line
    !********************************************************************************************************









    !>  Handles sending module-global 'line' string to a destination of either the screen, a file, or both.
    !!  Is is determined by the IO_DESTINATION variable from mod_constants.
    !!
    !!  Line wrapping is also handled, set by the module parameter 'max_msg_length'.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine send_line()
        integer :: delimiter_size
        integer :: line_size
        integer :: line_trim
        integer :: lstart, lend, section

        character(len=:), allocatable :: writeline
        integer                       :: section_length

        !
        ! Get line section length
        !
        if ( len(line) > max_msg_length ) then
            section_length = max_msg_length
        else
            section_length = len(line)
        end if


        !
        ! Get line/delimiter sizes
        !
        delimiter_size = len(current_delimiter)
        line_size      = len(line)
        line_trim      = line_size - delimiter_size



        !
        ! Remove trailing delimiter
        !
        line = line(1:line_trim)



        !
        ! Handle line IO. Writes the line in chunks for line-wrapping until the entire line has been processed.
        !
        writeline = line
        section = 1
        lend    = -1
        do while ( lend /= len(line) ) 
            lstart = 1 + ((section-1) * section_length)
            lend   = (section) * section_length

            !
            ! Dont go out-of-bounds
            !
            if (lstart > len(line) ) then
                exit
            end if

            if (lend > len(line)) then
                lend = len(line)
            end if


            !
            ! Get line for writing
            !
            writeline = line(lstart:lend)



            !
            ! Colorize output
            !
            !writeline = color_begin//writeline//color_end






            !
            ! Write to destination
            !
            if ( trim(IO_DESTINATION) == 'screen' ) then
                print*, writeline


            else if ( trim(IO_DESTINATION) == 'file' ) then
                if (log_initialized) then
                    write(unit,*) writeline
                else
                    stop "Trying to write a line, but log file not inititlized. Call chidg%init('env')"
                end if



            else if ( trim(IO_DESTINATION) == 'both' ) then
                if (log_initialized) then
                    print*, writeline
                    write(unit,*) writeline
                else
                    stop "Trying to write a line, but log file not inititlized. Call chidg%init('env')"
                end if



            else
                print*, "Error: value for IO_DESTINATION is invalid. Valid selections are 'screen', 'file', 'both'."

            end if


            !
            ! Next line chunk to write
            !
            section = section + 1

        end do ! len(line) > max_msg_length


        !
        ! Clear line
        !
        line = ''



    end subroutine send_line
    !**************************************************************************************************************











    !>  Set line color.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  color   String indicating the color to set
    !!
    !-------------------------------------------------------------------------------------------
    subroutine set_color(color)
        character(*),   intent(in)  :: color

        color_end = achar(27)//'[m'

        select case (color)
            case ('black')
                color_begin = achar(27)//'[30m'

            case ('red')
                color_begin = achar(27)//'[31m'

            case ('green')
                color_begin = achar(27)//'[32m'

            case ('yellow')
                color_begin = achar(27)//'[33m'

            case ('blue')
                color_begin = achar(27)//'[34m'

            case ('purple')
                color_begin = achar(27)//'[35m'

            case ('aqua')
                color_begin = achar(27)//'[36m'

            case ('pink')
                color_begin = achar(27)//'[95m'

            case ('none')
                color_begin = ''
                color_end   = ''

            case default
                color_begin = achar(27)//'[30m'
                call message(__FILE__,__LINE__,1, "set_color: unrecognized color string.",color) ! send warning

        end select

    end subroutine set_color
    !******************************************************************************************












end module messenger
