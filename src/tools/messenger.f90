module messenger
    use mod_kinds,  only: rk,ik
    implicit none


contains

    !> Message routine for handling warnings and errors.
    !! Reports file name, line number, and warn/error message
    !!
    !! 'level' controls the action.
    !! - Warn            :: 1
    !! - Non-fatal error :: 2
    !! - Fatal error     :: 3
    subroutine message(signal,msg,filename,linenum)
        integer(ik),  intent(in)    :: signal
        character(*), intent(in)    :: msg
        character(*), intent(in)    :: filename
        integer(ik),  intent(in)    :: linenum

        character(100)   :: warnstr, errstr, killstr, genstr, starstr, linechar, dashstr


        warnstr = '*********************************  Warning  *********************************'
        errstr  = '*****************************  Non-fatal error  *****************************'
        killstr = '*******************************  Fatal error  *******************************'
        starstr = '*****************************************************************************'
        dashstr = '-----------------------------------------------------------------------------'

        write(linechar, '(i10)') linenum
        genstr = trim(filename) // ' at ' // adjustl(trim(linechar))


        ! Print message header
        print*, '| '//trim(dashstr)
        print*, '| '


        ! Select signal message
        select case (signal)
            case (0)    ! Normal message -- Code continues
                print*, '| '//trim(msg)

            case (1)    ! Warning message -- Code continues
                print*, '| '//trim(warnstr)
                print*, '| '//trim(genstr)
                print*, '| '
                print*, '| '//'Message:'
                print*, '| '//trim(msg)

            case (2)    ! Non-Fatal Error -- Code continues
                print*, '| '//trim(errstr)
                print*, '| '//trim(genstr)
                print*, '| '
                print*, '| '//'Message:'
                print*, '| '//trim(msg)

            case (3)    ! Fatal Error -- Code terminates
                print*, '| '//trim(killstr)
                print*, '| '//trim(genstr)
                print*, '| '
                print*, '| '//'Message:'
                print*, '| '//trim(msg)

            case default
                print*, "Messenger:message -- message code not recognized"
        end select


        ! Print message footer
        print*, '| '
        print*, '| '//trim(dashstr)








        ! Select signal action
        select case (signal)
            case (3)    ! Fatal Error -- Code terminates
                stop

            case default

        end select


    end subroutine


end module messenger
