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
    subroutine warn(level,filename,linenum,messg)
        integer(ik),  intent(in)    :: level
        character(*), intent(in)    :: filename
        integer(ik),  intent(in)    :: linenum
        character(*), intent(in)    :: messg

        character(100)   :: warnstr
        character(100)   :: errstr
        character(100)   :: killstr
        character(100)  :: genstr
        character(100)  :: starstr
        character(100)   :: linechar

        warnstr = '**************  Warning  **************'
        errstr  = '**********  Non-fatal error  **********'
        killstr = '************  Fatal error  ************'
        starstr = '***************************************'

        write(linechar, '(i10)') linenum
        genstr = trim(filename) // ' at ' // adjustl(trim(linechar))



        ! Warning message - code continues
        if (level == 1) then
            print*, trim(warnstr)
            print*, trim(genstr)
            print*, trim(messg)
            print*, starstr


        ! Error message - code continues
        else if (level == 2) then
            print*, trim(errstr)
            print*, trim(genstr)
            print*, trim(messg)
            print*, starstr

        ! Kill message - stop code
        else if (level == 3) then
!            print*, errstr, filename, trim(' at '), linenum, messg
            print*, trim(killstr)
            print*, trim(genstr)
            print*, trim(messg)
            print*, starstr
            stop

        end if
    end subroutine


end module messenger
