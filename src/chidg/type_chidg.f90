module type_chidg
#include <messenger.h>
    use mod_equations,  only: initialize_equations
    use mod_grid,       only: initialize_grid
    use mod_io,         only: read_input

    implicit none


    type, public    :: chidg_t


    contains
        procedure   :: init
    end type



contains

    !> ChiDG environment initialization routine
    !!      - Call initiailization procedures for equations, grid data, reading input
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  level   Initialization level specification. 'env' or 'full'
    !---------------------------------------------------------------
    subroutine init(self,level)
        class(chidg_t),  intent(inout)   :: self
        character(*),    intent(in)      :: level

        logical :: valid


        !
        ! Check for a valid 'level' that designates which initialization calls get executed
        ! Valid strings are:
        !   - 'env'     Basic environment initialization. Equations and supporting grid data
        !   - 'full'    Run 'env' routines plus call read_input to pull data from chidg.nml file 
        !
        if ( (trim(level) == 'env') .or.  & 
             (trim(level) == 'full') ) then
             valid = .true.
        else
             valid = .false.
        end if
        if (.not. valid) call signal(FATAL,'chidg%init: Initialization level is not valid')




        !
        ! Run these commands for all initialization levels. 'env'
        !
        call initialize_equations()
        call initialize_grid()




        !
        ! Run these commands for initialization level 'full'
        !
        if (trim(level) == 'full') then
            call read_input()
        end if



    end subroutine






end module type_chidg
