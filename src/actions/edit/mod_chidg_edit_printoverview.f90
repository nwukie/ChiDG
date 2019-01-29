module mod_chidg_edit_printoverview
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use hdf5
    use h5lt

    use mod_hdf_utilities,  only: get_properties_hdf, get_ndomains_hdf, get_domain_names_hdf,   &
                                  get_domain_coordinate_orders_hdf, get_domain_field_orders_hdf,           &
                                  get_domain_equation_sets_hdf, get_contains_grid_hdf,          &
                                  get_contains_solution_hdf, get_domain_name_hdf

    implicit none










contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !-----------------------------------------------------------------------------------
    subroutine print_overview(fid,active_domain)
        integer(HID_T),     intent(in)              :: fid
        integer(HID_T),     intent(in), optional    :: active_domain


        integer(ik)                         :: ndom, idom, idom_hdf, ierr
        integer(ik),            allocatable :: dindices(:)
        character(len=1024),    allocatable :: dnames(:), eqnset(:)
        character(len=:),       allocatable :: dname_trim, active_domain_name
        integer(ik),            allocatable :: corder(:), sorder(:)
        logical                             :: contains_grid, contains_solution



        !
        ! Print header
        !
        call print_chidg_edit_header()



        !
        ! Get HDF5 ChiDG file information
        !
        ndom     = get_ndomains_hdf(fid)
        dnames   = get_domain_names_hdf(fid)
        eqnset   = get_domain_equation_sets_hdf(fid)



        !
        ! Check file contents
        !
        contains_grid     = get_contains_grid_hdf(fid)
        contains_solution = get_contains_solution_hdf(fid)

        !
        ! Handle contains_grid
        !
        if ( contains_grid ) then
            corder   = get_domain_coordinate_orders_hdf(fid,dnames)
        else
            allocate(corder(ndom), stat=ierr)
            if (ierr /= 0) call AllocationError

            corder = 0  ! 0-Order coordinate indicating no grid
        end if

        
        !
        ! Handle contains_solution
        !
        if ( contains_solution ) then
            sorder   = get_domain_field_orders_hdf(fid)
        else
            allocate(sorder(ndom), stat=ierr)
            if (ierr /= 0) call AllocationError

            sorder = 0  ! 0-Order solution inticating no solution
        end if




        
        
        !
        ! Print information
        !
        call write_line('Domain name', 'Coordinate expansion', 'Solution expansion', 'Equation set', delimiter='  :  ', columns=.True., column_width=20)

        !
        ! Find the correct hdf domain index to print so they are in order
        !
        do idom = 1,ndom

                dname_trim = trim(dnames(idom))

                if (present(active_domain)) then
                    active_domain_name = get_domain_name_hdf(active_domain)
                    if ( active_domain_name == dname_trim ) then
                        call write_line(dname_trim, corder(idom), sorder(idom), trim(eqnset(idom)), delimiter='  :  ', columns=.True., column_width=20, color='pink')
                    else
                        call write_line(dname_trim, corder(idom), sorder(idom), trim(eqnset(idom)), delimiter='  :  ', columns=.True., column_width=20)
                    end if

                else
                    call write_line(dname_trim, corder(idom), sorder(idom), trim(eqnset(idom)), delimiter='  :  ', columns=.True., column_width=20)
                end if


        end do !idom

        

        call write_line(" ")

    end subroutine print_overview
    !*****************************************************************************************************













    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------------
    subroutine print_chidg_edit_header()



        call write_line("______________________________________________________________________________________________________________________", color='blue', ltrim=.false.)
        call write_line("                                                                                                                      ", color='blue', ltrim=.false.)
        call write_line("                                                                                                                      ", color='blue', ltrim=.false.)
        call write_line("     _|_|_|  _|        _|  _|_|_|      _|_|_|                                                                         ", color='blue', ltrim=.false.)
        call write_line("   _|        _|_|_|        _|    _|  _|                                                                               ", color='blue', ltrim=.false.)
        call write_line("   _|        _|    _|  _|  _|    _|  _|  _|_|                               chidg edit                                ", color='blue', ltrim=.false.)
        call write_line("   _|        _|    _|  _|  _|    _|  _|    _|                                                                         ", color='blue', ltrim=.false.)
        call write_line("     _|_|_|  _|    _|  _|  _|_|_|      _|_|_|                                                                         ", color='blue', ltrim=.false.)
        call write_line("                                                                                                                      ", color='blue', ltrim=.false.)
        call write_line("______________________________________________________________________________________________________________________", color='blue', ltrim=.false.)









    end subroutine print_chidg_edit_header
    !*******************************************************************************************************






    !>
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------
    subroutine chidg_clear_screen()
        
        ! Calling execute_command_line with "clear"
        ! seemed to be causing an error on some systems.
        !call execute_command_line("clear")


        ! Clear screen by printing ANSI code
        print*, achar(27)//"[2J"

        ! Move cursor to top left
        print*, achar(27)//"[;H"


    end subroutine chidg_clear_screen
    !***************************************************************************










end module mod_chidg_edit_printoverview
