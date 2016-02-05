module mod_chidg_edit_printoverview
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use hdf5
    use h5lt

    use mod_hdf_utilities,  only: get_properties_hdf, get_ndomains_hdf, get_domain_names_hdf,       &
                                  get_order_coordinate_hdf, get_order_solution_hdf, get_eqnset_hdf, &
                                  get_domain_indices_hdf

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
        integer(ik),        intent(in), optional    :: active_domain


        integer(ik)                         :: ndom, idom, idom_hdf
        integer(ik),            allocatable :: dindices(:)
        character(len=1024),    allocatable :: dnames(:), eqnset(:)
        character(len=:),       allocatable :: dname_trim
        integer(ik),            allocatable :: corder(:), sorder(:)



        !
        ! Print header
        !
        call print_chidg_edit_header()



        !
        ! Get HDF5 ChiDG file information
        !
        ndom     = get_ndomains_hdf(fid)
        dnames   = get_domain_names_hdf(fid)
        dindices = get_domain_indices_hdf(fid)



        eqnset   = get_eqnset_hdf(fid,dnames)
        corder   = get_order_coordinate_hdf(fid,dnames)
        sorder   = get_order_coordinate_hdf(fid,dnames)



        
        !
        ! Print information
        !
        call write_line('Domain index', 'Domain name', 'Coordinate expansion', 'Solution expansion', 'Equation set', delimiter='  :  ', columns=.True., column_width=20)
        do idom_hdf = 1,ndom

            !
            ! Find the correct hdf domain index to print so they are in order
            !
            do idom = 1,ndom
                if ( (dindices(idom) == idom_hdf) ) then

                    !
                    ! Trim domain identifier. Use dname_trim(3:) to remove 'D_'
                    !
                    dname_trim = trim(dnames(idom))

                    if (present(active_domain)) then
                        if ( active_domain == idom_hdf ) then
                            call write_line(idom_hdf,dname_trim(3:), corder(idom), sorder(idom), trim(eqnset(idom)), delimiter='  :  ', columns=.True., column_width=20, color='pink')
                        else
                            call write_line(idom_hdf,dname_trim(3:), corder(idom), sorder(idom), trim(eqnset(idom)), delimiter='  :  ', columns=.True., column_width=20)
                        end if

                    else
                        call write_line(idom_hdf,dname_trim(3:), corder(idom), sorder(idom), trim(eqnset(idom)), delimiter='  :  ', columns=.True., column_width=20)
                    end if

                end if ! dindices

            end do !idom

        end do  ! idom_hdf
        

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


        call write_line(" ____________________________________________________________________________________________________________________", color='blue', ltrim=.false.)
        call write_line("        ______               ____     _____                                                                          ", color='blue', ltrim=.false.)
        call write_line("       |        |           |    \   |     |                                                                         ", color='blue', ltrim=.false.)
        call write_line("       |        |        .  |     |  |                                                                               ", color='blue', ltrim=.false.)
        call write_line("       |        |_____      |     |  |   __                              chidg edit                                  ", color='blue', ltrim=.false.)
        call write_line("       |        |     |  |  |     |  |     |                                                                         ", color='blue', ltrim=.false.)
        call write_line("       |        |     |  |  |     |  |     |                                                                         ", color='blue', ltrim=.false.)
        call write_line("       |______  |     |  |  |____/   |_____|                                                                         ", color='blue', ltrim=.false.)
        call write_line(" ____________________________________________________________________________________________________________________", color='blue', ltrim=.false.)
        call write_line(" ")


    end subroutine print_chidg_edit_header
    !*******************************************************************************************************














end module mod_chidg_edit_printoverview
