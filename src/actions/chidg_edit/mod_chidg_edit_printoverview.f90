module mod_chidg_edit_printoverview
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use hdf5
    use h5lt

    use mod_hdf_utilities,  only: get_properties_hdf, get_ndomains_hdf, get_domain_names_hdf, &
                                  get_order_coordinate_hdf, get_order_solution_hdf, get_eqnset_hdf

    implicit none










contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !-----------------------------------------------------------------------------------
    subroutine print_overview(fid)
        integer(HID_T),     intent(in)  :: fid


        integer(ik)                         :: ndom, idom
        character(len=1024),    allocatable :: dnames(:), eqnset(:)
        integer(ik),            allocatable :: corder(:), sorder(:)


        !
        ! Get HDF5 ChiDG file information
        !
        ndom   = get_ndomains_hdf(fid)
        dnames = get_domain_names_hdf(fid)
        eqnset = get_eqnset_hdf(fid,dnames)
        corder = get_order_coordinate_hdf(fid,dnames)
        sorder = get_order_coordinate_hdf(fid,dnames)



        
        !
        ! Print information
        !
        call write_line('Domain index', 'Domain name', 'Coordinate expansion', 'Solution expansion', 'Equation set', delimiter='  :  ', columns=.True., column_width=20)
        do idom = 1,ndom
            call write_line(idom,trim(dnames(idom)), corder(idom), sorder(idom), trim(eqnset(idom)), delimiter='  :  ', columns=.True., column_width=20)
        end do
        

    end subroutine print_overview
    !************************************************************************************







end module mod_chidg_edit_printoverview
