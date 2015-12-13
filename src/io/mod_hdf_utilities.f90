module mod_hdf_utilities
#include <messenger.h>
    use mod_kinds,              only: rk, ik

    use type_file_properties,   only: file_properties_t

    use hdf5
    use h5lt
    implicit none







contains







    !> Return a properties instance containing information about an hdf file
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!  @param[in]  filename    Character string containing a filename for a file that gets interrogated
    !!
    !!  @result     prop        file_properties_t instance that gets returned with file information
    !!
    !----------------------------------------------------------------------------------------------
    function get_properties_hdf(filename) result(prop)
        character(*),   intent(in)  :: filename


        integer(HID_T)              :: fid
        integer(ik)                 :: ierr
        integer(ik), allocatable    :: nterms_1d(:)
        logical                     :: fileexists = .false.

        type(file_properties_t)     :: prop

        ! 
        !  Check file exists
        !
        inquire(file=filename, exist=fileexists)
        if (.not. fileexists) then
            call chidg_signal(FATAL,'get_properties_hdf: Could not find grid file')
        end if


        !
        !  Initialize Fortran interface.
        !
        call h5open_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'get_properties_hdf - h5open_f: HDF5 Fortran interface had an error during initialization')



        !
        !  Open input file using default properties.
        !
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'get_properties_hdf - h5fopen_f: There was an error opening the grid file.')



        !
        ! Get number of domains
        !
        prop%ndomains = get_ndomains_hdf(fid)


        !
        ! Check for Grid and Solution contents
        !
        prop%contains_grid      = check_contains_grid_hdf(fid)
        prop%contains_solution  = check_contains_solution_hdf(fid)


        !
        ! Get domain names
        !
        prop%domain_names = get_domain_names_hdf(fid)




        !
        ! Get order of coordinate and solution expansions
        !
        if ( prop%contains_grid ) then
            prop%order_c = get_order_coordinate_hdf(fid,prop%domain_names)
        end if

        if ( prop%contains_solution ) then
            prop%order_s = get_order_solution_hdf(fid,prop%domain_names)
        end if


        !
        ! Compute number of terms in the polynomial expansions
        !
        nterms_1d = (prop%order_c + 1)
        prop%nterms_c = nterms_1d * nterms_1d * nterms_1d

        nterms_1d = (prop%order_s + 1)
        prop%nterms_s = nterms_1d * nterms_1d * nterms_1d



        !
        ! Get equation set for each domain
        !
        prop%eqnset = get_eqnset_hdf(fid, prop%domain_names)




        !
        ! Close file
        !
        call h5fclose_f(fid,ierr)

    end function get_properties_hdf
    !###############################################################################################








    !> Given a file identifier, return the number of domains in an hdf5 file.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------
    function get_ndomains_hdf(fid) result(ndomains)
        integer(HID_T), intent(in)  :: fid


        integer(ik)             :: ndomains, ierr
        integer, dimension(1)   :: buf

        !
        !  Get number of domains from attribute 'ndomains' in file root
        !
        call h5ltget_attribute_int_f(fid, "/", 'ndomains', buf, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'get_ndomains_hdf: h5ltget_attribute_int_f had a problem getting the number of domains')

        ndomains = int(buf(1), kind=ik)


    end function get_ndomains_hdf
    !###############################################################################################









    !> Test if the file contains a Grid.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    function check_contains_grid_hdf(fid) result(grid_status)
        integer(HID_T), intent(in)  :: fid

        logical             :: grid_status
        character(len=10)   :: contains_grid_attr
        integer             :: ierr


        !
        ! Get attribute for 'contains_grid
        !
        call h5ltget_attribute_string_f(fid, "/", 'contains_grid', contains_grid_attr, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"check_contains_grid_hdf - h5ltget_attribute_int_f")


        !
        ! Test grid attribute
        !
        if (trim(contains_grid_attr) == 'Yes') then
            grid_status = .true.
        else
            grid_status = .false.
        end if



    end function check_contains_grid_hdf
    !##################################################################################################








    !> Test if the file contains a Solution.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    function check_contains_solution_hdf(fid) result(solution_status)
        integer(HID_T), intent(in)  :: fid

        logical             :: solution_status
        character(len=10)   :: contains_solution_attr
        integer             :: ierr


        !
        ! Get attribute for 'contains_grid
        !
        call h5ltget_attribute_string_f(fid, "/", 'contains_grid', contains_solution_attr, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"check_contains_solution - h5ltget_attribute_string_f")


        !
        ! Test solution attribute
        !
        if (trim(contains_solution_attr) == 'Yes') then
            solution_status = .true.
        else
            solution_status = .false.
        end if



    end function check_contains_solution_hdf
    !##################################################################################################








    !> Return a list of domain names from an HDF5 file identifier.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    function get_domain_names_hdf(fid) result(names)
        integer(HID_T),     intent(in)  :: fid

        character(len=1024), allocatable    :: names(:)
        character(len=1024)                 :: gname
        integer                             :: ndomains, nmembers, type
        integer                             :: igrp, idom, ierr

        !
        ! Get number of domains
        !
        ndomains = get_ndomains_hdf(fid)


        !
        ! Allocate names
        !
        allocate(names(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        !  Get number of groups in the file root
        !
        call h5gn_members_f(fid, "/", nmembers, ierr)


        !
        !  Loop through groups and read domain names
        !
        idom = 1
        do igrp = 0,nmembers-1
            !
            ! Get group name
            !
            call h5gget_obj_info_idx_f(fid,"/", igrp, gname, type, ierr)

            !
            ! Test if group is a 'Domain'
            !
            if (gname(1:2) == 'D_') then

                !
                ! Store name
                !
                names(idom) = trim(gname)

                idom = idom + 1
            end if
        end do


    end function get_domain_names_hdf
    !##################################################################################################






    !> Returns an array of integer that specifies the order of the coordinate expansion for every domain.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    function get_order_coordinate_hdf(fid, dnames) result(orders)
        integer(HID_T),         intent(in)  :: fid
        character(len=1024),    intent(in)  :: dnames(:)

        integer(ik), allocatable    :: orders(:)
        integer                     :: ierr, idom, mapping
        integer, dimension(1)       :: buf


        !
        ! Allocate storage for orders
        !
        allocate(orders(size(dnames)), stat=ierr)
        if (ierr /= 0) call AllocationError

        !
        !  Loop through groups and read domains
        !
        do idom = 1,size(dnames)

            !
            !  Get coordinate mapping
            !
            call h5ltget_attribute_int_f(fid, trim(dnames(idom)), 'mapping', buf, ierr)
            if (ierr /= 0) stop "Error: get_order_coordinate_hdf - h5ltget_attribute_int_f"


            !
            ! Compute number of terms in coordinate expansion
            !
            mapping = buf(1)

            orders(idom) = int(mapping, kind=ik)

        end do

    end function get_order_coordinate_hdf
    !########################################################################################################







    !> Returns an array of integer that specifies the order of the solution expansion for every domain.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    function get_order_solution_hdf(fid, dnames) result(orders)
        integer(HID_T),         intent(in)  :: fid
        character(len=1024),    intent(in)  :: dnames(:)

        integer(ik), allocatable    :: orders(:)
        integer                     :: ierr, idom, order
        integer, dimension(1)       :: buf


        !
        ! Allocate storage for orders
        !
        allocate(orders(size(dnames)), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        !  Loop through groups and read domains
        !
        do idom = 1,size(dnames)

            !
            !  Get coordinate mapping
            !
            call h5ltget_attribute_int_f(fid, trim(dnames(idom)), 'order_solution', buf, ierr)
            if (ierr /= 0) stop "Error: get_order_coordinate_hdf - h5ltget_attribute_int_f"


            !
            ! Compute number of terms in coordinate expansion
            !
            order = buf(1)

            orders(idom) = int(order, kind=ik)

        end do

    end function get_order_solution_hdf
    !########################################################################################################










    !> Returns an array of integer that specifies the order of the solution expansion for every domain.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    function get_eqnset_hdf(fid, dnames) result(eqnsets)
        integer(HID_T),         intent(in)  :: fid
        character(len=1024),    intent(in)  :: dnames(:)

        character(len=1024), allocatable    :: eqnsets(:)
        integer                             :: ierr, idom


        !
        ! Allocate storage for eqnsets
        !
        allocate(eqnsets(size(dnames)), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        !  Loop through groups and read domains
        !
        do idom = 1,size(dnames)

            !
            !  Get eqnset string
            !
            call h5ltget_attribute_string_f(fid, trim(dnames(idom)), 'eqnset', eqnsets(idom), ierr)
            if (ierr /= 0) stop "Error: get_eqnset_hdf - h5ltget_attribute_int_f"


        end do

    end function get_eqnset_hdf
    !########################################################################################################










end module mod_hdf_utilities
