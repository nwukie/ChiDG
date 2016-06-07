module mod_hdf_utilities
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: NFACES, TWO_DIM, THREE_DIM
    use type_file_properties,   only: file_properties_t
    use hdf5
    use h5lt
    implicit none





contains





    !> Return a properties instance containing information about an hdf file
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!
    !!  @param[in]  filename    Character string containing a filename for a file that gets interrogated
    !!  @result     prop        file_properties_t instance that gets returned with file information
    !!
    !-----------------------------------------------------------------------------------------------------------
    function get_properties_hdf(filename) result(prop)
        character(*),   intent(in)  :: filename

        integer(HID_T)              :: fid
        integer(ik)                 :: ierr, idom
        integer(ik)                 :: nterms_1d
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
        call prop%set_ndomains(prop%ndomains)


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
        ! Get number of spatial dimensions
        !
        prop%spacedim = get_spacedim_hdf(fid,prop%domain_names)




        !
        ! Compute number of terms in the polynomial expansions for each domain
        !
        do idom = 1,prop%ndomains
            

            nterms_1d = (prop%order_c(idom) + 1)
            if ( prop%spacedim(idom) == THREE_DIM ) then
                prop%nterms_c(idom) = nterms_1d * nterms_1d * nterms_1d
            else if ( prop%spacedim(idom) == TWO_DIM ) then
                prop%nterms_c(idom) = nterms_1d * nterms_1d
            end if


 
            nterms_1d = (prop%order_s(idom) + 1)
            if ( prop%spacedim(idom) == THREE_DIM ) then
                prop%nterms_s(idom) = nterms_1d * nterms_1d * nterms_1d
            else if ( prop%spacedim(idom) == TWO_DIM ) then
                prop%nterms_s(idom) = nterms_1d * nterms_1d
            end if


        end do ! idom


        !
        ! Get equation set for each domain
        !
        prop%eqnset = get_eqnset_hdf(fid, prop%domain_names)




        !
        ! Close file
        !
        call h5fclose_f(fid,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_properties_hdf: h5fclose.")        
        call h5close_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_properties_hdf: h5close.")        

    end function get_properties_hdf
    !**************************************************************************************************************












    !>  Given a file identifier, return the number of domains in an hdf5 file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier
    !!
    !--------------------------------------------------------------------------------------------------------------
    function get_ndomains_hdf(fid) result(ndomains)
        integer(HID_T), intent(in)  :: fid


        integer                 :: ierr
        !integer(ik)             :: ndomains, ierr
        integer(ik)             :: ndomains
        integer, dimension(1)   :: buf


        !
        !  Get number of domains from attribute 'ndomains' in file root
        !
        call h5ltget_attribute_int_f(fid, "/", "ndomains", buf, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_ndomains_hdf: h5ltget_attribute_int_f had a problem getting the number of domains")

        ndomains = int(buf(1), kind=ik)


    end function get_ndomains_hdf
    !**************************************************************************************************************












    !>  Test if the file contains a Grid.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid             HDF5 file identifier.
    !!  @result     grid_status     Logical indicating if ChiDG grid exists in file, fid.
    !!
    !--------------------------------------------------------------------------------------------------------------
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
    !**************************************************************************************************************














    !>  Test if the file contains a Solution.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid                 HDF5 file identifier.
    !!  @result     solution_status     Logical indicating if ChiDG solution exists in file, fid.
    !!
    !---------------------------------------------------------------------------------------------------------------
    function check_contains_solution_hdf(fid) result(solution_status)
        integer(HID_T), intent(in)  :: fid

        logical             :: solution_status
        character(len=10)   :: contains_solution_attr
        integer             :: ierr


        !
        ! Get attribute for 'contains_grid
        !
        call h5ltget_attribute_string_f(fid, "/", 'contains_solution', contains_solution_attr, ierr)
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
    !***************************************************************************************************************












    !>  Return a list of domain names from an HDF5 file identifier.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier
    !!
    !---------------------------------------------------------------------------------------------------------------
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
    !*************************************************************************************************************














    !>  Return a domain name given a domain index.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid         HDF file identifier
    !!  @param[in]  idom_hdf    A specified domain index to be queried. This is an attribute of each domain in 
    !!                          the HDF file, per the ChiDG convention.
    !!
    !---------------------------------------------------------------------------------------------------------------
    function get_domain_name_hdf(fid,idom_hdf) result(dname)
        integer(HID_T),     intent(in)  :: fid
        integer(ik),        intent(in)  :: idom_hdf

        character(len=1024)                 :: dname
        character(len=1024), allocatable    :: dnames(:)
        integer(ik),         allocatable    :: dindices(:)
        integer                             :: iind, ndomains

        !
        ! Get number of domains, domain names, and domain indices
        !
        ndomains = get_ndomains_hdf(fid)
        dnames   = get_domain_names_hdf(fid)
        dindices = get_domain_indices_hdf(fid)



        do iind = 1,ndomains

            if ( dindices(iind) == idom_hdf ) then
               dname = dnames(iind) 
            end if

        end do



    end function get_domain_name_hdf
    !*************************************************************************************************************

















    !> Return a list of domain indices from an HDF5 file identifier. This is because, the current method of detecting
    !! domains by name can change the order they are detected in. So, each domain is given an idomain attribute that 
    !! is independent of the order of discovery from the file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier
    !!
    !---------------------------------------------------------------------------------------------------------------
    function get_domain_indices_hdf(fid) result(indices)
        integer(HID_T),     intent(in)  :: fid

        integer(HID_T)                          :: did
        integer(ik),            allocatable     :: indices(:)
        character(len=1024),    allocatable     :: names(:)
        integer(ik)                             :: idom, ndomains, ierr
        integer, dimension(1)                   :: buf
        integer(HSIZE_T)                        :: adim
        logical                                 :: attribute_exists

        !
        ! Get number of domains
        !
        ndomains = get_ndomains_hdf(fid)
        names    = get_domain_names_hdf(fid)

        !
        ! Allocate indices
        !
        allocate(indices(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        !  Loop through groups and read domain names
        !
        idom = 1
        do idom = 1,ndomains
            !
            ! Open domain group
            !
            call h5gopen_f(fid,trim(adjustl(names(idom))), did, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_indices_hdf: error opening domain group")

            !
            ! Get idomain attribute from fid/domain/idomain
            ! TODO: fix attribute discover. This causes a warning/error if not found.
            !
            call h5aexists_f(did, 'idomain', attribute_exists, ierr)
            !call h5ltget_attribute_int_f(fid,trim(adjustl(names(idom))),'idomain',buf,ierr)

            
            !
            ! If it doesn't exist, set to the current value of idom
            !
            adim = 1
            if ( .not. attribute_exists ) then

                ! Set value.
                call h5ltset_attribute_int_f(fid, trim(adjustl(names(idom))), 'idomain', [idom], adim, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_domain_indices_hdf: error writing an initial domain index")

            end if


            !
            ! Get value that was just set to be sure. 
            !
            call h5ltget_attribute_int_f(fid, trim(adjustl(names(idom))), 'idomain', buf, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_indices_hdf: error retrieving domain indices")

            !
            ! Set value detected to indices array that will be passed back from the function
            !
            indices(idom) = buf(1)


            !
            ! Close domain
            !
            call h5gclose_f(did,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_domain_indices: h5gclose")

        end do


    end function get_domain_indices_hdf
    !*************************************************************************************************************




















    !> Returns an array of integers that specifies the order of the coordinate expansion for every domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid         HDF file identifier.
    !!  @param[in]  dnames(:)   List of domain names to be interrogated. 
    !!
    !-------------------------------------------------------------------------------------------------------------
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
    !*************************************************************************************************************







    !> Returns an array of integer that specifies the order of the solution expansion for every domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier.
    !!  @param[in]  dnames  List of domain names to be interrogated.
    !!
    !-------------------------------------------------------------------------------------------------------------
    function get_order_solution_hdf(fid, dnames) result(orders)
        integer(HID_T),         intent(in)  :: fid
        character(len=1024),    intent(in)  :: dnames(:)

        integer(HID_T)              :: did
        integer(HSIZE_T)            :: adim
        logical                     :: order_exists
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
            ! Open domain group
            !
            call h5gopen_f(fid, trim(dnames(idom)), did, ierr)
            if (ierr /= 0) call chidg_signal_one(FATAL,"get_order_solution_hdf: error opening domain group.", trim(dnames(idom)) )

            
            !
            ! Check 'order_solution' attribute exists
            !
            call h5aexists_f(did, 'order_solution', order_exists, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_order_solution_hdf: error check attribue exists.")



            !
            ! Handle attribute does not exist
            !
            if ( .not. order_exists ) then
                adim = 1
                buf = 0
                call h5ltset_attribute_int_f(did, ".", 'order_solution', buf, adim, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_order_solution_hdf: error setting attribute.")
            end if



            !
            !  Get solution order
            !
            !call h5ltget_attribute_int_f(fid, trim(dnames(idom)), 'order_solution', buf, ierr)
            call h5ltget_attribute_int_f(did, ".", 'order_solution', buf, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_order_solution_hdf: error getting order_solution attribute.")


            !
            ! Compute number of terms in coordinate expansion
            !
            order = buf(1)
            orders(idom) = int(order, kind=ik)


            call h5gclose_f(did,ierr)

        end do

    end function get_order_solution_hdf
    !*************************************************************************************************************














    !>  Returns an array of integers that specifies the number of spatial dimensions to use for every domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!  @param[in]  fid         HDF file identifier.
    !!  @param[in]  dnames(:)   List of domain names to be interrogated. 
    !!
    !-------------------------------------------------------------------------------------------------------------
    function get_spacedim_hdf(fid, dnames) result(spacedims)
        integer(HID_T),         intent(in)  :: fid
        character(len=1024),    intent(in)  :: dnames(:)

        integer(ik), allocatable    :: spacedims(:)
        integer                     :: ierr, idom, spacedim
        integer, dimension(1)       :: buf


        !
        ! Allocate storage for orders
        !
        allocate(spacedims(size(dnames)), stat=ierr)
        if (ierr /= 0) call AllocationError

        !
        !  Loop through groups and read domains
        !
        do idom = 1,size(dnames)

            !
            !  Get coordinate mapping
            !
            call h5ltget_attribute_int_f(fid, trim(dnames(idom)), 'spacedim', buf, ierr)
            if (ierr /= 0) stop "Error: get_spacedim_hdf - h5ltget_attribute_int_f"


            !
            ! Compute number of terms in coordinate expansion
            !
            spacedim = buf(1)

            spacedims(idom) = int(spacedim, kind=ik)

        end do

    end function get_spacedim_hdf
    !*************************************************************************************************************


























    !> Returns an array of integer that specifies the order of the solution expansion for every domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid     HDF file identifier.
    !!  @param[in]  dnames  List of domain names to be interrogated.
    !!
    !------------------------------------------------------------------------------------------------------------
    function get_eqnset_hdf(fid, dnames) result(eqnsets)
        integer(HID_T),         intent(in)  :: fid
        character(len=1024),    intent(in)  :: dnames(:)

        integer(HID_T)                      :: did
        logical                             :: eqnset_exists
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
            ! Open domain
            !
            call h5gopen_f(fid, trim(dnames(idom)), did, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_eqnset_hdf: error opening domain group.")


            !
            ! Check eqnset attribute exists.
            !
            call h5aexists_f(did, 'eqnset', eqnset_exists, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_eqnset_hdf: error checking if 'eqnset' exists.")


            !
            ! If eqnset doesn't exists, create attribute and set to 'empty'
            !
            if ( .not. eqnset_exists ) then
                !call h5ltset_attribute_string_f(fid, trim(dnames(idom)), 'eqnset', 'empty', ierr)
                call h5ltset_attribute_string_f(did, ".", 'eqnset', 'empty', ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"get_eqnset_hdf: error setting empty 'eqnset' attribute.")
            end if


            !
            !  Get eqnset string from hdf attribute.
            !
            !call h5ltget_attribute_string_f(fid, trim(dnames(idom)), 'eqnset', eqnsets(idom), ierr)
            call h5ltget_attribute_string_f(did, ".", 'eqnset', eqnsets(idom), ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_eqnset_hdf - h5ltget_attribute_int_f.")


            call h5gclose_f(did,ierr)

        end do

    end function get_eqnset_hdf
    !*********************************************************************************************************











    !>  Return the boundary condition names from the HDF file that are set for each face of a particular domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  fid         ChiDG HDF5 file identifier.
    !!  @param[in]  dname       String indicating the boundary condition to query.
    !!  @result     bcnames     Array of boundary condition names for the specified domain. bcnames(ibc)
    !!
    !---------------------------------------------------------------------------------------------------------
    function get_bcnames_hdf(fid,dname) result(bcnames)
        integer(HID_T),         intent(in)  :: fid
        character(len=1024),    intent(in)  :: dname

        integer(HID_T)                      :: bc_id, bcface_id
        integer(ik)                         :: ndom, ierr, idom, iface
        character(len=1024), allocatable    :: bcnames(:)
        character(len=10)                   :: faces(NFACES)
        logical                             :: exists, bcname_exists



        !
        ! Get file information
        !
        ndom   = get_ndomains_hdf(fid)

        allocate(bcnames(NFACES), stat=ierr)
        if (ierr /= 0) call AllocationError




        !
        ! Check if 'BoundaryConditions' group exists
        !
        call h5lexists_f(fid, "/"//trim(dname)//"/BoundaryConditions", exists, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5lexists - checking BoundaryConditions")



        !
        ! Open the Domain/BoundaryConditions group
        !
        if (exists) then
            call h5gopen_f(fid, "/"//trim(dname)//"/BoundaryConditions", bc_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5gopen - BoundaryConditions")
        else
            call h5gcreate_f(fid, '/'//trim(dname)//"/BoundaryConditions", bc_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5gcreate - BoundaryConditions")
        end if


        !
        ! Loop through faces
        !
        faces = ["  XI_MIN","  XI_MAX"," ETA_MIN"," ETA_MAX","ZETA_MIN","ZETA_MAX"]
        do iface = 1,NFACES

            !
            ! Open boundary condition face group
            !
            call h5gopen_f(bc_id, trim(adjustl(faces(iface))), bcface_id, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5gopen_f - bcface (ex. XI_MIN)")

            !
            ! Check if boundary condition name attribute exists
            !
            call h5aexists_f(bcface_id, "bc_name", bcname_exists, ierr)


            if ( bcname_exists ) then
                call h5ltget_attribute_string_f(bc_id, trim(adjustl(faces(iface))), "bc_name", bcnames(iface), ierr)
            else
                call h5ltset_attribute_string_f(bc_id, trim(adjustl(faces(iface))), "bc_name", "empty", ierr)
                call h5ltget_attribute_string_f(bc_id, trim(adjustl(faces(iface))), "bc_name", bcnames(iface), ierr)
            end if
            if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5ltget_attribute_string")

    
            
            !
            ! Close boundary condition face group
            !
            call h5gclose_f(bcface_id,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5gclose.")

        end do


        !
        ! Close the boundary condition group
        ! 
        call h5gclose_f(bc_id,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"get_bcnames_hdf: h5gclose.")


    end function get_bcnames_hdf
    !**********************************************************************************************************














    !>  Delete all the attributes attached to a specified group identifier.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @param[in]  gid     HDF5 group identifier.
    !!
    !----------------------------------------------------------------------------------------------------------
    subroutine delete_group_attributes(gid)
        integer(HID_T),         intent(in)      :: gid

        integer(ik)                             :: nattr, ierr
        integer(HSIZE_T)                        :: iattr, idx
        type(h5o_info_t), target                :: h5_info


        !
        ! Get number of attributes attached to the group id
        !
        call h5oget_info_f(gid, h5_info, ierr)
        nattr = h5_info%num_attrs
        if (ierr /= 0) call chidg_signal(FATAL,"delete_group_attributes: error getting current number of attributes.")


        !
        ! Delete any existing attributes
        !
        if ( nattr > 0 ) then

            !
            ! Delete by index. h5adelete_by_idx_f returns idx with the next index so it doesn't need manually updated.
            !
            idx = 0
            do iattr = 1,nattr
                call h5adelete_by_idx_f(gid, ".", H5_INDEX_CRT_ORDER_F, H5_ITER_NATIVE_F, idx, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"delete_group_attributes: error deleting attribute")
            end do


        end if ! nattr


    end subroutine delete_group_attributes
    !**********************************************************************************************************















end module mod_hdf_utilities
