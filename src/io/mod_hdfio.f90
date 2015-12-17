module mod_hdfio
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO
    use type_meshdata,      only: meshdata_t
    use type_chidg_data,    only: chidg_data_t
    use mod_hdf_utilities,  only: get_ndomains_hdf
    use hdf5
    use h5lt
    use mod_io, only: nterms_s, eqnset
    implicit none



contains




    !>  Read HDF5 grid
    !!
    !!  Opens an hdf5 file, finds how many grid domains exist, allocates that number of
    !!  domains for initialization, and calls the geometry initialization routine, passing the
    !!  points read from the hdf5 file for initialization.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      filename    Character string of the file to be read
    !!  @param[inout]   domains     Allocatable array of domains. Allocated in this routine.
    !!
    !-----------------------------------------------------------------------------------------------------------------------------------
    subroutine read_grid_hdf(filename, meshdata)
        use mod_io, only: nterms_s
        character(*),                   intent(in)      :: filename
        type(meshdata_t), allocatable,  intent(inout)   :: meshdata(:)

        integer(HID_T)   :: fid, gid, sid, did_x, did_y, did_z      ! Identifiers
        integer(HSIZE_T) :: dims(3), maxdims(3)                     ! Dataspace dimensions

        type(c_ptr)                                     :: pts
        real(rk), dimension(:,:,:), allocatable, target :: xpts, ypts, zpts
        type(c_ptr)                                     :: cp_pts

        character(10)                           :: gname
        integer                                 :: nmembers, type, ierr, ndomains, igrp,    &
                                                   npts, izeta, ieta, ixi, idom, nterms_1d, &
                                                   mapping, nterms_c
        integer, dimension(1)                   :: buf
        logical                                 :: FileExists



        !
        !  Check file exists
        !
        inquire(file=filename, exist=FileExists)
        if (.not. FileExists) then
            call chidg_signal(FATAL,'read_grid_hdf5: Could not find grid file')
        end if


        !
        !  Initialize Fortran interface.
        !
        call h5open_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'read_grid_hdf5 - h5open_f: HDF5 Fortran interface had an error during initialization')



        !
        !  Open input file using default properties.
        !
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'read_grid_hdf5 - h5fopen_f: There was an error opening the grid file.')



        !
        !  Get number of domains from attribute 'ndomains' in file root
        !
        ndomains = get_ndomains_hdf(fid)
        !call h5ltget_attribute_int_f(fid, "/", 'ndomains', buf, ierr)
        !ndomains = buf(1)
        !if (ierr /= 0) call chidg_signal(FATAL,'read_grid_hdf5: h5ltget_attribute_int_f had a problem getting the number of domains')


        !
        !  Allocate number of domains
        !
        if (ndomains == 0) call chidg_signal(FATAL,'read_hdf5: No Domains were found in the file')
        allocate(meshdata(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError





        !  Get number of groups in the file root
        call h5gn_members_f(fid, "/", nmembers, ierr)

        !  Loop through groups and read domains
        idom = 1
        do igrp = 0,nmembers-1
            call h5gget_obj_info_idx_f(fid,"/", igrp, gname, type, ierr)

            if (gname(1:2) == 'D_') then

                ! Open the Domain/Grid group
                call h5gopen_f(fid, trim(gname)//"/Grid", gid, ierr, H5P_DEFAULT_F)
                if (ierr /= 0) stop "Error: read_grid_hdf5 -- h5gopen_f: Domain/Grid group did not open properly"


                !  Get number of terms in coordinate expansion
                call h5ltget_attribute_int_f(fid, trim(gname), 'mapping', buf, ierr)


                mapping = buf(1)
                if (ierr /= 0) stop "Error: read_grid_hdf5 - h5ltget_attribute_int_f"
                nterms_1d = (mapping + 1)
                nterms_c = nterms_1d * nterms_1d * nterms_1d

                meshdata(idom)%nterms_c = nterms_c
                meshdata(idom)%name     = gname


                !  Open the Coordinate datasets
                call h5dopen_f(gid, "CoordinateX", did_x, ierr, H5P_DEFAULT_F)
                call h5dopen_f(gid, "CoordinateY", did_y, ierr, H5P_DEFAULT_F)
                call h5dopen_f(gid, "CoordinateZ", did_z, ierr, H5P_DEFAULT_F)


                !  Get the dataspace id and dimensions
                call h5dget_space_f(did_x, sid, ierr)
                call h5sget_simple_extent_dims_f(sid, dims, maxdims, ierr)


                !  Read x-points
                allocate(xpts(dims(1),dims(2),dims(3)))
                cp_pts = c_loc(xpts(1,1,1))
                call h5dread_f(did_x, H5T_NATIVE_DOUBLE, cp_pts, ierr)
                if (ierr /= 0) stop "Error: read_grid_hdf5 -- h5dread_f"


!                allocate(ypts, mold=xpts)   ! bug in gcc
                allocate(ypts(dims(1),dims(2),dims(3)))
                cp_pts = c_loc(ypts(1,1,1))
                call h5dread_f(did_y, H5T_NATIVE_DOUBLE, cp_pts, ierr)
                if (ierr /= 0) stop "Error: read_grid_hdf5 -- h5dread_f"


!                allocate(zpts, mold=xpts)   ! bug in gcc
                allocate(zpts(dims(1),dims(2),dims(3)))
                cp_pts = c_loc(zpts(1,1,1))
                call h5dread_f(did_z, H5T_NATIVE_DOUBLE, cp_pts, ierr)
                if (ierr /= 0) stop "Error: read_grid_hdf5 -- h5dread_f"


                !  Accumulate pts into a single points_t matrix to initialize domain
                npts = dims(1)*dims(2)*dims(3)
                allocate(meshdata(idom)%points(dims(1),dims(2),dims(3)), stat=ierr)
                if (ierr /= 0) call AllocationError
                    
                do izeta = 1,dims(3)
                    do ieta = 1,dims(2)
                        do ixi = 1,dims(1)
                            call meshdata(idom)%points(ixi,ieta,izeta)%set(xpts(ixi,ieta,izeta),ypts(ixi,ieta,izeta),zpts(ixi,ieta,izeta))
                        end do
                    end do
                end do


                ! Read equation set attribute
!                call h5ltget_attribute_string_f(fid,trim(gname),'EquationSet',eqnset,ierr)
!                if (ierr /= 0) stop "Error: read_grid_hdf5: h5ltget_attribute_string_f -- EquationSet"




                !
                ! Close the Coordinate datasets
                !
                call h5dclose_f(did_x,ierr)
                call h5dclose_f(did_y,ierr)
                call h5dclose_f(did_z,ierr)


                !
                ! Close the dataspace id
                !
                call h5sclose_f(sid,ierr)


                !
                ! Close the Domain/Grid group
                !
                call h5gclose_f(gid,ierr)

                ! Deallocate points for the current domain
                deallocate(zpts,ypts,xpts)
                idom = idom + 1
            end if
        end do


        !
        !  Close file and Fortran interface
        !
        call h5fclose_f(fid, ierr)
        call h5close_f(ierr)

    end subroutine read_grid_hdf
    !#########################################################################################################################################


























    !> Read HDF5 variable
    !!
    !!  Opens a given hdf5 file. Loads the EquationSet and solution order and calls solution initialization
    !!  procedure for each domain. Searches for the given variable and time instance. If it finds it, load to a
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      fid         HDF5 file identifier.
    !!  @param[in]      varstring   Character string of the variable name to be read.
    !!  @param[in]      itime       Integer of the time instance for the current variable to be read.
    !!  @param[in]      dname       Character string of the domain to be read from.
    !!  @param[inout]   data        ChiDG data containing domains. Already allocated.
    !--------------------------------------------------------------------------------------------------------------------------------------
    subroutine read_variable_hdf(fid,varstring,itime,dname,data)
        use ISO_C_BINDING
        integer(HID_T),     intent(in)      :: fid
        character(*),       intent(in)      :: varstring
        integer(ik),        intent(in)      :: itime
        character(*),       intent(in)      :: dname
        type(chidg_data_t), intent(inout)   :: data


        integer(HID_T)   :: gid, sid, vid           ! Identifiers
        integer(HSIZE_T) :: maxdims(2)              ! Dataspace dimensions
        integer(HSIZE_T) :: dims(3)

        integer, dimension(1)           :: ibuf

        character(100)                  :: cbuf
        character(100)                  :: var_gqp

        real(rk), allocatable, target   :: var(:,:,:)
        real(rk), allocatable           :: bufferterms(:)
        type(c_ptr)                     :: cp_var

        integer                         :: type,    ierr,       igrp,               &
                                           npts,    nterms_1d,  nterms_s,   order,  &
                                           ivar,    ielem,      nterms_ielem,   idom
        logical                         :: ElementsEqual




        !
        ! Get name of current domain
        !
        idom = data%get_domain_index(dname)



        !
        ! Open the Domain/Variables group
        !
        call h5gopen_f(fid, trim(dname)//"/Variables", gid, ierr, H5P_DEFAULT_F)
        if (ierr /= 0) call chidg_signal(FATAL,"h5gopen_f -- Domain/Grid group did not open properly")


        !
        ! Get number of terms in solution expansion
        !
        call h5ltget_attribute_int_f(fid, trim(dname), 'order_solution', ibuf, ierr)

        order = ibuf(1)
        if (ierr /= 0) call chidg_signal(FATAL,"read_variable_hdf5 - h5ltget_attribute_int_f")
        nterms_1d = (order + 1) ! To be consistent with the definition of (Order = 'Order of the polynomial')
        nterms_s = nterms_1d*nterms_1d*nterms_1d


        
        !
        ! Open the Variable dataset
        !
        call h5dopen_f(gid, trim(varstring), vid, ierr, H5P_DEFAULT_F)
        if (ierr /= 0) call chidg_signal(FATAL,"read_variable_hdf5 -- variable does not exist or was not opened correctly")


        !
        ! Get the dataspace id and dimensions
        !
        call h5dget_space_f(vid, sid, ierr)
        call h5sget_simple_extent_dims_f(sid, dims, maxdims, ierr)


        !
        ! Read 'variable' dataset
        !
        allocate(var(dims(1),dims(2),dims(3)), stat=ierr)               ! Allocate variable buffer
        if ( ierr /= 0 ) call AllocationError

        cp_var = c_loc(var(1,1,1))                                      ! Get C-address for buffer

        call h5dread_f(vid, H5T_NATIVE_DOUBLE, cp_var, ierr)            ! Fortran 2003 interface
        if (ierr /= 0) call chidg_signal(FATAL,"read_variable_hdf5 -- h5dread_f")


        !
        !  Get variable index in EquationSet
        !
        ivar = data%eqnset(idom)%item%prop%get_eqn_index(trim(varstring))


        !
        !  Test to make sure the number of elements in the variable group
        !  and the current domain are conforming
        !
        ElementsEqual = (size(data%sdata%q%dom(idom)%lvecs) == size(var,2))
        if (ElementsEqual) then
            !
            !  Loop through elements and set 'variable' values
            !
            do ielem = 1,data%mesh(idom)%nelem
                !
                ! Get number of terms initialized for the current element
                !
                nterms_ielem = data%sdata%q%dom(idom)%lvecs(ielem)%nterms()


                !
                ! Allocate bufferterm storage that will be used to set variable data
                !
                if (allocated(bufferterms)) deallocate(bufferterms)
                allocate(bufferterms(nterms_ielem), stat=ierr)
                bufferterms = ZERO
                if (ierr /= 0) call AllocationError


                !
                ! Check for reading lower, higher, or same-order solution
                !
                if ( nterms_s < nterms_ielem ) then
                    bufferterms(1:nterms_s) = var(1:nterms_s, ielem, itime)             ! Reading a lower order solution
                else if ( nterms_s > nterms_ielem ) then
                    bufferterms(1:nterms_ielem) = var(1:nterms_ielem, ielem, itime)     ! Reading a higher-order solution
                else
                    bufferterms(1:nterms_ielem) = var(1:nterms_ielem, ielem, itime)     ! Reading a solution of same order
                end if



                call data%sdata%q%dom(idom)%lvecs(ielem)%setvar(ivar,bufferterms)
            end do

        else
            call chidg_signal(FATAL,"read_variable_hdf5 -- number of elements in file variable and domain do not match")
        end if




        call h5dclose_f(vid,ierr)       ! Close the variable dataset
        call h5gclose_f(gid,ierr)       ! Close the Domain/Variable group

    end subroutine read_variable_hdf
    !#######################################################################################################################################
























    !> Write HDF5 variable
    !!
    !!  Opens a given hdf5 file. Loads the EquationSet and solution order and calls solution initialization
    !!  procedure for each domain. Searches for the given variable and time instance. If it finds it, load to a
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]      fid         HDF5 file identifier.
    !!  @param[in]      varstring   Character string of the variable name to be read.
    !!  @param[in]      itime       Integer of the time instance for the current variable to be read.
    !!  @param[in]      dname       Character string of the domain name to be read from.
    !!  @param[inout]   data        chidg_data_t instance containing grid and solution.
    !-----------------------------------------------------------------------------------------------------------------------------------------
    subroutine write_variable_hdf(fid,varstring,itime,dname,data)
        integer(HID_T),     intent(in)      :: fid
        character(*),       intent(in)      :: varstring
        integer(ik),        intent(in)      :: itime
        character(*),       intent(in)      :: dname
        type(chidg_data_t), intent(inout)   :: data


        integer(HID_T)   :: gid, sid, did, crp_list         ! Identifiers
        integer(HSIZE_T) :: maxdims(2), adim                ! Dataspace dimensions
        integer(HSIZE_T) :: dims(3)                         ! Dataspace dimensions
        integer(HSIZE_T) :: dimsc(3)                        ! Chunk size for extendible data sets
        type(H5O_INFO_T) :: info                            ! Object info type

        integer                         :: ndims
        integer, dimension(1)           :: ibuf
        character(100)                  :: cbuf
        character(100)                  :: var_grp
        character(100)                  :: ctime

        real(rk), allocatable, target   :: var(:,:,:)
        type(c_ptr)                     :: cp_var

        integer(ik)                     :: nmembers,    type,   ierr,       ndomains,   igrp,   &
                                           npts,        nterms_1d,  nterms_s,   order,  &
                                           ivar,        ielem,      idom
        logical                         :: FileExists, VariablesExists, DataExists, ElementsEqual
        logical                         :: exists




        !
        ! Check if 'Variables' group exists
        !
        call h5lexists_f(fid, trim(dname)//"/Variables", exists, ierr)



        !
        ! Open the Domain/Variables group
        !
        if (exists) then

            ! If 'Variables' group exists then open the existing group
            call h5gopen_f(fid, trim(dname)//"/Variables", gid, ierr, H5P_DEFAULT_F)
            if (ierr /= 0) call chidg_signal(FATAL,"h5gopen_f -- Domain/Grid group did not open properly")
        else
            ! If 'Variables group does not exist, then create one.
            call h5gcreate_f(fid, trim(dname)//"/Variables", gid, ierr)
        end if




        !
        ! Set dimensions of dataspace to write
        !
        idom = data%get_domain_index(dname)
        ndims = 3

        dims(1) = data%mesh(idom)%nterms_s
        dims(2) = data%mesh(idom)%nelem
        dims(3) = itime                     ! TODO: Should probably better inform the dataspace dimension here. Probably set mesh_t%ntime
        maxdims(1) = H5S_UNLIMITED_F
        maxdims(2) = H5S_UNLIMITED_F
        maxdims(2) = H5S_UNLIMITED_F




        !
        ! Open the Variable dataset, check if Variable dataset already exists
        !
        call h5lexists_f(gid, trim(varstring), exists, ierr)



        !
        ! Modify dataset creation properties, i.e. enable chunking in order to append dataspace, if needed.
        !
        dimsc = [1, data%mesh(idom)%nelem, 1]  ! Chunk size

        call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, ierr)
        if (ierr /= 0) call chidg_signal(FATAL, "Error: write_variable_hdf5 -- h5pcreate_f error enabling chunking")

        call h5pset_chunk_f(crp_list, ndims, dimsc, ierr)
        if (ierr /= 0) call chidg_signal(FATAL, "Error: write_variable_hdf5 -- h5pset_chunk_f error setting chunk properties")



        !
        ! Reset dataspace size if necessary
        !
        if (exists) then
            ! Open the existing dataset
            call h5dopen_f(gid, trim(varstring), did, ierr, H5P_DEFAULT_F)
            if (ierr /= 0) call chidg_signal(FATAL,"Error: write_variable_hdf5 -- variable does not exist or was not opened correctly")


            ! Extend dataset if necessary
            call h5dset_extent_f(did, dims, ierr)
            if (ierr /= 0) call chidg_signal(FATAL, "Error: write_variable_hdf5 -- h5dset_extent_f")


            ! Update existing dataspace ID since it may have been expanded
            call h5dget_space_f(did, sid, ierr)
            if (ierr /= 0) call chidg_signal(FATAL, "Error: write_variable_hdf5 -- h5dget_space_f")

        else
            ! Create a new dataspace
            call h5screate_simple_f(ndims,dims, sid, ierr, maxdims)
            if (ierr /= 0) call chidg_signal(FATAL,"Error: write_variable_hdf5 - h5screate_simple_f")


            ! Create a new dataset
            call h5dcreate_f(gid, trim(varstring), H5T_NATIVE_DOUBLE, sid, did, ierr, crp_list)
            if (ierr /= 0) call chidg_signal(FATAL,"write_variable_hdf5 - h5dcreate_f")
        end if



        !
        ! Get variable integer index from variable character string
        !
        !ivar = data%eqnset(idom)%item%prop%get_eqn_index(cvar)
        ivar = data%eqnset(idom)%item%prop%get_eqn_index(varstring)



        !
        ! Assemble variable buffer matrix that gets written to file
        !
        allocate(var(dims(1),dims(2),dims(3)))

        do ielem = 1,data%mesh(idom)%nelem
                var(:,ielem,itime) = data%sdata%q%dom(idom)%lvecs(ielem)%getvar(ivar)
        end do



        !
        ! Write variable buffer
        !
        cp_var = c_loc(var(1,1,1))

        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, cp_var, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"write_variable_hdf5 - h5dwrite_f")


        !
        ! Update 'contains_solution' attribute
        !
        call h5ltset_attribute_string_f(fid, "/", 'contains_solution', 'Yes', ierr)


        call h5pclose_f(crp_list, ierr) ! Close dataset creation property
        call h5dclose_f(did,ierr)       ! Close Variable datasets
        call h5sclose_f(sid,ierr)       ! Close Variable dataspaces
        call h5gclose_f(gid,ierr)       ! Close Domain/Variable group



    end subroutine write_variable_hdf
    !####################################################################################################################################
































    !> Read solution modes from HDF file.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      filename    Character string of the file to be read from
    !!  @param[inout]   data        chidg_data_t that will accept the solution modes
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine read_solution_hdf(filename,data)
        use ISO_C_BINDING
        character(*),       intent(in)      :: filename
        type(chidg_data_t), intent(inout)   :: data

        integer(HID_T)                  :: fid
        integer                         :: ierr
        character(:),       allocatable :: dname

        integer(ik)                     :: idom, ndomains
        integer(ik)                     :: ieqn, neqns
        integer(ik)                     :: time
        character(len=:),   allocatable :: cvar
        logical                         :: fileexists = .false.


        !
        ! Get number of domains contained in the ChiDG data instance
        !
        ndomains = data%ndomains()


        !
        ! Set default time instance
        !
        time = 1




        !
        ! Check file exists
        !
        inquire(file=filename, exist=fileexists)
        if (.not. fileexists) call chidg_signal_one(FATAL,"Error: read_solution_hdf5 - file not found: ", filename)


        !
        ! Initialize Fortran interface.
        !
        call h5open_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"read_solution_hdf5 - h5open_f")


        !
        ! Open input file using default properties.
        !
        call h5fopen_f(filename, H5F_ACC_RDWR_F, fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"read_solution_hdf5 - h5fopen_f")














        !
        ! Read solution for each domain
        !
        do idom = 1,ndomains
            !
            ! Get name of current domain
            !
            dname = data%info(idom)%name
            

            !
            ! Get number of equations for the current domain
            !
            neqns = data%eqnset(idom)%item%neqns



            !
            ! For each equation specified for the domain
            ! 
            do ieqn = 1,neqns
                !
                ! Get variable character string
                !
                cvar = trim(data%eqnset(idom)%item%prop%eqns(ieqn)%name)

                !
                ! Read variable
                !
                !call read_variable_hdf(filename,cvar,time,idom,data)
                call read_variable_hdf(fid,cvar,time,trim(dname),data)
            end do ! ieqn



        end do ! idom




        !
        ! Close HDF5 file and Fortran interface
        !
        call h5fclose_f(fid, ierr)      ! Close the HDF5 file
        call h5close_f(ierr)            ! Close the HDF5 Fortran interface

    end subroutine read_solution_hdf
    !################################################################################################


















    !> Write solution modes to HDF file.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!  @param[in]      filename    Character string of the file to be written to
    !!  @param[inout]   data        chidg_data_t containing solution to be written
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine write_solution_hdf(filename,data)
        use ISO_C_BINDING
        character(*),       intent(in)      :: filename
        type(chidg_data_t), intent(inout)   :: data


        integer(HID_T)                  :: fid
        integer(HSIZE_T)                :: adim
        integer(ik)                     :: idom, ndomains
        integer(ik)                     :: ieqn, neqns
        integer(ik)                     :: time
        character(len=:),   allocatable :: cvar
        character(len=:),   allocatable :: dname
        integer                         :: ierr, order_s
        logical                         :: fileexists


        !
        ! Get number of domains contained in the ChiDG data instance
        !
        ndomains = data%ndomains()


        !
        ! Set default time instance
        !
        time = 1


        !
        ! Check file exists
        !
        inquire(file=filename, exist=fileexists)
        if (.not. fileexists) call chidg_signal_one(FATAL, "Error: write_solution_hdf - file not found: ", filename)


        !
        ! Initialize Fortran interface.
        !
        call h5open_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"write_solution_hdf - h5open_f")

        
        !
        ! Open input file using default properties.
        !
        call h5fopen_f(filename, H5F_ACC_RDWR_F, fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"write_solution_hdf - h5fopen_f")

        



        !
        ! Read solution for each domain
        !
        do idom = 1,ndomains
            !
            ! Get name of current domain
            !
            dname = data%info(idom)%name
            

            !
            ! Write domain attributes: solution order, equation set
            !
            adim = 1
            order_s = 0
            do while ( order_s*order_s*order_s /= data%mesh(idom)%nterms_s )
               order_s = order_s + 1 
            end do
            order_s = order_s - 1 ! to be consistent with he definition of 'Order of the polynomial'

            call h5ltset_attribute_int_f(fid, trim(dname), 'order_solution', [order_s], adim, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"write_variable_hdf5 - h5ltset_attribute_int_f")

            call h5ltset_attribute_string_f(fid, trim(dname), 'eqnset', trim(data%eqnset(idom)%item%name), ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"write_variable_hdf5 - h5ltset_attribute_int_f")


            !
            ! Get number of equations for the current domain
            !
            neqns = data%eqnset(idom)%item%neqns


            !
            ! For each equation specified for the domain
            ! 
            do ieqn = 1,neqns
                !
                ! Get variable character string
                !
                cvar = trim(data%eqnset(idom)%item%prop%eqns(ieqn)%name)

                !
                ! Write variable
                !
                call write_variable_hdf(fid,cvar,time,dname,data)
            end do ! ieqn

        end do ! idom



        !
        ! Close HDF5 file and Fortran interface
        !
        call h5fclose_f(fid, ierr)      ! Close HDF5 File
        call h5close_f(ierr)            ! Close HDF5 Fortran interface

    end subroutine write_solution_hdf
    !##################################################################################################















end module mod_hdfio
