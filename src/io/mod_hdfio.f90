module mod_hdfio
#include <messenger.h>
    use mod_kinds,                  only: rk,ik,rdouble
    use mod_constants,              only: ZERO, NFACES, TWO_DIM, THREE_DIM, NO_PROC
    use mod_bc,                     only: create_bc
    use mod_hdf_utilities,          only: get_ndomains_hdf, get_domain_names_hdf, get_eqnset_hdf, get_domain_indices_hdf, get_domain_name_hdf
    use mod_io,                     only: nterms_s
    use mod_chidg_mpi,              only: IRANK

    use type_chidg_data,            only: chidg_data_t
    use type_meshdata,              only: meshdata_t
    use type_bcdata,                only: bcdata_t
    use type_bc,                    only: bc_t
    use type_domain_connectivity,   only: domain_connectivity_t
    use type_partition,             only: partition_t
    use hdf5
    use h5lt
    implicit none



contains




    !>  Read HDF5 grid
    !!
    !!  Opens an hdf5 file, finds how many grid domains exist, allocates that number of
    !!  domains for initialization, and calls the geometry initialization routine, passing the
    !!  points read from the hdf5 file for initialization.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      filename    Character string of the file to be read
    !!  @param[inout]   domains     Allocatable array of domains. Allocated in this routine.
    !!
    !---------------------------------------------------------------------------------------------------------------------
    subroutine read_grid_hdf(filename, meshdata)
        use mod_io, only: nterms_s
        character(*),                   intent(in)      :: filename
        type(meshdata_t), allocatable,  intent(inout)   :: meshdata(:)

        integer(HID_T)   :: fid, gid, sid, did_x, did_y, did_z, did_e               ! Identifiers
        integer(HSIZE_T) :: rank_one_dims(1), rank_two_dims(2), dims(3), maxdims(3) ! Dataspace dimensions

        type(c_ptr)                                         :: pts
        real(rdouble), dimension(:), allocatable, target    :: xpts, ypts, zpts
        integer,                     allocatable, target    :: connectivity(:,:)
        type(c_ptr)                                         :: cp_pts, cp_conn

        character(len=1024),    allocatable     :: dnames(:), eqnset(:)
        character(1024)                         :: gname
        integer                                 :: nmembers, type, ierr, ndomains, igrp,    &
                                                   npts, izeta, ieta, ixi, idom, nterms_1d, &
                                                   mapping, nterms_c, spacedim, ipt, idomain, ielem, nelements, nnodes
        integer, dimension(1)                   :: mapping_buf, spacedim_buf, idomain_buf
        logical                                 :: FileExists



        !
        !  Check file exists
        !
        inquire(file=filename, exist=FileExists)
        if (.not. FileExists) then
            call chidg_signal_one(FATAL,'read_grid_hdf5: Could not find grid file',filename)
        end if


        !
        !  Initialize Fortran interface.
        !
        call h5open_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'read_grid_hdf5 - h5open_f: HDF5 Fortran interface had an error during initialization')



        !
        !  Open input file using default properties.
        !
        !call h5fopen_f(filename, H5F_ACC_RDONLY_F, fid, ierr)
        call h5fopen_f(filename, H5F_ACC_RDWR_F, fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'read_grid_hdf5 - h5fopen_f: There was an error opening the grid file.')



        !
        !  Get number of domains from attribute 'ndomains' in file root
        !
        ndomains = get_ndomains_hdf(fid)


        !
        !  Allocate number of domains
        !
        if (ndomains == 0) call chidg_signal(FATAL,'read_grid_hdf5: No Domains were found in the file')
        allocate(meshdata(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError




!        !
!        !  Get number of groups in the file root
!        !
!        call h5gn_members_f(fid, "/", nmembers, ierr)
!        if (ierr /= 0) call chidg_signal(FATAL,"read_grid_hdf5: error getting number of groups in file root.")



        !
        ! Get equationset strings.
        !
        dnames   = get_domain_names_hdf(fid)
        eqnset   = get_eqnset_hdf(fid,dnames)




        !
        !  Loop through groups and read domains
        !
        do idom = 1,size(dnames)


                gname = dnames(idom)

                !
                ! Open the Domain/Grid group
                !
                call h5gopen_f(fid, trim(gname)//"/Grid", gid, ierr, H5P_DEFAULT_F)
                if (ierr /= 0) call chidg_signal_one(FATAL,"read_grid_hdf: Domagin/Grid group did not open properly.", trim(gname)//'/Grid')



!                !
!                !  Get number of terms in coordinate expansion
!                !
!                call h5ltget_attribute_int_f(fid, trim(gname), 'idomain', idomain_buf, ierr)
!                meshdata(idom)%idomain = idomain_buf(1)
!                if (ierr /= 0) call chidg_signal(FATAL,"Error: read_grid_hdf5 - h5ltget_attribute_int_f")



                !
                !  Get number of terms in coordinate expansion
                !
                call h5ltget_attribute_int_f(fid, trim(gname), 'mapping', mapping_buf, ierr)
                mapping = mapping_buf(1)
                if (ierr /= 0) call chidg_signal(FATAL,"Error: read_grid_hdf5 - h5ltget_attribute_int_f")
                nterms_1d = (mapping + 1)


                call h5ltget_attribute_int_f(fid, trim(gname), 'spacedim', spacedim_buf, ierr)
                spacedim = spacedim_buf(1)


                if ( spacedim == THREE_DIM ) then
                    nterms_c = nterms_1d * nterms_1d * nterms_1d
                else if ( spacedim == TWO_DIM ) then
                    nterms_c = nterms_1d * nterms_1d
                end if

                meshdata(idom)%nterms_c = nterms_c
                meshdata(idom)%name     = gname
                meshdata(idom)%spacedim = spacedim


                !
                !  Open the Coordinate datasets
                !
                call h5dopen_f(gid, "CoordinateX", did_x, ierr, H5P_DEFAULT_F)
                call h5dopen_f(gid, "CoordinateY", did_y, ierr, H5P_DEFAULT_F)
                call h5dopen_f(gid, "CoordinateZ", did_z, ierr, H5P_DEFAULT_F)


                !
                !  Get the dataspace id and dimensions
                !
                call h5dget_space_f(did_x, sid, ierr)
                call h5sget_simple_extent_dims_f(sid, rank_one_dims, maxdims, ierr)
                npts = rank_one_dims(1)


                !
                !  Read x-points
                !
                allocate(xpts(npts),stat=ierr)
                cp_pts = c_loc(xpts(1))
                call h5dread_f(did_x, H5T_NATIVE_DOUBLE, cp_pts, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"Error: read_grid_hdf5 -- h5dread_f")


!                allocate(ypts, mold=xpts)   ! bug in gcc
                allocate(ypts(npts))
                cp_pts = c_loc(ypts(1))
                call h5dread_f(did_y, H5T_NATIVE_DOUBLE, cp_pts, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"Error: read_grid_hdf5 -- h5dread_f")


!                allocate(zpts, mold=xpts)   ! bug in gcc
                allocate(zpts(npts))
                cp_pts = c_loc(zpts(1))
                call h5dread_f(did_z, H5T_NATIVE_DOUBLE, cp_pts, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"Error: read_grid_hdf5 -- h5dread_f")


                !
                !  Accumulate pts into a single points_t matrix to initialize domain
                !
                allocate(meshdata(idom)%points(npts), stat=ierr)
                if (ierr /= 0) call AllocationError
                    

                do ipt = 1,rank_one_dims(1)
                    call meshdata(idom)%points(ipt)%set(real(xpts(ipt),rk),real(ypts(ipt),rk),real(zpts(ipt),rk))
                end do





                !
                ! Open Elements connectivity dataset
                !
                call h5dopen_f(gid, "Elements", did_e, ierr, H5P_DEFAULT_F)

                !
                !  Get the dataspace id and dimensions
                !
                call h5dget_space_f(did_e, sid, ierr)
                call h5sget_simple_extent_dims_f(sid, rank_two_dims, maxdims, ierr)
                if (allocated(connectivity)) deallocate(connectivity)
                allocate(connectivity(rank_two_dims(1),rank_two_dims(2)))


                !
                ! Read connectivity
                ! 
                cp_conn = c_loc(connectivity(1,1))
                call h5dread_f(did_e, H5T_NATIVE_INTEGER, cp_conn, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"Error: read_grid_hdf5 -- h5dread_f")
                !meshdata(idom)%connectivity = connectivity




                ! Initialize domain connectivity structure
!                idomain   = meshdata(idom)%idomain  !prob don't need this
                nelements = size(connectivity,1)
                nnodes    = size(meshdata(idom)%points)
                call meshdata(idom)%connectivity%init(nelements, nnodes)


                !connectivities(idom)%data = connectivity
                do ielem = 1,nelements
                    meshdata(idom)%connectivity%data(ielem)%data = connectivity(ielem,:)
!                    call meshdata(idom)%connectivity%data(ielem)%set_element_nodes(connectivity(ielem,:))
                    call meshdata(idom)%connectivity%data(ielem)%set_element_partition(IRANK)
                end do






                !
                ! Read equation set attribute
                !
                meshdata(idom)%eqnset = eqnset(idom)


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


        end do  ! igrp


        !
        !  Close file and Fortran interface
        !
        call h5fclose_f(fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"read_grid_hdf: error closing file.")
        call h5close_f(ierr)

    end subroutine read_grid_hdf
    !***************************************************************************************************************















   
   
   
   
   
    !>  Read HDF5 grid
    !!
    !!  Opens an hdf5 file, finds how many grid domains exist, allocates that number of
    !!  domains for initialization, and calls the geometry initialization routine, passing the
    !!  points read from the hdf5 file for initialization.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      filename    Character string of the file to be read
    !!  @param[inout]   domains     Allocatable array of domains. Allocated in this routine.
    !!
    !---------------------------------------------------------------------------------------------------------------------
    subroutine read_grid_partition_hdf(filename, partition, meshdata)
        use mod_io, only: nterms_s
        character(*),                   intent(in)      :: filename
        type(partition_t),              intent(in)      :: partition
        type(meshdata_t), allocatable,  intent(inout)   :: meshdata(:)

        integer(HID_T)   :: fid, gid, sid, did_x, did_y, did_z, did_e               ! Identifiers
        integer(HSIZE_T) :: rank_one_dims(1), rank_two_dims(2), dims(3), maxdims(3) ! Dataspace dimensions

        type(c_ptr)                                         :: pts
        real(rdouble), dimension(:), allocatable, target    :: xpts, ypts, zpts
        type(c_ptr)                                         :: cp_pts, cp_conn

        character(len=1024),    allocatable     :: dnames(:), eqnset(:)
        character(1024)                         :: gname
        integer                                 :: nmembers, type, ierr, ndomains, igrp,    &
                                                   npts, izeta, ieta, ixi, idom, nterms_1d, &
                                                   mapping, nterms_c, spacedim, ipt, iconn, nconn, nelements, nnodes
        integer, dimension(1)                   :: mapping_buf, spacedim_buf
        logical                                 :: FileExists



        !
        !  Check file exists
        !
        inquire(file=filename, exist=FileExists)
        if (.not. FileExists) then
            call chidg_signal_one(FATAL,'read_grid_hdf5: Could not find grid file',filename)
        end if


        !
        !  Initialize Fortran interface.
        !
        call h5open_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'read_grid_hdf5 - h5open_f: HDF5 Fortran interface had an error during initialization')



        !
        !  Open input file using default properties.
        !
        !call h5fopen_f(filename, H5F_ACC_RDONLY_F, fid, ierr)
        call h5fopen_f(filename, H5F_ACC_RDWR_F, fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'read_grid_hdf5 - h5fopen_f: There was an error opening the grid file.')


        !
        !  Allocate a meshdata structure for each domain connectivity in the partition
        !
        nconn = size(partition%connectivities)
        allocate(meshdata(nconn), stat=ierr)
        if (ierr /= 0) call AllocationError




        !
        ! Get equationset strings.
        !
        dnames   = get_domain_names_hdf(fid)
        eqnset   = get_eqnset_hdf(fid,dnames)


        !
        !  Loop through groups and read domains
        !
        do iconn = 1,nconn


            ! Get connectivity domain index 
            idom = partition%connectivities(iconn)%get_domain_index()
            !meshdata(iconn)%idomain = idom

            
            ! Get the name of the current domain from the HDF file
            gname = get_domain_name_hdf(fid,idom)


            !gname = dnames(idom)

            !
            ! Open the Domain/Grid group
            !
            call h5gopen_f(fid, trim(gname)//"/Grid", gid, ierr, H5P_DEFAULT_F)
            if (ierr /= 0) call chidg_signal_one(FATAL,"read_grid_hdf: Domagin/Grid group did not open properly.", trim(gname)//'/Grid')


            !
            !  Get number of terms in coordinate expansion
            !
            call h5ltget_attribute_int_f(fid, trim(gname), 'mapping', mapping_buf, ierr)
            mapping = mapping_buf(1)
            if (ierr /= 0) stop "Error: read_grid_hdf5 - h5ltget_attribute_int_f"
            nterms_1d = (mapping + 1)





            call h5ltget_attribute_int_f(fid, trim(gname), 'spacedim', spacedim_buf, ierr)
            spacedim = spacedim_buf(1)


            if ( spacedim == THREE_DIM ) then
                nterms_c = nterms_1d * nterms_1d * nterms_1d
            else if ( spacedim == TWO_DIM ) then
                nterms_c = nterms_1d * nterms_1d
            end if

            meshdata(iconn)%nterms_c = nterms_c
            meshdata(iconn)%name     = gname
            meshdata(iconn)%spacedim = spacedim


            !
            !  Open the Coordinate datasets
            !
            call h5dopen_f(gid, "CoordinateX", did_x, ierr, H5P_DEFAULT_F)
            call h5dopen_f(gid, "CoordinateY", did_y, ierr, H5P_DEFAULT_F)
            call h5dopen_f(gid, "CoordinateZ", did_z, ierr, H5P_DEFAULT_F)


            !
            !  Get the dataspace id and dimensions
            !
            call h5dget_space_f(did_x, sid, ierr)
            call h5sget_simple_extent_dims_f(sid, rank_one_dims, maxdims, ierr)
            npts = rank_one_dims(1)


            !
            !  Read x-points
            !
            allocate(xpts(npts),stat=ierr)
            cp_pts = c_loc(xpts(1))
            call h5dread_f(did_x, H5T_NATIVE_DOUBLE, cp_pts, ierr)
            if (ierr /= 0) stop "Error: read_grid_hdf5 -- h5dread_f"


!                allocate(ypts, mold=xpts)   ! bug in gcc
            allocate(ypts(npts))
            cp_pts = c_loc(ypts(1))
            call h5dread_f(did_y, H5T_NATIVE_DOUBLE, cp_pts, ierr)
            if (ierr /= 0) stop "Error: read_grid_hdf5 -- h5dread_f"


!                allocate(zpts, mold=xpts)   ! bug in gcc
            allocate(zpts(npts))
            cp_pts = c_loc(zpts(1))
            call h5dread_f(did_z, H5T_NATIVE_DOUBLE, cp_pts, ierr)
            if (ierr /= 0) stop "Error: read_grid_hdf5 -- h5dread_f"


            !
            !  Accumulate pts into a single points_t matrix to initialize domain
            !
            allocate(meshdata(iconn)%points(npts), stat=ierr)
            if (ierr /= 0) call AllocationError
                

            do ipt = 1,rank_one_dims(1)
                call meshdata(iconn)%points(ipt)%set(real(xpts(ipt),rk),real(ypts(ipt),rk),real(zpts(ipt),rk))
            end do





!            !
!            ! Open Elements connectivity dataset
!            !
!            call h5dopen_f(gid, "Elements", did_e, ierr, H5P_DEFAULT_F)
!
!            !
!            !  Get the dataspace id and dimensions
!            !
!            call h5dget_space_f(did_e, sid, ierr)
!            call h5sget_simple_extent_dims_f(sid, rank_two_dims, maxdims, ierr)
!            if (allocated(connectivity)) deallocate(connectivity)
!            allocate(connectivity(rank_two_dims(1),rank_two_dims(2)))
!
!
!            !
!            ! Read connectivity
!            ! 
!            cp_conn = c_loc(connectivity(1,1))
!            call h5dread_f(did_e, H5T_NATIVE_INTEGER, cp_conn, ierr)
!            if (ierr /= 0) stop "Error: read_grid_hdf5 -- h5dread_f"
!            meshdata(idom)%connectivity = connectivity

            !meshdata(idom)%connectivity = partition%connectivities(iconn)%data
            nelements = partition%connectivities(iconn)%get_nelements()
            nnodes    = partition%connectivities(iconn)%get_nnodes()
            call meshdata(iconn)%connectivity%init(nelements,nnodes)
            meshdata(iconn)%connectivity%data = partition%connectivities(iconn)%data


            !
            ! Read equation set attribute
            !
            meshdata(iconn)%eqnset = eqnset(idom)


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


        end do  ! iconn


        !
        !  Close file and Fortran interface
        !
        call h5fclose_f(fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"read_grid_hdf: error closing file.")
        call h5close_f(ierr)

    end subroutine read_grid_partition_hdf
    !***************************************************************************************************************
  
  
  
   
   
   
   
   
   
  
  
  
   
   
   
   
  
  
   
   
       !> Read HDF5 variable
       !!
       !!  Opens a given hdf5 file. Loads the EquationSet and solution order and calls solution initialization
       !!  procedure for each domain. Searches for the given variable and time instance. If it finds it, load to a
       !!
       !!  Note: Convention is that all floating-point data is double precision format. Conversion to working-precision
       !!        should happen after reading the data from the HDF file.
       !!
       !!  @author Nathan A. Wukie
       !!  @date   2/3/2016
       !!
       !!  @param[in]      fid         HDF5 file identifier.
       !!  @param[in]      varstring   Character string of the variable name to be read.
       !!  @param[in]      itime       Integer of the time instance for the current variable to be read.
       !!  @param[in]      dname       Character string of the domain to be read from.
       !!  @param[inout]   data        ChiDG data containing domains. Already allocated.
       !!
       !------------------------------------------------------------------------------------------------------------------
       subroutine read_variable_hdf(fid,varstring,itime,dname,data)
           use ISO_C_BINDING
           integer(HID_T),     intent(in)      :: fid
           character(*),       intent(in)      :: varstring
           integer(ik),        intent(in)      :: itime
           character(*),       intent(in)      :: dname
           type(chidg_data_t), intent(inout)   :: data
   
   
           integer(HID_T)   :: gid, sid, vid           ! Identifiers
           integer(HSIZE_T) :: maxdims(3)              ! Dataspace dimensions
           integer(HSIZE_T) :: dims(3)
  
           integer, dimension(1)           :: ibuf
  
           character(100)                  :: cbuf
           character(100)                  :: var_gqp
   
           real(rdouble), allocatable, target   :: var(:,:,:)
           real(rdouble), allocatable           :: bufferterms(:)
           type(c_ptr)                     :: cp_var
  
           integer(ik)                     :: spacedim, ielem_g
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
  
   
   
           spacedim = data%mesh(idom)%spacedim
   
           if ( spacedim == THREE_DIM ) then
               nterms_s = nterms_1d*nterms_1d*nterms_1d
           else if ( spacedim == TWO_DIM ) then
               nterms_s = nterms_1d*nterms_1d
           end if
   
   
           
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
           !  Loop through elements and set 'variable' values
           !
           do ielem = 1,data%mesh(idom)%nelem
               !
               ! Get number of terms initialized for the current element
               !
               nterms_ielem = data%sdata%q%dom(idom)%vecs(ielem)%nterms()
               ielem_g      = data%mesh(idom)%elems(ielem)%ielement_g


               !
               ! Allocate bufferterm storage that will be used to set variable data
               !
               if (allocated(bufferterms)) deallocate(bufferterms)
               allocate(bufferterms(nterms_ielem), stat=ierr)
               if (ierr /= 0) call AllocationError
               bufferterms = ZERO


               !
               ! Check for reading lower, higher, or same-order solution
               !
               if ( nterms_s < nterms_ielem ) then
               bufferterms(1:nterms_s) = var(1:nterms_s, ielem_g, itime)             ! Reading a lower order solution
               else if ( nterms_s > nterms_ielem ) then
               bufferterms(1:nterms_ielem) = var(1:nterms_ielem, ielem_g, itime)     ! Reading a higher-order solution
               else
               bufferterms(1:nterms_ielem) = var(1:nterms_ielem, ielem_g, itime)     ! Reading a solution of same order
               end if



               call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar,real(bufferterms,rk))


           end do


   
   
           !
           ! Close variable dataset, domain/variable group.
           !
           call h5dclose_f(vid,ierr)       ! Close the variable dataset
           call h5gclose_f(gid,ierr)       ! Close the Domain/Variable group
  
   
       end subroutine read_variable_hdf
       !*****************************************************************************************************************
   
   
   
   
   
   
  
   
   








    !>  Write HDF5 variable
    !!
    !!  Opens a given hdf5 file. Loads the EquationSet and solution order and calls solution initialization
    !!  procedure for each domain. Searches for the given variable and time instance. If it finds it, load to a
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      fid         HDF5 file identifier.
    !!  @param[in]      varstring   Character string of the variable name to be read.
    !!  @param[in]      itime       Integer of the time instance for the current variable to be read.
    !!  @param[in]      dname       Character string of the domain name to be read from.
    !!  @param[inout]   data        chidg_data_t instance containing grid and solution.
    !!
    !-----------------------------------------------------------------------------------------------------------------
    subroutine write_variable_hdf(fid,varstring,itime,dname,data)
        integer(HID_T),     intent(in)      :: fid
        character(*),       intent(in)      :: varstring
        integer(ik),        intent(in)      :: itime
        character(*),       intent(in)      :: dname
        type(chidg_data_t), intent(inout)   :: data


        integer(HID_T)   :: gid, sid, did, crp_list         ! Identifiers
        integer(HID_T)   :: grid_id, sid_e, did_e           ! Identifiers
        integer(HID_T)   :: memspace, filespace             ! Identifiers
        integer(HSIZE_T) :: edims(2), maxdims(3), adim      ! Dataspace dimensions
        integer(HSIZE_T) :: dims(3), dimsm(3)               ! Dataspace dimensions
        integer(HSIZE_T) :: dimsc(3)                        ! Chunk size for extendible data sets
        integer(HSIZE_T) :: start(3), count(3)
        type(H5O_INFO_T) :: info                            ! Object info type

        integer                             :: ndims
        integer, dimension(1)               :: ibuf
        character(100)                      :: cbuf
        character(100)                      :: var_grp
        character(100)                      :: ctime

        real(rdouble), allocatable, target  :: var(:,:,:)
        type(c_ptr)                         :: cp_var

        integer(ik)                         :: nmembers,    type,       ierr,       ndomains,   igrp,   &
                                               npts,        order,              &
                                               ivar,        ielem,      idom,       nelem_g,    ielement_g
        logical                             :: FileExists, VariablesExists, DataExists, ElementsEqual
        logical                             :: exists



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
            if (ierr /= 0) call chidg_signal(FATAL,"write_variable_hdf: Domain/Variables group did not open properly")
        else
            ! If 'Variables group does not exist, then create one.
            call h5gcreate_f(fid, trim(dname)//"/Variables", gid, ierr)
        end if



        !
        ! Get total number of elements in the domain from the grid file
        !
        call h5gopen_f(fid, trim(dname)//"/Grid", grid_id, ierr, H5P_DEFAULT_F)
        if (ierr /= 0) call chidg_signal_one(FATAL,"write_variable_hdf: Domagin/Grid group did not open properly.", trim(dname)//'/Grid')
        call h5dopen_f(grid_id, "Elements", did_e, ierr, H5P_DEFAULT_F)
        call h5dget_space_f(did_e, sid_e, ierr)
        call h5sget_simple_extent_dims_f(sid_e, edims, maxdims, ierr)
        nelem_g  = edims(1)
        
        call h5dclose_f(did_e,ierr)
        call h5sclose_f(sid_e,ierr)
        call h5gclose_f(grid_id,ierr)





        !
        ! Set dimensions of dataspace to write
        !
        idom = data%get_domain_index(dname)
        ndims = 3

        dims(1) = data%mesh(idom)%nterms_s
        dims(2) = nelem_g
        dims(3) = itime                     ! TODO: Should probably better inform the dataspace dimension here. Probably set mesh_t%ntime
        maxdims(1) = H5S_UNLIMITED_F
        maxdims(2) = H5S_UNLIMITED_F
        maxdims(3) = H5S_UNLIMITED_F




        !
        ! Open the Variable dataset, check if Variable dataset already exists
        !
        call h5lexists_f(gid, trim(varstring), exists, ierr)



        !
        ! Modify dataset creation properties, i.e. enable chunking in order to append dataspace, if needed.
        !
        dimsc = [1, nelem_g, 1]  ! Chunk size

        call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, ierr)
        if (ierr /= 0) call chidg_signal(FATAL, "write_variable_hdf: h5pcreate_f error enabling chunking")

        call h5pset_chunk_f(crp_list, ndims, dimsc, ierr)
        if (ierr /= 0) call chidg_signal(FATAL, "write_variable_hdf: h5pset_chunk_f error setting chunk properties")



        !
        ! Reset dataspace size if necessary
        !
        if (exists) then
            ! Open the existing dataset
            call h5dopen_f(gid, trim(varstring), did, ierr, H5P_DEFAULT_F)
            if (ierr /= 0) call chidg_signal(FATAL,"write_variable_hdf: variable does not exist or was not opened correctly")


            ! Extend dataset if necessary
            call h5dset_extent_f(did, dims, ierr)
            if (ierr /= 0) call chidg_signal(FATAL, "write_variable_hdf: h5dset_extent_f")


            ! Update existing dataspace ID since it may have been expanded
            call h5dget_space_f(did, sid, ierr)
            if (ierr /= 0) call chidg_signal(FATAL, "write_variable_hdf: h5dget_space_f")

        else
            ! Create a new dataspace
            call h5screate_simple_f(ndims,dims,sid,ierr,maxdims)
            if (ierr /= 0) call chidg_signal(FATAL,"write_variable_hdf: h5screate_simple_f")


            ! Create a new dataset
            call h5dcreate_f(gid, trim(varstring), H5T_NATIVE_DOUBLE, sid, did, ierr, crp_list)
            if (ierr /= 0) call chidg_signal(FATAL,"write_variable_hdf: h5dcreate_f")
        end if



        !
        ! Get variable integer index from variable character string
        !
        ivar = data%eqnset(idom)%item%prop%get_eqn_index(varstring)



        !
        ! Assemble variable buffer matrix that gets written to file
        !
        allocate(var(dims(1),1,1))


        do ielem = 1,data%mesh(idom)%nelem

            ! get domain-global element index
            ielement_g = data%mesh(idom)%elems(ielem)%ielement_g

            start    = [1-1,ielement_g-1,itime-1]   ! actually offset, so 0-based
            count(1) = dims(1)
            count(2) = 1
            count(3) = 1

            ! Select subset of dataspace - sid
            call h5sselect_hyperslab_f(sid, H5S_SELECT_SET_F, start, count, ierr)


            ! Create a memory dataspace
            dimsm(1) = size(var,1)
            dimsm(2) = size(var,2)
            dimsm(3) = size(var,3)
            call h5screate_simple_f(ndims,dimsm,memspace,ierr)


            var(:,1,1) = real(data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar),rdouble)
            cp_var = c_loc(var(1,1,1))

            call h5dwrite_f(did, H5T_NATIVE_DOUBLE, cp_var, ierr, memspace, sid)


            call h5sclose_f(memspace,ierr)


        end do



        !
        ! Update 'contains_solution' attribute
        !
        call h5ltset_attribute_string_f(fid, "/", 'contains_solution', 'Yes', ierr)


        call h5pclose_f(crp_list, ierr) ! Close dataset creation property
        call h5dclose_f(did,ierr)       ! Close Variable datasets
        call h5sclose_f(sid,ierr)       ! Close Variable dataspaces
        call h5gclose_f(gid,ierr)       ! Close Domain/Variable group



    end subroutine write_variable_hdf
    !************************************************************************************************************



















    !> Read solution modes from HDF file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      filename    Character string of the file to be read from
    !!  @param[inout]   data        chidg_data_t that will accept the solution modes
    !!
    !-------------------------------------------------------------------------------------------------------------
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
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, fid, ierr)
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

                ! Get variable character string
                cvar = trim(data%eqnset(idom)%item%prop%eqns(ieqn)%name)

                ! Read variable
                call read_variable_hdf(fid,cvar,time,trim(dname),data)

            end do ! ieqn

        end do ! idom



        ! Close HDF5 file and Fortran interface
        call h5fclose_f(fid, ierr)
        call h5close_f(ierr)

    end subroutine read_solution_hdf
    !*************************************************************************************************************




















    !> Write solution modes to HDF file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      filename    Character string of the file to be written to
    !!  @param[inout]   data        chidg_data_t containing solution to be written
    !!
    !!  @TODO   Allow for creation of a new solution file. Currently, incoming filename needs to exist already.
    !!
    !-------------------------------------------------------------------------------------------------------
    subroutine write_solution_hdf(filename,data)
        use ISO_C_BINDING
        character(*),       intent(in)      :: filename
        type(chidg_data_t), intent(inout)   :: data


        integer(HID_T)                  :: fid
        integer(HSIZE_T)                :: adim
        integer(ik)                     :: idom, ndomains
        integer(ik)                     :: ieqn, neqns, spacedim
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
            spacedim = data%mesh(idom)%spacedim

            if ( spacedim == THREE_DIM ) then

                do while ( order_s*order_s*order_s /= data%mesh(idom)%nterms_s )
                   order_s = order_s + 1 
                end do
                order_s = order_s - 1 ! to be consistent with he definition of 'Order of the polynomial'

            else if ( spacedim == TWO_DIM ) then
                do while ( order_s*order_s /= data%mesh(idom)%nterms_s )
                   order_s = order_s + 1 
                end do
                order_s = order_s - 1 ! to be consistent with he definition of 'Order of the polynomial'

            end if

            call h5ltset_attribute_int_f(fid, trim(dname), 'order_solution', [order_s], adim, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"write_solution_partition_hdf5 - h5ltset_attribute_int_f")

            call h5ltset_attribute_string_f(fid, trim(dname), 'eqnset', trim(data%eqnset(idom)%item%name), ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"write_solution_partition_hdf5 - h5ltset_attribute_int_f")




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
        if (ierr /= 0) call chidg_signal(FATAL,"write_solution_hdf: error closing file.")
        call h5close_f(ierr)            ! Close HDF5 Fortran interface
        if (ierr /= 0) call chidg_signal(FATAL,"write_solution_hdf: h5close.")

    end subroutine write_solution_hdf
    !*********************************************************************************************************












    !>  Read boundary conditions from HDF5 file in ChiDG format and return data in bcdata_t container.
    !!  The calling procedure can then use the returned bcdata_t to initialize boundary conditions.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!  @param[in]  filename        String of the HDF5 file to be read.
    !!  @param[inout]   bcdata(:)   Array of bcdata_t instances, one for each domain. These will be returned with
    !!                              data about the boundary conditions that can be used for initialization.
    !!
    !---------------------------------------------------------------------------------------------------------
    subroutine read_boundaryconditions_hdf(filename, bcdata)
        character(*),   intent(in)                  :: filename
        type(bcdata_t), intent(inout), allocatable  :: bcdata(:)

        integer(HID_T)                          :: fid, bcgroup, bcface, bcprop, faces_did, faces_sid
        integer(HSIZE_T)                        :: adim
        integer(HSIZE_T)                        :: rank_two_dims(2), maxdims(2)                    ! Dataspace dimensions
        logical                                 :: FileExists

        class(bc_t),            allocatable     :: bc
        character(len=1024)                     :: bcname, pname, oname, fname
        real(rk)                                :: ovalue
        real(rdouble),   dimension(1)           :: rbuf
        character(len=10)                       :: faces(NFACES)
        character(1024)                         :: gname
        integer                                 :: nmembers, type, ierr, ndomains, igrp, &
                                                   idom, iface, iopt, noptions, nprop, iprop, &
                                                   npts_face, nbcfaces, ibc_face
        integer, dimension(1)                   :: buf
        integer,            allocatable, target :: bc_connectivity(:,:)
        type(c_ptr)                             :: cp_conn




        faces = ["  XI_MIN","  XI_MAX"," ETA_MIN"," ETA_MAX","ZETA_MIN","ZETA_MAX"]



        !
        !  Check file exists
        !
        inquire(file=filename, exist=FileExists)
        if (.not. FileExists) then
            call chidg_signal(FATAL,'read_boundaryconditions_hdf: Could not find grid file')
        end if


        !
        ! Open hdf interface
        !
        call h5open_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'read_boundaryconditions_hdf - h5open_f: HDF5 Fortran interface had an error during initialization')


        !
        ! Open file
        !
        !call h5fopen_f(filename, H5F_ACC_RDONLY_F, fid, ierr)
        call h5fopen_f(filename, H5F_ACC_RDWR_F, fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'read_boundaryconditions_hdf - h5fopen_f: There was an error opening the grid file.')



        !
        !  Get number of domains from attribute 'ndomains' in file root
        !
        ndomains = get_ndomains_hdf(fid)


        !
        !  Allocate number of domains
        !
        if (ndomains == 0) call chidg_signal(FATAL,'read_boundaryconditions_hdf: No Domains were found in the file')
        allocate(bcdata(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Loop domains
        !
        !  Get number of groups in the file root
        call h5gn_members_f(fid, "/", nmembers, ierr)


        !  Loop through groups and read domains
        idom = 1
        do igrp = 0,nmembers-1
            call h5gget_obj_info_idx_f(fid,"/", igrp, gname, type, ierr)

            if (gname(1:2) == 'D_') then

                !
                ! Set domain name.
                !
                bcdata(idom)%domain_ = gname


                !
                ! Open the Domain/BoundaryConditions group
                !
                call h5gopen_f(fid, trim(gname)//"/BoundaryConditions", bcgroup, ierr, H5P_DEFAULT_F)
                if (ierr /= 0) stop "Error: read_boundaryconditions_hdf -- h5gopen_f: Domain/BoundaryConditions group did not open properly"



                !
                ! Allocation bcs for current domain
                !
                allocate( bcdata(idom)%bcs(NFACES), bcdata(idom)%bc_connectivity(NFACES), stat=ierr )
                if (ierr /= 0) call AllocationError


                !
                ! Loop faces and get boundary condition for each
                !
                ! TODO: should probably turn this into a loop over bcs instead of faces.
                do iface = 1,NFACES


                    !
                    ! Open face boundary condition group
                    !
                    call h5gopen_f(bcgroup, trim(adjustl(faces(iface))), bcface, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_hdf: error opening boundary face group")
        

                    !
                    ! Get face associated with the 
                    !
                    ! TODO: WARNING, should replace with XI_MIN, XI_MAX, etc. somehow.
                    call h5dopen_f(bcface, "Faces", faces_did, ierr, H5P_DEFAULT_F)


                    !
                    !  Get the dataspace id and dimensions
                    !
                    call h5dget_space_f(faces_did, faces_sid, ierr)
                    call h5sget_simple_extent_dims_f(faces_sid, rank_two_dims, maxdims, ierr)
                    nbcfaces  = rank_two_dims(1)
                    npts_face = rank_two_dims(2)



                    !
                    ! Read boundary condition connectivity
                    !
                    if ( allocated(bc_connectivity) ) deallocate(bc_connectivity)
                    allocate(bc_connectivity(nbcfaces,npts_face),stat=ierr)
                    if (ierr /= 0) call AllocationError
                    cp_conn = c_loc(bc_connectivity(1,1))
                    call h5dread_f(faces_did, H5T_NATIVE_INTEGER, cp_conn, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_hdf5 -- h5dread_f")

                    call h5dclose_f(faces_did,ierr)
                    call h5sclose_f(faces_sid,ierr)

                    !
                    ! Store boundary condition connectivity
                    !
                    !bcdata(idom)%bc_connectivity(iface)%data = bc_connectivity
                    call bcdata(idom)%bc_connectivity(iface)%init(nbcfaces)
                    do ibc_face = 1,nbcfaces
                        bcdata(idom)%bc_connectivity(iface)%data(ibc_face)%data = bc_connectivity(ibc_face,:)
                    end do




                    !
                    ! Get boundary condition name string
                    !
                    call h5ltget_attribute_string_f(bcface, ".", "bc_name", bcname, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_hdf: error reading boundary condition name, bc_name")


                    if ( trim(bcname) == 'empty' ) then

                        cycle

                    else

                        !
                        ! Create boundary condition from face data
                        !
                        call create_bc(trim(bcname),bc)



                        !
                        ! Get available boundary condition properties to search for.
                        !
                        nprop = bc%get_nproperties()

                        
                        !
                        ! Loop through properties
                        !
                        do iprop = 1,nprop


                            !
                            ! Get property name
                            !
                            pname = bc%get_property_name(iprop)


                            !
                            ! Open property in HDF file
                            !
                            call h5gopen_f(bcface, "BCP_"//trim(pname), bcprop, ierr)
                            if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions: error opening bcproperty group.")


                            !
                            ! Read the function name set for the property.
                            !
                            call h5ltget_attribute_string_f(bcprop, ".", "function", fname, ierr)
                            if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions: error getting function name.")

                            
                            !
                            ! Set/Create the function for the current property
                            !
                            call bc%set_fcn(trim(pname), trim(fname))

                            
                            !
                            ! Get number of options for the function
                            !
                            noptions = bc%get_noptions(iprop)


                            !
                            ! Get each option value
                            !
                            do iopt = 1,noptions

                                !
                                ! Get option name
                                !
                                oname = bc%get_option_key(iprop,iopt)


                                !
                                ! Get option value from file
                                !
                                adim = 1
                                call h5ltget_attribute_double_f(bcprop, ".", trim(oname), rbuf, ierr)
                                if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions: error getting option value")
                                ovalue = real(rbuf(1),rk)


                                !
                                ! Set boundary condition option
                                !
                                call bc%set_fcn_option(trim(pname), trim(oname), ovalue)


                            end do ! iopt


                            !
                            ! Close current property group
                            !
                            call h5gclose_f(bcprop,ierr)
                            if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_hdf: h5gclose")



                        end do !iprop



                        !
                        ! Boundary condition is defined for the current face, save to bcdata
                        !
                        allocate(bcdata(idom)%bcs(iface)%bc, source = bc, stat=ierr)
                        if (ierr /= 0) call AllocationError


                    end if ! bcname == empty


                    !
                    ! Close face boundary condition group
                    !
                    call h5gclose_f(bcface, ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_hdf: h5gclose")


                end do ! iface



                idom = idom + 1



                !
                ! Close BoundaryCondition group
                !
                call h5gclose_f(bcgroup, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_hdf: h5gclose")


            end if !gname, idom


        end do  ! igrp


        !
        !  Close file and Fortran interface
        !
        call h5fclose_f(fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_hdf: error closing file.")
        call h5close_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_hdf: error closing hdf interface.")


    end subroutine read_boundaryconditions_hdf
    !*********************************************************************************************************




















    !>  Read boundary conditions from HDF5 file in ChiDG format and return data in bcdata_t container.
    !!  The calling procedure can then use the returned bcdata_t to initialize boundary conditions.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/9/2016
    !!
    !!  @param[in]      filename    String of the HDF5 file to be read.
    !!  @param[inout]   bcdata(:)   Array of bcdata_t instances, one for each domain. These will be returned with
    !!                              data about the boundary conditions that can be used for initialization.
    !!  @param[in]      partition   Partition information to only read boundary conditions for the domains in the partition
    !!
    !---------------------------------------------------------------------------------------------------------
    subroutine read_boundaryconditions_partition_hdf(filename, bcdata, partition)
        character(*),       intent(in)                  :: filename
        type(bcdata_t),     intent(inout), allocatable  :: bcdata(:)
        type(partition_t),  intent(in)                  :: partition

        integer(HID_T)                          :: fid, bcgroup, bcface, bcprop, faces_did, faces_sid
        integer(HSIZE_T)                        :: adim
        integer(HSIZE_T)                        :: rank_two_dims(2), maxdims(2)                    ! Dataspace dimensions
        logical                                 :: FileExists

        class(bc_t),            allocatable     :: bc
        character(len=1024)                     :: bcname, pname, oname, fname
        real(rk)                                :: ovalue
        real(rdouble),   dimension(1)           :: rbuf
        character(len=10)                       :: faces(NFACES)
        character(1024)                         :: gname
        integer                                 :: nmembers, type, ierr, ndomains, igrp, &
                                                   idom, iface, iopt, noptions, nprop, iprop, &
                                                   npts_face, nbcfaces, iconn, nconn, ibc_face
        integer, dimension(1)                   :: buf
        integer,            allocatable, target :: bc_connectivity(:,:)
        type(c_ptr)                             :: cp_conn




        faces = ["  XI_MIN","  XI_MAX"," ETA_MIN"," ETA_MAX","ZETA_MIN","ZETA_MAX"]



        !
        !  Check file exists
        !
        inquire(file=filename, exist=FileExists)
        if (.not. FileExists) then
            call chidg_signal(FATAL,'read_boundaryconditions_partition_hdf: Could not find grid file')
        end if


        !
        ! Open hdf interface
        !
        call h5open_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'read_boundaryconditions_partition_hdf - h5open_f: HDF5 Fortran interface had an error during initialization')


        !
        ! Open file
        !
        !call h5fopen_f(filename, H5F_ACC_RDONLY_F, fid, ierr)
        call h5fopen_f(filename, H5F_ACC_RDWR_F, fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'read_boundaryconditions_partition_hdf - h5fopen_f: There was an error opening the grid file.')




        !
        !  Allocate for number of connectivities in the partition
        !
        nconn = size(partition%connectivities)
        allocate(bcdata(nconn), stat=ierr)
        if (ierr /= 0) call AllocationError




        !
        !  Loop through connectivities and read boundary conditions
        !
        do iconn = 1,nconn

            ! Get domain index of current connectivity
            idom = partition%connectivities(iconn)%get_domain_index()

            ! Get name of current domain
            gname = get_domain_name_hdf(fid,idom)

            ! Set domain name.
            bcdata(iconn)%domain_ = gname


            !
            ! Open the Domain/BoundaryConditions group
            !
            call h5gopen_f(fid, trim(gname)//"/BoundaryConditions", bcgroup, ierr, H5P_DEFAULT_F)
            if (ierr /= 0) stop "Error: read_boundaryconditions_partition_hdf -- h5gopen_f: Domain/BoundaryConditions group did not open properly"



            !
            ! Allocation bcs for current domain
            !
            allocate( bcdata(iconn)%bcs(NFACES), bcdata(iconn)%bc_connectivity(NFACES), stat=ierr )
            if (ierr /= 0) call AllocationError


            !
            ! Loop faces and get boundary condition for each
            !
            ! TODO: should probably turn this into a loop over bcs instead of faces.
            do iface = 1,NFACES


                !
                ! Open face boundary condition group
                !
                call h5gopen_f(bcgroup, trim(adjustl(faces(iface))), bcface, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_partition_hdf: error opening boundary face group")
    

                !
                ! Get boundary condition name string
                !
                call h5ltget_attribute_string_f(bcface, ".", "bc_name", bcname, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_partition_hdf: error reading boundary condition name, bc_name")



                !
                ! Get face associated with the 
                !
                ! TODO: WARNING, should replace with XI_MIN, XI_MAX, etc. somehow. Maybe not...
                call h5dopen_f(bcface, "Faces", faces_did, ierr, H5P_DEFAULT_F)


                !
                !  Get the dataspace id and dimensions
                !
                call h5dget_space_f(faces_did, faces_sid, ierr)
                call h5sget_simple_extent_dims_f(faces_sid, rank_two_dims, maxdims, ierr)
                nbcfaces  = rank_two_dims(1)
                npts_face = rank_two_dims(2)



                !
                ! Read boundary condition connectivity
                !
                if ( allocated(bc_connectivity) ) deallocate(bc_connectivity)
                allocate(bc_connectivity(nbcfaces,npts_face),stat=ierr)
                if (ierr /= 0) call AllocationError
                cp_conn = c_loc(bc_connectivity(1,1))
                call h5dread_f(faces_did, H5T_NATIVE_INTEGER, cp_conn, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_partition_hdf5 -- h5dread_f")


                call h5dclose_f(faces_did,ierr)
                call h5sclose_f(faces_sid,ierr)


                !
                ! Store boundary condition connectivity
                !
                call bcdata(iconn)%bc_connectivity(iface)%init(nbcfaces)
                do ibc_face = 1,nbcfaces
                    bcdata(iconn)%bc_connectivity(iface)%data(ibc_face)%data = bc_connectivity(ibc_face,:)
                end do








!                if ( trim(bcname) == 'empty' ) then
!
!                    cycle
!
!                else


                    !
                    ! Create boundary condition from face data
                    !
                    call create_bc(trim(bcname),bc)



                    !
                    ! Get available boundary condition properties to search for.
                    !
                    nprop = bc%get_nproperties()

                    
                    !
                    ! Loop through properties
                    !
                    do iprop = 1,nprop


                        !
                        ! Get property name
                        !
                        pname = bc%get_property_name(iprop)


                        !
                        ! Open property in HDF file
                        !
                        call h5gopen_f(bcface, "BCP_"//trim(pname), bcprop, ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_partition_hdf: error opening bcproperty group.")


                        !
                        ! Read the function name set for the property.
                        !
                        call h5ltget_attribute_string_f(bcprop, ".", "function", fname, ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_partition_hdf: error getting function name.")

                        
                        !
                        ! Set/Create the function for the current property
                        !
                        call bc%set_fcn(trim(pname), trim(fname))

                        
                        !
                        ! Get number of options for the function
                        !
                        noptions = bc%get_noptions(iprop)


                        !
                        ! Get each option value
                        !
                        do iopt = 1,noptions

                            !
                            ! Get option name
                            !
                            oname = bc%get_option_key(iprop,iopt)


                            !
                            ! Get option value from file
                            !
                            adim = 1
                            call h5ltget_attribute_double_f(bcprop, ".", trim(oname), rbuf, ierr)
                            if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_partition_hdf: error getting option value")
                            ovalue = real(rbuf(1),rk)


                            !
                            ! Set boundary condition option
                            !
                            call bc%set_fcn_option(trim(pname), trim(oname), ovalue)


                        end do ! iopt


                        !
                        ! Close current property group
                        !
                        call h5gclose_f(bcprop,ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_partition_hdf: h5gclose")



                    end do !iprop



                    !
                    ! Boundary condition is defined for the current face, save to bcdata
                    !
                    allocate(bcdata(iconn)%bcs(iface)%bc, source=bc, stat=ierr)
                    if (ierr /= 0) call AllocationError


!                end if ! bcname == empty


                !
                ! Close face boundary condition group
                !
                call h5gclose_f(bcface, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_partition_hdf: h5gclose")


            end do ! iface


            ! Close BoundaryCondition group
            call h5gclose_f(bcgroup, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_partition_hdf: h5gclose")


        end do  ! iconn


        !
        !  Close file and Fortran interface
        !
        call h5fclose_f(fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_partition_hdf: error closing file.")
        call h5close_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"read_boundaryconditions_partition_hdf: error closing hdf interface.")


    end subroutine read_boundaryconditions_partition_hdf
    !*********************************************************************************************************

























    !>  This reads an HDF ChiDG grid file and returns an array of connectivities, one for each domain.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/9/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------------
    subroutine read_connectivity_hdf(filename, connectivities)
        use mod_io, only: nterms_s
        character(*),                               intent(in)      :: filename
        type(domain_connectivity_t), allocatable,   intent(inout)   :: connectivities(:)

        integer(HID_T)   :: fid, gid, sid, did_x, did_e               ! Identifiers
        integer(HSIZE_T) :: rank_one_dims(1), rank_two_dims(2), dims(3), maxdims(3) ! Dataspace dimensions

        integer,                     allocatable, target    :: connectivity(:,:)
        type(c_ptr)                                         :: cp_conn

        integer(ik),            allocatable     :: domain_indices(:)
        character(len=1024),    allocatable     :: dnames(:), eqnset(:)
        character(1024)                         :: gname
        integer                                 :: nmembers, type, ierr, ndomains, igrp,    &
                                                   idom, idomain, nelements, ielem, nnodes, mapping
        logical                                 :: FileExists



        !
        !  Check file exists
        !
        inquire(file=filename, exist=FileExists)
        if (.not. FileExists) then
            call chidg_signal_one(FATAL,'read_connectivity_hdf5: Could not find grid file',filename)
        end if


        !
        !  Initialize Fortran interface.
        !
        call h5open_f(ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'read_connectivity_hdf5 - h5open_f: HDF5 Fortran interface had an error during initialization')



        !
        !  Open input file using default properties.
        !
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'read_connectivity_hdf5 - h5fopen_f: There was an error opening the grid file.')



        !
        !  Get number of domains from attribute 'ndomains' in file root
        !
        ndomains = get_ndomains_hdf(fid)


        !
        !  Allocate number of domains
        !
        if (ndomains == 0) call chidg_signal(FATAL,'read_connectivity_hdf5: No Domains were found in the file')
        allocate(connectivities(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError




        !
        !  Get number of groups in the file root
        !
        call h5gn_members_f(fid, "/", nmembers, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"read_connectivity_hdf5: error getting number of groups in file root.")



        !
        ! Get equationset strings.
        !
        dnames         = get_domain_names_hdf(fid)
        domain_indices = get_domain_indices_hdf(fid)
        eqnset         = get_eqnset_hdf(fid,dnames)



        !
        !  Loop through groups and read domain connectivities
        !
        idom = 1
        do idom = 1,size(dnames)

                gname = dnames(idom)

                !
                ! Open the Domain/Grid group
                !
                call h5gopen_f(fid, trim(gname)//"/Grid", gid, ierr, H5P_DEFAULT_F)
                if (ierr /= 0) call chidg_signal_one(FATAL,"read_connectivity_hdf: Domagin/Grid group did not open properly.", trim(gname)//'/Grid')


                !
                ! Get number of nodes in the domain
                !
                call h5dopen_f(gid, "CoordinateX", did_x, ierr, H5P_DEFAULT_F)


                !
                !  Get the dataspace id and dimensions
                !
                call h5dget_space_f(did_x, sid, ierr)
                call h5sget_simple_extent_dims_f(sid, rank_one_dims, maxdims, ierr)
                nnodes = rank_one_dims(1)


                !
                ! Open Elements connectivity dataset
                !
                call h5dopen_f(gid, "Elements", did_e, ierr, H5P_DEFAULT_F)

                !
                !  Get the dataspace id and dimensions
                !
                call h5dget_space_f(did_e, sid, ierr)
                call h5sget_simple_extent_dims_f(sid, rank_two_dims, maxdims, ierr)
                if (allocated(connectivity)) deallocate(connectivity)
                allocate(connectivity(rank_two_dims(1),rank_two_dims(2)))


                !
                ! Read connectivity
                ! 
                cp_conn = c_loc(connectivity(1,1))
                call h5dread_f(did_e, H5T_NATIVE_INTEGER, cp_conn, ierr)
                if (ierr /= 0) stop "Error: read_connectivity_hdf5 -- h5dread_f"




                ! Initialize domain connectivity structure
                idomain   = domain_indices(idom)    ! prob don't need this
                nelements = size(connectivity,1)
                call connectivities(idom)%init(nelements, nnodes)


                !connectivities(idom)%data = connectivity
                do ielem = 1,nelements
                    mapping = connectivity(ielem,3)
                    call connectivities(idom)%data(ielem)%init(mapping)
                    connectivities(idom)%data(ielem)%data = connectivity(ielem,:)
                    call connectivities(idom)%data(ielem)%set_element_partition(NO_PROC)
                end do


                !
                ! Close the Elements datasets
                !
                call h5dclose_f(did_e, ierr)
                call h5dclose_f(did_x, ierr)


                !
                ! Close the dataspace id
                !
                call h5sclose_f(sid,ierr)


                !
                ! Close the Domain/Grid group
                !
                call h5gclose_f(gid,ierr)


        end do  ! igrp


        !
        !  Close file and Fortran interface
        !
        call h5fclose_f(fid, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"read_connectivity_hdf: error closing file.")
        call h5close_f(ierr)

    end subroutine read_connectivity_hdf
    !********************************************************************************************************
























end module mod_hdfio
