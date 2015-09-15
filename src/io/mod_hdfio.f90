module mod_hdfio
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use type_point,     only: point_t
    use type_domain,    only: domain_t
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
    !!
    !!  @param[in]      filename    Character string of the file to be read
    !!  @param[inout]   domains     Allocatable array of domains. Allocated in this routine.
    !------------------------------------------------------------------------------------------
    subroutine read_grid_hdf(filename, domains)
        character(*),                   intent(in)    :: filename
        type(domain_t), allocatable,    intent(inout) :: domains(:)

        integer(HID_T)   :: fid, gid, sid, did_x, did_y, did_z      ! Identifiers
        integer(HSIZE_T) :: dims(3), maxdims(3)                     ! Dataspace dimensions

        type(c_ptr)                                     :: pts
        type(point_t), allocatable                      :: points(:,:,:)
        real(rk), dimension(:,:,:), allocatable, target :: xpts, ypts, zpts
        type(c_ptr)                                     :: cp_pts

        character(10)                           :: gname
        integer                                 :: nmembers, type, ierr, ndomains, igrp,    &
                                                   npts, izeta, ieta, ixi, idom, nterms_1d, &
                                                   nterms_c, mapping
        integer, dimension(1)                   :: buf
        logical                                 :: FileExists

        !  Check file exists
        inquire(file=filename, exist=FileExists)
        if (.not. FileExists) then
            call signal(FATAL,'read_grid_hdf5: Could not find grid file')
        end if


        !  Initialize Fortran interface.
        call h5open_f(ierr)
        if (ierr /= 0) call signal(FATAL,'read_grid_hdf5 - h5open_f: HDF5 Fortran interface had an error during initialization')



        !  Open input file using default properties.
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, fid, ierr)
        if (ierr /= 0) call signal(FATAL,'read_grid_hdf5 - h5fopen_f: There was an error opening the grid file.')



        !  Get number of domains from attribute 'ndomains' in file root
        call h5ltget_attribute_int_f(fid, "/", 'ndomains', buf, ierr)
        ndomains = buf(1)
        if (ierr /= 0) call signal(FATAL,'read_grid_hdf5: h5ltget_attribute_int_f had a problem getting the number of domains')



        !  Allocate number of domains
        if (ndomains == 0) then
            call signal(FATAL,'read_hdf5: No Domains were found in the file')
        else
            allocate(domains(ndomains), stat=ierr)
            if (ierr /= 0) call AllocationError
        end if


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


                !  Get number of terms
                call h5ltget_attribute_int_f(fid, "/", 'mapping', buf, ierr)
                mapping = buf(1)
                if (ierr /= 0) stop "Error: read_grid_hdf5 - h5ltget_attribute_int_f"
                nterms_1d = (mapping + 1)
                nterms_c = nterms_1d * nterms_1d * nterms_1d


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
                allocate(points(dims(1),dims(2),dims(3)), stat=ierr)
                if (ierr /= 0) call AllocationError
                    
                do izeta = 1,dims(3)
                    do ieta = 1,dims(2)
                        do ixi = 1,dims(1)
                            call points(ixi,ieta,izeta)%set(xpts(ixi,ieta,izeta),ypts(ixi,ieta,izeta),zpts(ixi,ieta,izeta))
                        end do
                    end do
                end do


                ! Call domain geometry initialization
                call domains(idom)%init_geom(nterms_c,points)


                ! Close the Coordinate datasets
                call h5dclose_f(did_x,ierr)
                call h5dclose_f(did_y,ierr)
                call h5dclose_f(did_z,ierr)


                ! Close the dataspace id
                call h5sclose_f(sid,ierr)


                ! Close the Domain/Grid group
                call h5gclose_f(gid,ierr)

                ! Deallocate points for the current domain
                deallocate(zpts,ypts,xpts,points)
                idom = idom + 1
            end if
        end do

        !  Close file
        call h5fclose_f(fid, ierr)

        !  Close Fortran interface
        call h5close_f(ierr)

    end subroutine






    !> Read HDF5 variable
    !!
    !!  Opens a given hdf5 file. Loads the EquationSet and solution order and calls solution initialization
    !!  procedure for each domain. Searches for the given variable and time instance. If it finds it, load to a
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      filename    Character string of the file to be read
    !!  @param[in]      cvar        Character string of the variable name to be read
    !!  @param[in]      time        Integer of the time instance for the current variable to be read
    !!  @param[inout]   domains     Array of domains. Already allocated
    !-----------------------------------------------------------------------------------------------------------
    subroutine read_var_hdf(filename,cvar,time,domains)
        use ISO_C_BINDING
        character(*),   intent(in)      :: filename
        character(*),   intent(in)      :: cvar
        integer(ik),    intent(in)      :: time
        type(domain_t), intent(inout)   :: domains(:)


        integer(HID_T)   :: fid, gid, sid, vid          !> Identifiers
        integer(HSIZE_T) :: dims(2), maxdims(2)         !> Dataspace dimensions

        integer, dimension(1)           :: ibuf
        character(100)                  :: cbuf, eqnstring, varstring, var_grp, ctime

!        real(rk), allocatable           :: var(:,:)
        real(rk), allocatable, target   :: var(:,:)
        type(c_ptr)                     :: cp_var

        character(10)                   :: gname
        integer                         :: nmembers,    type,   ierr,       ndomains,   igrp,   &
                                           npts,        idom,   nterms_1d,  nterms_s,   order,  &
                                           ivar,        ielem
        logical                         :: FileExists, ElementsEqual



        !>  Check file exists
        inquire(file=filename, exist=FileExists)
        if (.not. FileExists) then
            print*, "Error: read_grid_hdf5 - file not found: ", filename
            stop
        end if


        !>  Initialize Fortran interface.
        call h5open_f(ierr)
        if (ierr /= 0) stop "Error: read_var_hdf5 - h5open_f"

        !>  Open input file using default properties.
        call h5fopen_f(filename, H5F_ACC_RDWR_F, fid, ierr)
        if (ierr /= 0) stop "Error: read_var_hdf5 - h5fopen_f"

        !>  Get number of domains from attribute 'ndomains' in file root
        call h5ltget_attribute_int_f(fid, "/", 'ndomains', ibuf, ierr)
        ndomains = ibuf(1)
        if (ierr /= 0) stop "Error: read_var_hdf5 - h5ltget_attribute_int_f"
        if (ndomains /= size(domains)) stop "Error: read_var_hdf5 - number of domains in file and passed domain list to not match"



        !>  Get number of groups in the file root
        call h5gn_members_f(fid, "/", nmembers, ierr)

        !>  Loop through groups and read domains
        idom = 1
        do igrp = 0,nmembers-1
            call h5gget_obj_info_idx_f(fid,"/", igrp, gname, type, ierr)

            if (gname(1:2) == 'D_') then

                !>  Get EquationSet
                call h5ltget_attribute_string_f(fid, trim(gname), 'EquationSet', cbuf, ierr)
                eqnstring = cbuf
                if (ierr /= 0) stop "Error: read_var_hdf5 - h5ltget_attribute_string_f"


                !> Open the Domain/Variables group
                call h5gopen_f(fid, trim(gname)//"/Variables", gid, ierr, H5P_DEFAULT_F)
                if (ierr /= 0) stop "Error: h5gopen_f -- Domain/Grid group did not open properly"


                !>  Get number of terms in solution expansion
                call h5ltget_attribute_int_f(gid, "/", 'Order', ibuf, ierr)
                order = ibuf(1)
                if (ierr /= 0) stop "Error: read_var_hdf5 - h5ltget_attribute_int_f"
                nterms_1d = (order + 1)
                nterms_s = nterms_1d*nterms_1d*nterms_1d


                !> Initialize numerics if they are NOT already initialized
                if (.not. domains(idom)%numInitialized) then
                    call domains(idom)%init_sol(eqnstring,nterms_s)
                end if

                !>  Open the Variable dataset
                write(ctime, '(I0.3)') time                     !> write time as character string
                varstring = trim(cvar)//'_'//trim(ctime)        !> compose variable name as 'var_time'
                call h5dopen_f(gid, trim(varstring), vid, ierr, H5P_DEFAULT_F)
                if (ierr /= 0) stop "Error: read_var_hdf5 -- variable does not exist or was not opened correctly"

                !>  Get the dataspace id and dimensions
                call h5dget_space_f(vid, sid, ierr)
                call h5sget_simple_extent_dims_f(sid, dims, maxdims, ierr)

                !>  Read 'variable' dataset
                allocate(var(dims(1),dims(2)))                          !> Allocate variable buffer
                cp_var = c_loc(var(1,1))                                !> Get C-address for buffer
                call h5dread_f(vid, H5T_NATIVE_DOUBLE, cp_var, ierr)    !> Fortran 2003 interface
                if (ierr /= 0) stop "Error: read_var_hdf5 -- h5dread_f"


                !>  Get variable index in EquationSet
                ivar = domains(idom)%eqnset%prop%get_eqn_index(trim(cvar))


                !>  Test to make sure the number of elements in the variable group
                !!  and the current domain are conforming
                ElementsEqual = (size(domains(idom)%sdata%q%lvecs) == size(var,2))
                if (ElementsEqual) then
                    !>  Loop through elements and assign 'variable' values
                    do ielem = 1,domains(idom)%mesh%nelem
                        !domains(idom)%sdata%q%lvecs(ielem)%mat(:,ivar) = var(:,ielem)
                        call domains(idom)%sdata%q%lvecs(ielem)%setvar(ivar,var(:,ielem))
                    end do
                else
                    stop "Error: read_var_hdf5 -- number of elements in file variable and domain do not match"
                end if


                !> Close the Variable datasets
                call h5dclose_f(vid,ierr)


                !> Close the Domain/Variable group
                call h5gclose_f(gid,ierr)

                !> Deallocate variable buffer storage for the current domain
                deallocate(var)
                idom = idom + 1
            end if
        end do

        !>  Close file
        call h5fclose_f(fid, ierr)

        !>  Close Fortran interface
        call h5close_f(ierr)

    end subroutine




    !> Write HDF5 variable
    !!
    !!  Opens a given hdf5 file. Loads the EquationSet and solution order and calls solution initialization
    !!  procedure for each domain. Searches for the given variable and time instance. If it finds it, load to a
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]      filename    Character string of the file to be read
    !!  @param[in]      cvar        Character string of the variable name to be read
    !!  @param[in]      time        Integer of the time instance for the current variable to be read
    !!  @param[inout]   domains     Array of domains. Already allocated
    !-----------------------------------------------------------------------------------------------------------
    subroutine write_var_hdf(filename,cvar,time,domains)
        character(*),   intent(in)      :: filename
        character(*),   intent(in)      :: cvar
        integer(ik),    intent(in)      :: time
        type(domain_t), intent(inout)   :: domains(:)


        integer(HID_T)   :: fid, gid, sid, did           !> Identifiers
        integer(HSIZE_T) :: dims(2), maxdims(2), adim    !> Dataspace dimensions
        type(H5O_INFO_T) :: info                         !> Object info type

        integer                         :: ndims
        integer, dimension(1)           :: ibuf
        character(100)                  :: cbuf, eqnstring, varstring, var_grp, ctime
        real(rk), allocatable, target   :: var(:,:)
        type(c_ptr)                     :: cp_var
        character(10)                   :: gname
        integer(ik)                     :: nmembers,    type,   ierr,       ndomains,   igrp,   &
                                           npts,        idom,   nterms_1d,  nterms_s,   order,  &
                                           ivar,        ielem
        logical                         :: FileExists, VariablesExists, DataExists, ElementsEqual
        logical                         :: exists




        !>  Check file exists
        !------------------------------------------------------------------
        inquire(file=filename, exist=FileExists)
        if (.not. FileExists) then
            print*, "Error: read_grid_hdf5 - file not found: ", filename
            stop
        end if


        !>  Initialize Fortran interface.
        !------------------------------------------------------------------
        call h5open_f(ierr)
        if (ierr /= 0) stop "Error: write_var_hdf5 - h5open_f"

        !>  Open input file using default properties.
        call h5fopen_f(filename, H5F_ACC_RDWR_F, fid, ierr)
        if (ierr /= 0) stop "Error: write_var_hdf5 - h5fopen_f"

        !>  Get number of domains from attribute 'ndomains' in file root
        !------------------------------------------------------------------
        call h5ltget_attribute_int_f(fid, "/", 'ndomains', ibuf, ierr)
        ndomains = ibuf(1)
        if (ierr /= 0) stop "Error: write_var_hdf5 - h5ltget_attribute_int_f"
        if (ndomains /= size(domains)) stop "Error: write_var_hdf5 - number of domains in file and passed domain list to not match"



        !>  Get number of groups in the file root
        call h5gn_members_f(fid, "/", nmembers, ierr)

        !>  Loop through groups and read domains
        !------------------------------------------------------------------
        idom = 1
        do igrp = 0,nmembers-1
            call h5gget_obj_info_idx_f(fid,"/", igrp, gname, type, ierr)

            if (gname(1:2) == 'D_') then

                !>  Get EquationSet
                cbuf = Domains(idom)%eqnset%name
                call h5ltset_attribute_string_f(fid, trim(gname), 'EquationSet', cbuf, ierr)
                if (ierr /= 0) stop "Error: write_var_hdf5 - h5ltset_attribute_string_f"


                !> Check if 'Variables' group exists
                call h5lexists_f(fid, trim(gname)//"/Variables", exists, ierr)


                !> Open the Domain/Variables group
                if (exists) then
                    ! If 'Variables' group exists then open the existing group
                    call h5gopen_f(fid, trim(gname)//"/Variables", gid, ierr, H5P_DEFAULT_F)
                    if (ierr /= 0) stop "Error: h5gopen_f -- Domain/Grid group did not open properly"
                else
                    ! If 'Variables group does not exist, then create one.
                    call h5gcreate_f(fid, trim(gname)//"/Variables", gid, ierr)
                end if



                !>  Set number of terms in solution expansion
                adim = 1
                ibuf = domains(idom)%mesh%nterms_s
                call h5ltset_attribute_int_f(gid, "/", 'Order', ibuf, adim, ierr)
                if (ierr /= 0) stop "Error: write_var_hdf5 - h5ltset_attribute_int_f"


                !> Compose variable string
                write(ctime, '(I0.3)') time                     !> write time as character string
                varstring = trim(cvar)//'_'//trim(ctime)        !> compose variable name as 'var_time'

                !> Set dimensions of dataspace to write
                ndims = 2
                dims(1) = domains(idom)%mesh%nterms_s
                dims(2) = domains(idom)%mesh%nelem
                maxdims(1) = H5S_UNLIMITED_F
                maxdims(2) = H5S_UNLIMITED_F

                !>  Open the Variable dataset
                ! Check if variable dataset already exists
                call h5lexists_f(gid, trim(varstring), exists, ierr)


                !> Reset dataspace size if necessary
                if (exists) then
                    ! Open the existing dataset
                    call h5dopen_f(gid, trim(varstring), did, ierr, H5P_DEFAULT_F)
                    if (ierr /= 0) stop "Error: write_var_hdf5 -- variable does not exist or was not opened correctly"

                    ! Get the existing dataspace id
                    call h5dget_space_f(did, sid, ierr)

                    ! Resize the existing dataspace
                    call h5sset_extent_simple_f(sid,ndims,dims,maxdims,ierr)
                    if (ierr /= 0) stop "Error: write_var_hdf5 -- dataspace resizing"
                else
                    !> Create a new dataspace
                    call h5screate_simple_f(ndims,dims, sid, ierr)
                    if (ierr /= 0) stop "Error: write_var_hdf5 - h5screate_simple_f"

                    !> Create a new dataset
                    call h5dcreate_f(gid, trim(varstring), H5T_NATIVE_DOUBLE, sid, did, ierr)
                    if (ierr /= 0) stop "Error: write_var_hdf5 - h5dcreate_f"
                end if



                !> Get variable integer index from variable character string
                ivar = domains(idom)%eqnset%prop%get_eqn_index(cvar)

                !> Assemble variable buffer matrix
                allocate(var(dims(1),dims(2)))
                do ielem = 1,Domains(idom)%mesh%nelem
                        !var(:,ielem) = Domains(idom)%sdata%q%lvecs(ielem)%var(ivar)
                        var(:,ielem) = Domains(idom)%sdata%q%lvecs(ielem)%getvar(ivar)
                end do


                !> Write variable buffer
!                call h5dwrite_f(did, H5T_NATIVE_DOUBLE, var, dims, ierr)
                cp_var = c_loc(var(1,1))
                call h5dwrite_f(did, H5T_NATIVE_DOUBLE, cp_var, ierr)
                if (ierr /= 0) stop "Error: write_var_hdf5 - h5dwrite_f"



                !> Close the Variable datasets
                call h5dclose_f(did,ierr)

                !> Close Variable dataspaces
                call h5sclose_f(sid,ierr)

                !> Close the Domain/Variable group
                call h5gclose_f(gid,ierr)

                !> Deallocate variable buffer storage for the current domain
                deallocate(var)
                idom = idom + 1
            end if
        end do

        !>  Close file
        call h5fclose_f(fid, ierr)

        !>  Close Fortran interface
        call h5close_f(ierr)



    end subroutine











end module mod_hdfio
