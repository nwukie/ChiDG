module mod_hdfio
#include <messenger.h>
    use mod_kinds,                  only: rk,ik,rdouble
    use mod_constants,              only: ZERO, NFACES, TWO_DIM, THREE_DIM, NO_PROC
    use mod_bc,                     only: create_bc
    use mod_chidg_mpi,              only: IRANK, NRANK, ChiDG_COMM

    use mod_kinds,                  only: rk,ik,rdouble
    use mod_constants,              only: ZERO, NFACES, TWO_DIM, THREE_DIM, NO_PROC
    use mod_bc,                     only: create_bc
    use mod_chidg_mpi,              only: IRANK, NRANK, ChiDG_COMM
    use mod_hdf_utilities,          only: get_ndomains_hdf, get_domain_names_hdf,                        &
                                          get_domain_equation_set_hdf, set_coordinate_order_hdf,         &
                                          set_solution_order_hdf, get_solution_order_hdf,                &
                                          get_domain_mapping_hdf, get_domain_dimensionality_hdf,         &
                                          set_contains_solution_hdf, set_domain_equation_set_hdf,        &
                                          get_domain_coordinates_hdf, get_domain_coordinate_system_hdf,  &
                                          get_domain_connectivity_hdf, get_domain_nnodes_hdf,            &
                                          check_file_storage_version_hdf, check_file_exists_hdf,         &
                                          get_contains_solution_hdf, get_contains_grid_hdf,              &
                                          get_bc_state_names_hdf, get_bc_state_hdf,                      &
                                          get_nbc_state_groups_hdf, get_bc_state_group_names_hdf,        &
                                          get_bc_patch_group_hdf, get_bc_state_group_family_hdf,         &
                                          get_bc_patch_hdf, open_file_hdf, close_file_hdf,               &
                                          open_domain_hdf, close_domain_hdf, initialize_file_hdf,        &
                                          initialize_file_structure_hdf, open_bc_group_hdf,              &
                                          close_bc_group_hdf, get_domain_nelements_hdf,                  &
                                          get_domain_name_hdf, set_ntimes_hdf, get_ntimes_hdf

    use type_svector,               only: svector_t
    use mod_string,                 only: string_t
    use type_chidg_data,            only: chidg_data_t
    use type_meshdata,              only: meshdata_t
    use type_bc_patch_data,         only: bc_patch_data_t
    use type_bc_group,              only: bc_group_t
    use type_bc_state,              only: bc_state_t
    use type_domain_connectivity,   only: domain_connectivity_t
    use type_partition,             only: partition_t
    use iso_c_binding,              only: c_ptr
    use hdf5
    use h5lt
    use mpi_f08
    implicit none



contains

    !----------------------------------------------------------------------------------------
    !!
    !!  High-Level API for ChiDG HDF File Format
    !!
    !!  Procedures:
    !!  -----------
    !!
    !!  read_grid_hdf
    !!
    !!  read_solution_hdf
    !!      read_domain_field_hdf
    !!
    !!  write_solution_hdf
    !!      write_variable_hdf
    !!      write_domain_field_hdf
    !!
    !!  read_boundaryconditions_hdf
    !!      read_bc_patches_hdf
    !!      read_bc_state_groups_hdf
    !!
    !!  read_connectivity_hdf
    !!
    !!  TODO:
    !!  -----------
    !!      - Relocate solution_order to a particular variable, instead of the domain.
    !!        This is in case a variable is writted from another analysis like a wall
    !!        distance calculation, but then a Navier Stokes solution is started.
    !!      
    !!
    !****************************************************************************************






   
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
    !----------------------------------------------------------------------------------------
    subroutine read_grid_hdf(filename, partition, meshdata)
        character(*),                   intent(in)      :: filename
        type(partition_t),              intent(in)      :: partition
        type(meshdata_t), allocatable,  intent(inout)   :: meshdata(:)

        integer(HID_T)   :: fid, gid, domain_id

        logical                                 :: contains_grid
        character(:),           allocatable     :: user_msg, domain_name
        integer                                 :: ierr, nterms_1d, mapping, &
                                                   nterms_c, spacedim, iconn, &
                                                   nconn, nelements, nnodes


        !
        ! Open file
        !
        fid = open_file_hdf(filename)


        !
        ! Check contains grid
        !
        contains_grid = get_contains_grid_hdf(fid)
        user_msg = "We didn't find a grid to read in the file that was specified. &
                    The file could be a bare ChiDG file or maybe was generated by &
                    an incompatible version of the ChiDG library."
        if (.not. contains_grid) call chidg_signal(FATAL,user_msg)


        !
        !  Allocate a meshdata structure for each domain connectivity in the partition
        !
        nconn = size(partition%connectivities)
        allocate(meshdata(nconn), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        !  Loop through groups and read domains
        !
        do iconn = 1,nconn

            
            !
            ! Open domain. Get name, get hdf identifier
            !
            domain_name = partition%connectivities(iconn)%get_domain_name()
            domain_id = open_domain_hdf(fid,trim(domain_name))


            !
            !  Get number of terms in coordinate expansion
            !
            mapping = get_domain_mapping_hdf(domain_id)
            nterms_1d = (mapping + 1)


            !
            ! Get dimension of the current block: 3D
            !
            spacedim = get_domain_dimensionality_hdf(domain_id)
            nterms_c = nterms_1d * nterms_1d * nterms_1d

            meshdata(iconn)%nterms_c = nterms_c
            meshdata(iconn)%name     = domain_name
            meshdata(iconn)%spacedim = spacedim


            !
            ! Get coordinates
            !
            meshdata(iconn)%points       = get_domain_coordinates_hdf(domain_id)
            meshdata(iconn)%coord_system = get_domain_coordinate_system_hdf(domain_id)


            !
            ! Store connectivity
            !
            nelements = partition%connectivities(iconn)%get_nelements()
            nnodes    = partition%connectivities(iconn)%get_nnodes()
            call meshdata(iconn)%connectivity%init(domain_name,nelements,nnodes)
            meshdata(iconn)%connectivity%data = partition%connectivities(iconn)%data


            !
            ! Read equation set attribute
            !
            meshdata(iconn)%eqnset = get_domain_equation_set_hdf(domain_id)



            !
            ! Close the Domain group
            !
            call close_domain_hdf(domain_id)



        end do  ! iconn



        !  Close file and Fortran interface
        call close_file_hdf(fid)

    end subroutine read_grid_hdf
    !****************************************************************************************
  
  
  
   







    !> Read solution modes from HDF file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      filename    Character string of the file to be read from
    !!  @param[inout]   data        chidg_data_t that will accept the solution modes
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine read_solution_hdf(filename,data)
        character(*),       intent(in)      :: filename
        type(chidg_data_t), intent(inout)   :: data

        integer(HID_T)                  :: fid, domain_id
        integer                         :: ierr

        integer(ik)                     :: idom, ndomains, ieqn, neqns, itime, ntime, eqn_ID
        character(:),       allocatable :: field_name, user_msg, domain_name
        logical                         :: file_exists, contains_solution


        !
        ! Get number of domains contained in the ChiDG data instance
        !
        ndomains = data%ndomains()

        !
        ! Open file
        !
        fid = open_file_hdf(filename)


        !
        ! Check file contains solution
        !
        contains_solution = get_contains_solution_hdf(fid)
        user_msg = "We didn't find a solution to read in the file that was specified. &
                    You could set solutionfile_in = 'none' to initialize a solution &
                    to default values instead."
        if (.not. contains_solution) call chidg_signal(FATAL,user_msg)



        !
        ! Read solution for each time step
        !
        ntime = get_ntimes_hdf(fid)

        do itime = 1, ntime

            !
            ! Read solution for each domain
            !
            do idom = 1,ndomains


                ! Get domain name and number of primary fields
                domain_name = data%info(idom)%name



                ! For each primary field in the domain, get the field name and read from file.
                domain_id = open_domain_hdf(fid,domain_name)

                eqn_ID = data%mesh(idom)%eqn_ID
                !do ieqn = 1,data%eqnset(idom)%prop%nprimary_fields()
                do ieqn = 1,data%eqnset(eqn_ID)%prop%nprimary_fields()
                    !field_name = trim(data%eqnset(idom)%prop%get_primary_field_name(ieqn))
                    field_name = trim(data%eqnset(eqn_ID)%prop%get_primary_field_name(ieqn))
                    call read_domain_field_hdf(data,domain_id,field_name,itime,'Primary')
                end do ! ieqn

    !            do ieqn = 1,data%eqnset(eqn_ID)%prop%nauxiliary_fields()
    !                field_name = trim(data%eqnset(eqn_ID)%prop%get_primary_field_name(ieqn))
    !                call read_field_domain_hdf(data,domain_id,field_name,itime,'Auxiliary')
    !            end do ! ieqn

                call close_domain_hdf(domain_id)

            end do ! idom

        end do ! itime

        ! Close file
        call close_file_hdf(fid)


    end subroutine read_solution_hdf
    !****************************************************************************************












    !> Write solution modes to HDF file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      filename    Character string of the file to be written to
    !!  @param[inout]   data        chidg_data_t containing solution to be written
    !!
    !!
    !!  @autor Matteo Ugolotti
    !!  @date   2/22/2017
    !!
    !----------------------------------------------------------------------------------------
    subroutine write_solution_hdf(data,file_name,field)
        type(chidg_data_t), intent(in)              :: data
        character(*),       intent(in)              :: file_name
        character(*),       intent(in), optional    :: field


        character(:),   allocatable     :: field_name, domain_name
        integer(HID_T)                  :: fid, domain_id
        integer(HSIZE_T)                :: adim
        integer(ik)                     :: idom, ieqn, neqns, iwrite, spacedim, time, field_index, iproc, eqn_ID
        integer                         :: ierr, order_s
        logical                         :: file_exists
        integer(ik)                     :: itime

        !
        ! Check for file existence
        !
        file_exists = check_file_exists_hdf(file_name)


        !
        ! Create new file if necessary
        !
        if (.not. file_exists) then

                ! Create a new file
                if (IRANK == GLOBAL_MASTER) then
                    call initialize_file_hdf(file_name)
                end if
                call MPI_Barrier(ChiDG_COMM,ierr)

                ! Initialize the file structure.
                do iproc = 0,NRANK-1
                    if (iproc == IRANK) then
                        fid = open_file_hdf(file_name)
                        call initialize_file_structure_hdf(fid,data)
                        call close_file_hdf(fid)
                    end if
                    call MPI_Barrier(ChiDG_COMM,ierr)
                end do

            else


        end if
        call MPI_Barrier(ChiDG_COMM,ierr)



        !
        ! Each process, write its own portion of the solution
        !
        do iwrite = 0,NRANK-1
            if ( iwrite == IRANK ) then


                fid = open_file_hdf(file_name)
                !
                ! Set the attribute times to the hdf file
                !
                call set_ntimes_hdf(fid,data%ntime())

                !
                ! Write solution for each domain
                !
                do idom = 1,data%ndomains()

                    domain_name = data%info(idom)%name
                    domain_id   = open_domain_hdf(fid,trim(domain_name))
                    eqn_ID      = data%mesh(idom)%eqn_ID
                    

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



                    !
                    ! Set some data about the block: solution order + equation set
                    !
                    call set_solution_order_hdf(domain_id,order_s)
                    
                    !
                    ! Loop through time levels
                    !
                    do itime = 1,data%ntime()


                        !
                        ! If specified, only write specified field.
                        !
                        if (present(field)) then

                            !field_index = data%eqnset(idom)%prop%get_primary_field_index(trim(field))
                            field_index = data%eqnset(eqn_ID)%prop%get_primary_field_index(trim(field))

                            if (field_index /= 0) then
                                call write_domain_field_hdf(domain_id,data,field,itime)
                            end if


                        !
                        ! Else, write each field in the file.
                        !
                        else

                            !
                            ! For each field: get the name, write to file
                            ! 
                            !neqns = data%eqnset(idom)%prop%nprimary_fields()
                            neqns = data%eqnset(eqn_ID)%prop%nprimary_fields()
                            do ieqn = 1,neqns
                                !field_name = trim(data%eqnset(idom)%prop%get_primary_field_name(ieqn))
                                field_name = trim(data%eqnset(eqn_ID)%prop%get_primary_field_name(ieqn))
                                call write_domain_field_hdf(domain_id,data,field_name,itime)
                            end do ! ieqn

                        end if

                    end do ! itime

                    call close_domain_hdf(domain_id)


                end do ! idom


                call set_contains_solution_hdf(fid,"True")
                call close_file_hdf(fid)

            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
        end do

    end subroutine write_solution_hdf
    !*****************************************************************************************

   
   
   
   
   
  
  
   
    !>  Read HDF5 variable
    !!
    !!  Opens a given ChiDG-formatted HDF file. Loads the equation set and solution order
    !!  and calls solution initialization procedure for each domain. Searches for the given
    !!  variable and time instance. If it finds it, load to a
    !!
    !!  Note: Convention is that all floating-point data is double precision format.
    !!        Conversion to working-precision should happen after reading the data from 
    !!        the HDF file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[inout]   data        ChiDG data containing domains. Already allocated.
    !!  @param[in]      domain_id   HDF5 Domain identifier.
    !!  @param[in]      field_name  Character string of the field to be read.
    !!  @param[in]      itime       Integer of the time instance for the current variable 
    !!                              to be read.
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine read_domain_field_hdf(data,domain_id,field_name,itime,field_type)
        type(chidg_data_t),         intent(inout)   :: data
        integer(HID_T),             intent(in)      :: domain_id
        character(*),               intent(in)      :: field_name
        integer(ik),                intent(in)      :: itime
        character(*),               intent(in)      :: field_type


        integer(HID_T)          :: gid, sid, vid, memspace
        integer(HSIZE_T)        :: start(3), count(3), dimsm(3)
        integer, dimension(1)   :: ibuf


        character(:),   allocatable         :: user_msg, domain_name
        character(100)                      :: cbuf, var_gqp

        real(rdouble),  allocatable, target :: var(:,:,:)
        real(rdouble),  allocatable         :: bufferterms(:)
        type(c_ptr)                         :: cp_var

        integer(ik)                         :: spacedim, ielem_g, aux_vector_index, eqn_ID
        integer                             :: type, ierr, nterms_1d, nterms_s, order,  &
                                               ivar, ielem, nterms_ielem, idom, ndims
        logical                             :: ElementsEqual, variables_exists


        !
        ! Check valid field_type input
        !
        if ( (trim(field_type) /= 'Primary') .and. &
             (trim(field_type) /= 'Auxiliary') ) then
             user_msg = "read_field_domain_hdf: An invalid field type was passed to the routine. &
                         valid field types are 'Primary' and 'Auxiliary'."
             call chidg_signal_one(FATAL,user_msg,trim(field_type))
        end if



        !
        ! Check if 'Variables' group exists
        !
        call h5lexists_f(domain_id, "Variables", variables_exists, ierr)
        if (.not. variables_exists) call chidg_signal(FATAL,"read_field_domain_hdf: Variables group does not exist")


        !
        ! Open the Domain/Variables group
        !
        call h5gopen_f(domain_id, "Variables", gid, ierr, H5P_DEFAULT_F)
        if (ierr /= 0) call chidg_signal(FATAL,"read_field_domain_hdf: h5gopen_f -- Variables group did not open properly")


        !
        ! Get number of terms in solution expansion
        !
        order = get_solution_order_hdf(domain_id)
        nterms_1d = (order + 1) ! To be consistent with the definition of (Order = 'Order of the polynomial')


        domain_name = get_domain_name_hdf(domain_id)
        idom     = data%get_domain_index(domain_name)
        spacedim = data%mesh(idom)%spacedim

        if ( spacedim == THREE_DIM ) then
            nterms_s = nterms_1d*nterms_1d*nterms_1d
        else if ( spacedim == TWO_DIM ) then
            nterms_s = nterms_1d*nterms_1d
        end if


        
        !
        ! Open the Variable dataset
        !
        call h5dopen_f(gid, trim(field_name), vid, ierr, H5P_DEFAULT_F)
        if (ierr /= 0) call chidg_signal(FATAL,"read_field_domain_hdf: variable does not exist or was not opened correctly")


        !
        ! Get the dataspace id and dimensions
        !
        call h5dget_space_f(vid, sid, ierr)





        if (field_type == 'Auxiliary') then 
            aux_vector_index = data%sdata%get_auxiliary_field_index(trim(field_name))
        end if




        !
        ! Allocate storage for incoming element modes, create a memory dataspace (memspace)
        !
        allocate(var(nterms_s,1,1))
        cp_var = c_loc(var(1,1,1))

        ndims    = 3
        dimsm(1) = size(var,1)
        dimsm(2) = size(var,2)
        dimsm(3) = size(var,3)
        call h5screate_simple_f(ndims,dimsm,memspace,ierr)




        !
        !  Loop through elements and set 'variable' values
        !
        eqn_ID = data%mesh(idom)%eqn_ID
        do ielem = 1,data%mesh(idom)%nelem


            !
            ! Get number of terms initialized for the current element
            !
            if (field_type == 'Primary') then
                !nterms_ielem = data%sdata%q%dom(idom)%vecs(ielem)%nterms()
                nterms_ielem = data%sdata%q_in%dom(idom)%vecs(ielem)%nterms()
            else if (field_type == 'Auxiliary') then
                nterms_ielem = data%sdata%auxiliary_field(aux_vector_index)%dom(idom)%vecs(ielem)%nterms() 
            end if


            !
            ! get domain-global element index
            !
            ielem_g = data%mesh(idom)%elems(ielem)%ielement_g
            start = [1-1,ielem_g-1,itime-1] ! 0-based 
            count = [nterms_s, 1, 1]
            
!            !
!            ! Select the the correct time levels in the HDF5 file based on the time_integrator
!            ! used in the previous solution
!            !
!            if ( switch == 'steady' ) then
!                start = [1-1,ielem_g-1,hdf_time]   ! 0-based
!            else if ( switch == 'HB' ) then
!                start = [1-1,ielem_g-1,itime-1]   ! 0-based
!            else !( switch == 'time_marching' )
!                start = [1-1,ielem_g-1,(hdf_time + itime)-1]   ! 0-based
!            end if
!
!            !
!            ! Chunk of data in domain/Variable dataspace to be read
!            !
!            count = [nterms_s, 1, 1]


            !
            ! Allocate bufferterm storage that will be used to set variable data
            !
            if (allocated(bufferterms)) deallocate(bufferterms)
            allocate(bufferterms(nterms_ielem), stat=ierr)
            if (ierr /= 0) call AllocationError




            !
            ! Select subset of dataspace - sid, read selected modes into cp_var 
            !
            call h5sselect_hyperslab_f(sid, H5S_SELECT_SET_F, start, count, ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"read_domain_field_hdf: h5sselect_hyperslab_f.")
            call h5dread_f(vid, H5T_NATIVE_DOUBLE, cp_var, ierr, memspace, sid)
            if (ierr /= 0) call chidg_signal(FATAL,"read_domain_field_hdf: h5d_read_f_f.")



            !
            ! Check for reading lower, higher, or same-order solution
            !
            bufferterms = ZERO
            if ( nterms_s < nterms_ielem ) then
                ! Reading a lower-order solution
                bufferterms(1:nterms_s) = var(1:nterms_s, 1, 1)
            else if ( nterms_s > nterms_ielem ) then
                ! Reading a higher-order solution
                bufferterms(1:nterms_ielem) = var(1:nterms_ielem, 1, 1)
            else
                ! Reading a solution of same order
                bufferterms(1:nterms_ielem) = var(1:nterms_ielem, 1, 1)
            end if


            ! Store modes in ChiDG Vector
            if (field_type == 'Primary') then
                !ivar = data%eqnset(idom)%prop%get_primary_field_index(trim(field_name))
                !call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar,itime,real(bufferterms,rk))
                ivar = data%eqnset(eqn_ID)%prop%get_primary_field_index(trim(field_name))
                call data%sdata%q_in%dom(idom)%vecs(ielem)%setvar(ivar,itime,real(bufferterms,rk))
            else if (field_type == 'Auxiliary') then
                ! Implicitly assuming that an auxiliary field chidgVector contains only one field.
                ivar = 1
                call data%sdata%auxiliary_field(aux_vector_index)%dom(idom)%vecs(ielem)%setvar(ivar,itime,real(bufferterms,rk))
            end if



        end do







        !
        ! Close variable dataset, domain/variable group.
        !
        call h5sclose_f(memspace,ierr)  ! Close memory space
        call h5dclose_f(vid,ierr)       ! Close variable dataset
        call h5sclose_f(sid,ierr)       ! Close Variable dataspaces
        call h5gclose_f(gid,ierr)       ! Close Domain/Variable group


    end subroutine read_domain_field_hdf
    !*************************************************************************************
   
   
   
   
   
   
  






    !>  Write HDF5 field to a Domain.
    !!
    !!  Loads the equation set and solution order and calls solution
    !!  initialization procedure for each domain. Searches for the given variable and time
    !!  instance.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  data        chidg_data_t instance containing grid and solution.
    !!  @param[in]  domain_id   HDF5 file identifier.
    !!  @param[in]  field_name  Character string of the variable name to be read.
    !!  @param[in]  itime       Integer of the time instance for the current variable 
    !!                          to be read.
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   02/20/2017
    !!
    !!  Expanded the dataspace to include all the time levels
    !!
    !----------------------------------------------------------------------------------------
    subroutine write_domain_field_hdf(domain_id,data,field_name,itime,attribute_name,attribute_value)
        integer(HID_T),     intent(in)              :: domain_id
        type(chidg_data_t), intent(in)              :: data
        character(*),       intent(in)              :: field_name
        integer(ik),        intent(in)              :: itime
        character(*),       intent(in), optional    :: attribute_name
        real(rk),           intent(in), optional    :: attribute_value


        type(H5O_INFO_T) :: info
        integer(HID_T)   :: gid, sid, did, crp_list, memspace, filespace
        integer(HSIZE_T) :: edims(2), maxdims(3), dims(3), dimsm(3), dimsc(3), &
                            start(3), count(3)

        integer                             :: ndims, ibuf(1)
        character(100)                      :: cbuf, var_grp, ctime
        character(:),   allocatable         :: domain_name
        real(rdouble),  allocatable, target :: var(:,:,:)
        type(c_ptr)                         :: cp_var

        logical     :: DataExists, ElementsEqual, exists
        integer(ik) :: type, ierr, ndomains,    &
                       order, ivar, ielem, idom, nelem_g, &
                       ielement_g, nterms_s, ntime, eqn_ID


        !
        ! Open the Domain/Variables group
        !
        call h5lexists_f(domain_id, "Variables", exists, ierr)
        if (exists) then
            call h5gopen_f(domain_id, "Variables", gid, ierr, H5P_DEFAULT_F)
        else
            call h5gcreate_f(domain_id, "Variables", gid, ierr)
        end if
        if (ierr /= 0) call chidg_signal(FATAL,"write_field_domain_hdf: Domain/Variables group did not open properly.")



        !
        ! Set dimensions of dataspace to write
        !
        domain_name = get_domain_name_hdf(domain_id)
        idom        = data%get_domain_index(domain_name)
        nelem_g     = data%mesh(idom)%get_nelements_global()
        nterms_s    = data%mesh(idom)%nterms_s
        ndims       = 3
        ntime       = data%ntime()

        !
        ! Set the dimensions of the dataspace to write in
        !
        dims(1:3)    = [nterms_s, nelem_g, ntime]
        maxdims(1:3) = H5S_UNLIMITED_F


        !
        ! Modify dataset creation properties, i.e. enable chunking in order to append
        ! dataspace, if needed.
        !
        dimsc = [1, nelem_g, 1]  ! Chunk size

        call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, ierr)
        if (ierr /= 0) call chidg_signal(FATAL, "write_field_domain_hdf: h5pcreate_f error enabling chunking.")

        call h5pset_chunk_f(crp_list, ndims, dimsc, ierr)
        if (ierr /= 0) call chidg_signal(FATAL, "write_field_domain_hdf: h5pset_chunk_f error setting chunk properties.")



        !
        ! Reset dataspace size if necessary
        !
        call h5lexists_f(gid, trim(field_name), exists, ierr)
        if (exists) then
            ! Open the existing dataset
            call h5dopen_f(gid, trim(field_name), did, ierr, H5P_DEFAULT_F)
            if (ierr /= 0) call chidg_signal(FATAL,"write_field_domain_hdf: variable does not exist or was not opened correctly.")


            ! Extend dataset if necessary
            call h5dset_extent_f(did, dims, ierr)
            if (ierr /= 0) call chidg_signal(FATAL, "write_field_domain_hdf: h5dset_extent_f.")


            ! Update existing dataspace ID since it may have been expanded
            call h5dget_space_f(did, sid, ierr)
            if (ierr /= 0) call chidg_signal(FATAL, "write_field_domain_hdf: h5dget_space_f.")

        else
            ! Create a new dataspace
            call h5screate_simple_f(ndims,dims,sid,ierr,maxdims)
            if (ierr /= 0) call chidg_signal(FATAL,"write_field_domain_hdf: h5screate_simple_f.")


            ! Create a new dataset
            call h5dcreate_f(gid, trim(field_name), H5T_NATIVE_DOUBLE, sid, did, ierr, crp_list)
            if (ierr /= 0) call chidg_signal(FATAL,"write_field_domain_hdf: h5dcreate_f.")
        end if



        !
        ! Get variable integer index from variable character string
        !
        eqn_ID = data%mesh(idom)%eqn_ID
        ivar = data%eqnset(eqn_ID)%prop%get_primary_field_index(field_name)



        !
        ! Assemble variable buffer matrix that gets written to file
        !
        allocate(var(nterms_s,1,1))


        do ielem = 1,data%mesh(idom)%nelem

            !
            ! get domain-global element index
            !
            ielement_g = data%mesh(idom)%elems(ielem)%ielement_g
            start = [1-1,ielement_g-1,itime-1]   ! 0-based
            count = [nterms_s, 1, 1]

            !
            ! Select subset of dataspace - sid
            !
            call h5sselect_hyperslab_f(sid, H5S_SELECT_SET_F, start, count, ierr)

            !
            ! Create a memory dataspace
            !
            dimsm(1) = size(var,1)
            dimsm(2) = size(var,2)
            dimsm(3) = size(var,3)
            call h5screate_simple_f(ndims,dimsm,memspace,ierr)


            !
            ! Write modes
            !
            var(:,1,1) = real(data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar,itime),rdouble)
            cp_var = c_loc(var(1,1,1))
            call h5dwrite_f(did, H5T_NATIVE_DOUBLE, cp_var, ierr, memspace, sid)


            call h5sclose_f(memspace,ierr)


        end do


        !
        ! Write order of the domain variable
        !
        !call h5set_attribute



        call h5pclose_f(crp_list, ierr) ! Close dataset creation property
        call h5dclose_f(did,ierr)       ! Close Variable datasets
        call h5sclose_f(sid,ierr)       ! Close Variable dataspaces
        call h5gclose_f(gid,ierr)       ! Close Domain/Variable group


    end subroutine write_domain_field_hdf
    !****************************************************************************************

















    !>  Read boundary conditions from HDF5 file in ChiDG format and return data in bcdata_t
    !!  container. The calling procedure can then use the returned bcdata_t to initialize
    !!  boundary conditions.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/9/2016
    !!
    !!  @param[in]      filename    String of the HDF5 file to be read.
    !!  @param[inout]   bcdata(:)   Array of bcdata_t instances, one for each domain. 
    !!                              These will be returned with data about the boundary
    !!                              conditions that can be used for initialization.
    !!  @param[in]      partition   Partition information to only read boundary conditions 
    !!                              for the domains in the partition
    !!
    !----------------------------------------------------------------------------------------
    subroutine read_boundaryconditions_hdf(filename, bc_patch_data, bc_groups, partition)
        character(*),           intent(in)                  :: filename
        type(bc_patch_data_t),  intent(inout), allocatable  :: bc_patch_data(:)
        type(bc_group_t),       intent(inout), allocatable  :: bc_groups(:)
        type(partition_t),      intent(in)                  :: partition

        character(len=10)       :: faces(NFACES)
        integer(HID_T)          :: fid
        integer                 :: ierr, nconn


        faces = ["  XI_MIN","  XI_MAX"," ETA_MIN"," ETA_MAX","ZETA_MIN","ZETA_MAX"]


        fid = open_file_hdf(filename)


        !
        !  Allocate for number of domains in the partition
        !
        nconn = size(partition%connectivities)
        allocate(bc_patch_data(nconn), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Read boundary condition patches
        !
        call read_bc_patches_hdf(fid,bc_patch_data,partition)


        !
        ! Read boundary condition state groups
        !
        call read_bc_state_groups_hdf(fid,bc_groups,partition)



        call close_file_hdf(fid)

    end subroutine read_boundaryconditions_hdf
    !****************************************************************************************










    !>  Read the boundary condition patch connectivity data from file and store in
    !!  bc_patch_data.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/31/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine read_bc_patches_hdf(fid, bc_patch_data, partition)
        integer(HID_T),         intent(in)      :: fid
        type(bc_patch_data_t),  intent(inout)   :: bc_patch_data(:)
        type(partition_t),      intent(in)      :: partition

        integer(ik)                 :: iconn, nconn, iface, ierr
        integer(ik),    allocatable :: bc_patch(:,:)
        character(:),   allocatable :: bc_state_group
        integer                     :: ibc_face, nbcfaces
        character(1024)             :: domain
        character(len=10)           :: patches(NFACES)

        integer(HID_T)              :: patch_id


        patches = ["  XI_MIN","  XI_MAX"," ETA_MIN"," ETA_MAX","ZETA_MIN","ZETA_MAX"]


        !
        !  Loop through connectivities and read boundary conditions
        !
        nconn = size(partition%connectivities)
        do iconn = 1,nconn


            !
            ! Get name of current domain
            !
            domain = partition%connectivities(iconn)%get_domain_name()
            bc_patch_data(iconn)%domain_name = domain


            !
            ! Allocation bcs for current domain
            !
            allocate(bc_patch_data(iconn)%bc_connectivity(NFACES), stat=ierr)
            if (ierr /= 0) call AllocationError


            !
            ! Loop faces and get boundary condition for each
            !
            ! TODO: should probably turn this into a loop over bcs instead of faces.
            do iface = 1,NFACES


                ! Open face boundary condition group
                call h5gopen_f(fid, "D_"//trim(domain)//"/BoundaryConditions/"//trim(adjustl(patches(iface))), patch_id, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"read_bc_patches_hdf: error opening boundary face group")
    

                ! Get bc patch connectivity for current face
                bc_patch = get_bc_patch_hdf(patch_id)
                nbcfaces = size(bc_patch,1)


                ! Store boundary condition connectivity
                call bc_patch_data(iconn)%bc_connectivity(iface)%init(nbcfaces)
                do ibc_face = 1,nbcfaces
                    bc_patch_data(iconn)%bc_connectivity(iface)%data(ibc_face)%data = bc_patch(ibc_face,:)
                end do

                
                ! Read Boundary State Group
                bc_state_group = get_bc_patch_group_hdf(patch_id)
                call bc_patch_data(iconn)%bc_group_name%push_back(string_t(bc_state_group))


                ! Close face boundary condition group
                call h5gclose_f(patch_id, ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"read_bc_patches_hdf: h5gclose")


            end do ! iface


        end do  ! iconn



    end subroutine read_bc_patches_hdf
    !****************************************************************************************














    !>  Read boundary condition state groups from file and return.
    !!
    !!  Reads and returns all boundary condition state groups, regardless of if they
    !!  are used on the current partition or not.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/31/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine read_bc_state_groups_hdf(fid, bc_groups, partition)
        integer(HID_T),         intent(in)                  :: fid
        type(bc_group_t),       intent(inout), allocatable  :: bc_groups(:)
        type(partition_t),      intent(in)                  :: partition

        type(svector_t)                     :: bc_group_names, bc_state_names
        type(string_t)                      :: group_name, state_name
        class(bc_state_t),  allocatable     :: bc

        integer(HID_T)  :: group_id
        integer(ik)     :: igroup, ngroups, istate, ierr


        ngroups        = get_nbc_state_groups_hdf(fid)
        bc_group_names = get_bc_state_group_names_hdf(fid)


        if (allocated(bc_groups)) deallocate(bc_groups)
        allocate(bc_groups(ngroups), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Read each group of bc_state's
        !
        do igroup = 1,ngroups

            ! Open face boundary condition group
            group_name = bc_group_names%at(igroup)
            group_id = open_bc_group_hdf(fid,group_name%get())

            !
            ! Get bc_group Family attribute.
            !
            bc_groups(igroup)%family = get_bc_state_group_family_hdf(group_id)

            !
            ! Loop through and read states + their properties
            !
            bc_state_names = get_bc_state_names_hdf(group_id)
            do istate = 1,bc_state_names%size()

                ! Get bc_state name, return bc_state from file and source-allocate
                state_name = bc_state_names%at(istate)
                if (allocated(bc)) deallocate(bc)
                allocate(bc, source = get_bc_state_hdf(group_id,state_name%get()))

                ! Save to bc_group_data_t
                bc_groups(igroup)%name = group_name%get()
                call bc_groups(igroup)%bc_states%push_back(bc)

            end do !istate


            ! Close face boundary condition group
            call close_bc_group_hdf(group_id)

        end do !igroup



    end subroutine read_bc_state_groups_hdf
    !****************************************************************************************












    !>  This reads an HDF ChiDG grid file and returns an array of connectivities, one 
    !!  for each domain.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/9/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine read_connectivity_hdf(filename, connectivities)
        character(*),                               intent(in)      :: filename
        type(domain_connectivity_t), allocatable,   intent(inout)   :: connectivities(:)

        integer(HID_T)   :: fid, domain_id

        integer(ik),            allocatable :: connectivity(:,:)
        character(len=1024),    allocatable :: domain_names(:)
        character(:),           allocatable :: user_msg, domain_name
        integer                             :: type, ierr, ndomains,    &
                                               idom, idomain, nelements, ielem, nnodes, mapping
        logical                             :: contains_grid



        fid = open_file_hdf(filename)


        ! Check contains grid
        contains_grid = get_contains_grid_hdf(fid)
        user_msg = "We didn't find a grid to read in the file that was specified. &
                    The file could be a bare ChiDG file or maybe was generated by &
                    an incompatible version of the ChiDG library."
        if (.not. contains_grid) call chidg_signal(FATAL,user_msg)



        !
        !  Allocate number of domains
        !
        ndomains = get_ndomains_hdf(fid)
        user_msg = "read_connectivity_hdf: No domains were found in the file."
        if (ndomains == 0) call chidg_signal(FATAL,user_msg)
        allocate(connectivities(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError




        !
        !  Loop through groups and read domain connectivities
        !
        domain_names = get_domain_names_hdf(fid)
        do idom = 1,size(domain_names)

                domain_name = domain_names(idom)
                domain_id = open_domain_hdf(fid,trim(domain_name))


                !
                ! Get connectivity and total number of nodes in the domain
                !
                connectivity = get_domain_connectivity_hdf(domain_id)
                nnodes       = get_domain_nnodes_hdf(domain_id)


                !
                ! Initialize domain connectivity structure
                !
                nelements = size(connectivity,1)
                call connectivities(idom)%init(domain_name,nelements, nnodes)


                do ielem = 1,nelements
                    mapping = connectivity(ielem,3)
                    call connectivities(idom)%data(ielem)%init(mapping)
                    connectivities(idom)%data(ielem)%data = connectivity(ielem,:)
                    call connectivities(idom)%data(ielem)%set_element_partition(NO_PROC)
                end do


                ! Close domain
                call close_domain_hdf(domain_id)

        end do  ! idom


        ! Close file
        call close_file_hdf(fid)

    end subroutine read_connectivity_hdf
    !****************************************************************************************








    !>  Read the weights of each domain.
    !!
    !!  The weights here are defined as relative weight of compute intensity.
    !!  So, for example, a domain solving the Euler equations might be weighted 1, 
    !!  whereas a domain solving the Navier Stokes equations might be weighted 5.
    !!  Additionally, the weight might be based on solution order. So, a P1 domain might
    !!  be weighted 1 and a P2 domain might be weighted 8. 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/25/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine read_weights_hdf(chidg_file,weights)
        character(*),               intent(in)      :: chidg_file
        real(rk),   allocatable,    intent(inout)   :: weights(:)






    end subroutine read_weights_hdf
    !*****************************************************************************************





end module mod_hdfio
