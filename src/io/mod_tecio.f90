!--------------------------------------------------------------------------------------
!!
!!  API for writing tecplot files using the TecIO library.
!!
!!  Right now, these capabilities exist strictly to support post-processing
!!  and visualization of solution data. Modal data is interpolated to 
!!  discrete points and elements/faces are sub-divided to resolve polynomial
!!  variation within these entities.
!!
!!  Procedure hierarchy:
!!  ----------------------
!!  write_tecio_file
!!      : init_tecio_file
!!      : write_tecio_domains
!!          : init_tecio_volume_zones
!!          : init_tecio_volume_partition
!!      : write_tecio_surfaces
!!          : init_tecio_surface_zones
!!          : init_tecio_surface_partition
!!      : close_tecio_file
!!
!!
!!  @author Nathan A. Wukie
!!
!--------------------------------------------------------------------------------------
module mod_tecio
#include <messenger.h>
    use mod_kinds,              only: rk,ik,rdouble,TEC
    use mod_constants,          only: ONE, HALF, TWO, OUTPUT_RES, NO_ID, CYLINDRICAL, NO_PROC, NO_DIFF
    use mod_chidg_mpi,          only: GLOBAL_MASTER, IRANK
    use mpi_f08,                only: MPI_AllReduce, MPI_INTEGER4, MPI_LOGICAL, MPI_MAX, MPI_CHARACTER, &
                                      mpi_comm, MPI_LOR, MPI_Send, MPI_Recv, MPI_STATUS_IGNORE

    use type_chidg_data,        only: chidg_data_t
    use type_domain,            only: domain_t
    use type_bc_patch_group,    only: bc_patch_group_t
    use type_bc_state_group,    only: bc_state_group_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_chidg_cache,       only: chidg_cache_t
    use type_cache_handler,     only: cache_handler_t
    use type_element_info,      only: element_info_t, element_info
    use type_timer,             only: timer_t

    use mod_string,     only: string_t
    use type_svector,   only: svector_t
    use type_ivector,   only: ivector_t
    use DNAD_D
    use iso_c_binding
    implicit none

#include "tecio.f90"


    ! This gets incremented every time a new zone is initialized. 
    ! This get zeroed every time write_tecio is called.
    integer(TEC) :: strandID = 0

contains



    !>  Write ChiDG data to tecplot file for visualization.
    !!
    !!  NOTE: This processes solution data from data%sdata%q_out. As such, q_out shall
    !!        have been initialized appropriately.
    !!
    !!  @param[inout]   data            chidg_data instance containing initialized mesh/solution.
    !!  @param[in]      filename        String for the name of the file being written to.
    !!  @param[in]      write_domains   logical controlling if volume domain data is written.
    !!  @param[in]      write_surfaces  logical controlling if surface bc data is written.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/7/2016
    !!
    !-----------------------------------------------------------------------------------
    subroutine write_tecio_file(data,filename,write_domains,write_surfaces)
        type(chidg_data_t),     intent(inout)   :: data
        character(*),           intent(in)      :: filename
        logical,                intent(in)      :: write_domains
        logical,                intent(in)      :: write_surfaces

        integer(ik)     :: ierr, ieq, eqn_ID, eqn_ID_global, ndomains_g
        type(c_ptr)     :: handle
        character(1000) :: varstring
        type(timer_t)   :: timer
        logical         :: surface_switch, domain_switch


        call timer%start()

        ! Default switch overrides to true
        domain_switch  = .true.
        surface_switch = .true.

!        ! Turn of write_surfaces if writing in parallel. Current TecIO library does not 
!        ! allow 2D FE zones to be partitioned.
!        if (write_surfaces .and. NRANK > 1) then
!            surface_switch = .false.
!            if (IRANK == GLOBAL_MASTER) then
!                call chidg_signal(MSG,'write_tecio_file: NOTE: surface zones can only be written in serial due to current restrictions in TecIO library. Writing surface zones has been disabled.')
!            end if
!        end if


        ! Assemble variables string.
        !
        !   Default: Grid coordinates
        ieq = 1
        varstring = "X,Y,Z"

        ! Synchronize on equation set to get io_fields from
        if (data%mesh%ndomains() > 0) eqn_ID = data%mesh%domain(1)%elems(1)%eqn_ID
        call MPI_AllReduce(eqn_ID,eqn_ID_global,1,MPI_INTEGER4,MPI_MAX,ChiDG_COMM,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'write_tecio_file: error reducing eqn_ID to synchronize across ranks.')

        do while (ieq <= data%eqnset(eqn_ID_global)%prop%nio_fields())
            varstring = trim(varstring)//","//trim(data%eqnset(eqn_ID_global)%prop%get_io_field_name(ieq))
            ieq = ieq + 1
        end do

        
        ! Open and initialize TecIO file
        handle = init_tecio_file('solnfile',trim(varstring),filename,0)


        ! Write volume data from 'mesh%domains'
        if (write_domains .and. domain_switch) call write_tecio_domains(handle,data,eqn_ID_global,ndomains_g)


        ! Write surface data from 'mesh%domains'
        if (write_surfaces .and. surface_switch) call write_tecio_surfaces(handle,data,eqn_ID_global,ndomains_g)


        ! Close the current TecIO file context
        call close_tecio_file(handle)


        call timer%stop()
        call timer%report('Time:')

    end subroutine write_tecio_file
    !***********************************************************************************










    !>  Write tecio domains.
    !!
    !!  NOTE: This routine assumes that init_tecio_file has already been called.
    !!
    !!  For each domain instance in chidg_data, write a volume zone containing the
    !!  solution for that domain.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/7/2016
    !!
    !-----------------------------------------------------------------------------------
    subroutine write_tecio_domains(handle,data,eqn_ID,ndomains_g)
        type(c_ptr),            intent(in)      :: handle
        type(chidg_data_t),     intent(inout)   :: data
        integer(ik),            intent(in)      :: eqn_ID
        integer(ik),            intent(inout)   :: ndomains_g


        integer(ik)        :: nelem_xi, nelem_eta, nelem_zeta,  &
                              ielem_xi, ielem_eta, ielem_zeta,  &
                              ielem, ielem_global, ielem_offset,&
                              npts_element, nsub_per_element,   &
                              nsub_elements, ierr,              &
                              nelem, istart, ielem_start,       &
                              idom, itime, icoord, inode, ifield

        integer(c_int32_t)                          :: ipartition, zone_index, tecstat
        integer(TEC),   allocatable                 :: connectivity(:,:)
        real(rdouble),  allocatable, dimension(:)   :: val, r, theta, z

        type(AD_D),     allocatable :: var(:)
        character(:),   allocatable :: var_string

        type(chidg_worker_t)        :: worker
        type(chidg_cache_t)         :: cache
        type(cache_handler_t)       :: cache_handler
        type(element_info_t)        :: elem_info



        call write_line("   TecIO: Writing domains...",io_proc=GLOBAL_MASTER)

        ! Initialize Chidg Worker references
        call worker%init(data%mesh, data%eqnset(:)%prop, data%sdata, data%time_manager, cache)

        ! Loop time/domains/elements/fields
        do itime = 1,data%sdata%q_out%get_ntime()
            worker%itime = itime

            ! Initialize another set of zones for the next time instance
            ndomains_g = init_tecio_volume_zones(handle,data,eqn_ID,itime)

            ! Write proc-local domain partitions
            do idom = 1,data%mesh%ndomains()


                ! Initialize new domain partition in the TecIO file
                zone_index = int(data%mesh%domain(idom)%idomain_g + (itime-1)*ndomains_g,c_int32_t)
                ipartition = init_tecio_volume_partition(handle,data%mesh%domain(idom),zone_index)


                ! For each coordinate, compute it's value pointwise and save
                ! For each actual element, create a sub-sampling of elements to resolve solution variation
                nelem = data%mesh%domain(idom)%nelem
                do ielem = 1,nelem
                    do icoord = 1,3

                        ! Get coordinate value at point
                        if ( data%mesh%domain(idom)%elems(ielem)%coordinate_system == CYLINDRICAL ) then
                            r     = real(data%mesh%domain(idom)%elems(ielem)%interp_coords_def(:,1),rdouble)
                            theta = real(data%mesh%domain(idom)%elems(ielem)%interp_coords_def(:,2),rdouble)
                            z     = real(data%mesh%domain(idom)%elems(ielem)%interp_coords_def(:,3),rdouble)
                            if (icoord == 1) val = r*cos(theta)
                            if (icoord == 2) val = r*sin(theta)
                            if (icoord == 3) val = z

                        else
                            if (allocated(val)) deallocate(val)
                            allocate(val(size(data%mesh%domain(idom)%elems(ielem)%interp_coords_def,1)))
                            val = real(data%mesh%domain(idom)%elems(ielem)%interp_coords_def(:,icoord),rdouble)
                        end if

                        tecstat = tecZoneVarWriteDoubleValues(handle, zone_index, int(icoord,c_int32_t), ipartition, int(size(val),c_int64_t), val)
                        if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_domains: Error in call to tecZoneVarWriteDoubleValues")


                    end do ! coords
                end do !ielem



                ! For each variable in equation set, compute value pointwise and save
                do ielem = 1,nelem

                    ! Update location
                    elem_info = data%mesh%get_element_info(idom,ielem)
                    call worker%set_element(elem_info)

                    ! Update the element cache
                    call cache_handler%update(worker,data%eqnset, data%bc_state_group, components    = 'element',   &
                                                                                       face          = NO_ID,       &
                                                                                       differentiate = NO_DIFF,     &
                                                                                       lift          = .false.)

                    ! Retrieve name of current field, retrieve interpolation, write interpolation to file
                    do ifield = 1,data%eqnset(eqn_ID)%prop%nio_fields()
                        var_string = data%eqnset(eqn_ID)%prop%get_io_field_name(ifield)
                        var = worker%get_field(var_string, 'value', 'element')
                        tecstat = tecZoneVarWriteDoubleValues(handle, zone_index, int(3+ifield,c_int32_t), ipartition, int(size(var),c_int64_t), real(var(:)%x_ad_,rdouble))
                        if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_domains: Error in call to tecZoneVarWriteDoubleValues")
                    end do ! ifield

                end do ! ielem


                !
                ! Write element connectivity
                !
                npts_element     = (OUTPUT_RES+1)*(OUTPUT_RES+1)*(OUTPUT_RES+1)
                nsub_per_element = (OUTPUT_RES*OUTPUT_RES*OUTPUT_RES)
                nsub_elements    = nsub_per_element * nelem

                if (allocated(connectivity)) deallocate(connectivity)
                allocate(connectivity(8,nsub_elements), stat=ierr)
                if (ierr /= 0) call AllocationError

                nelem_zeta = OUTPUT_RES
                nelem_eta  = OUTPUT_RES
                nelem_xi   = OUTPUT_RES


                ! For each element
                do ielem = 1,nelem

                    ! Within each element, sample into sub-elements to account for variation within element
                    do ielem_zeta = 1,nelem_zeta
                        do ielem_eta = 1,nelem_eta
                            do ielem_xi = 1,nelem_xi

                                ielem_start  = ielem_xi  +  (ielem_eta-1)*(nelem_xi)  +  (ielem_zeta-1)*(nelem_xi*nelem_eta)
                                ielem_global = ielem_start   +   (ielem-1)*(nelem_zeta*nelem_eta*nelem_xi)
                                ielem_offset = (npts_element)*(ielem-1)
                                istart       = ielem_xi  +  (ielem_eta-1)*(nelem_xi+1)  +  (ielem_zeta-1)*((nelem_xi+1)*(nelem_eta+1))

                                connectivity(1,ielem_global) = istart                                                               +  ielem_offset
                                connectivity(2,ielem_global) = istart+1                                                             +  ielem_offset
                                connectivity(3,ielem_global) = istart + (OUTPUT_RES+1) + 1                                          +  ielem_offset
                                connectivity(4,ielem_global) = istart + (OUTPUT_RES+1)                                              +  ielem_offset
                                connectivity(5,ielem_global) = istart                         +  (OUTPUT_RES+1)*(OUTPUT_RES+1)      +  ielem_offset
                                connectivity(6,ielem_global) = istart+1                       +  (OUTPUT_RES+1)*(OUTPUT_RES+1)      +  ielem_offset
                                connectivity(7,ielem_global) = istart + (OUTPUT_RES+1) + 1    +  (OUTPUT_RES+1)*(OUTPUT_RES+1)      +  ielem_offset
                                connectivity(8,ielem_global) = istart + (OUTPUT_RES+1)        +  (OUTPUT_RES+1)*(OUTPUT_RES+1)      +  ielem_offset

                            end do !ielem_xi
                        end do !ielem_eta
                    end do !ielem_zeta

                end do ! ielem

                tecstat = tecZoneNodeMapWrite32(handle,zone_index,ipartition,int(1,c_int32_t),int(size(connectivity),c_int64_t),reshape(int(connectivity,c_int32_t),[size(connectivity)]))
                if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_domains: Error in call to tecZoneNodeMapWrite32")

            end do ! idom
        end do ! itime

    end subroutine write_tecio_domains
    !***********************************************************************************









    !>  Write tecio surfaces.
    !!
    !!  NOTE: This routine assumes that init_tecio_file has already been called.
    !!
    !!  For each boundary condition patch group in chidg_data, write a surface zone 
    !!  containing the solution for that surface.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/24/2017
    !!
    !-----------------------------------------------------------------------------------
    subroutine write_tecio_surfaces(handle,data,eqn_ID,ndomains_g)
        type(c_ptr),        intent(in)      :: handle
        type(chidg_data_t), intent(inout)   :: data
        integer(ik),        intent(in)      :: eqn_ID
        integer(ik),        intent(in)      :: ndomains_g


        integer(ik)        :: nelem_xi, nelem_eta, nelem_zeta,  &
                              ielem_xi, ielem_eta, ielem_zeta,  &
                              ielem, ielem_global, ielem_offset,&
                              npts_element, nsub_per_element,   &
                              nsub_elements, ierr,              &
                              nelem, istart, ielem_start,       &
                              idom, isurface, itime, icoord,    &
                              nfaces, ibc_face, current_face,   &
                              iface, ipatch, inode, ifield, nsurfaces_g, isurface_g, ntime

        type(AD_D),     allocatable                 :: var(:)
        integer(TEC),   allocatable                 :: connectivity(:,:)
        real(rdouble),  allocatable, dimension(:)   :: val, r, theta, z

    
        integer(c_int32_t) :: ipartition, zone_index, tecstat, iproc, numvars


        integer(ik)                 :: bc_ID
        character(:),   allocatable :: zone_string, var_string

        type(chidg_worker_t)        :: worker
        type(chidg_cache_t)         :: cache
        type(cache_handler_t)       :: cache_handler
        type(element_info_t)        :: elem_info

        type(svector_t)     :: surface_names

        call write_line("   TecIO: Writing surfaces...",io_proc=GLOBAL_MASTER)

        !
        ! Initialize Chidg Worker references
        !
        call worker%init(data%mesh, data%eqnset(:)%prop, data%sdata, data%time_manager, cache)

        numvars = 3 + data%eqnset(eqn_ID)%prop%nio_fields()

        !
        ! Loop time instances
        !
        ntime = data%sdata%q_out%get_ntime()
        do itime = 1,ntime
            worker%itime = itime

            ! Initialize another set of zones for the next time instance
            surface_names = init_tecio_surface_zones(handle,data,eqn_ID,itime)
            nsurfaces_g = surface_names%size()

            ! Loop surfaces
            do isurface = 1,data%mesh%nbc_patch_groups()

                ! Check if valid IO surface
                nfaces      = data%mesh%bc_patch_group(isurface)%nfaces()
                zone_string = data%mesh%bc_patch_group(isurface)%name
                isurface_g  = surface_names%loc(trim(zone_string))

                if ( (nfaces /= 0) .and. (trim(zone_string) /= 'empty') ) then


                    !! Initialize new domain partition in the TecIO file
                    bc_ID = data%get_bc_state_group_id(trim(zone_string))
                    zone_index = int(isurface_g + (itime-1)*nsurfaces_g + ntime*ndomains_g,c_int32_t)
                    ipartition = init_tecio_surface_partition(handle,data%mesh%bc_patch_group(isurface),data%bc_state_group(bc_ID), zone_index)


                    ! For each coordinate, compute it's value pointwise and save
                    ! For each face in each patch, create a sub-sampling of faces to resolve solution variation
                    do ipatch = 1,data%mesh%bc_patch_group(isurface)%npatches()
                        do ibc_face = 1,data%mesh%bc_patch_group(isurface)%patch(ipatch)%nfaces()
                            do icoord = 1,3
                                
                                ! Get face index from the current element
                                idom  = data%mesh%bc_patch_group(isurface)%patch(ipatch)%idomain_l()
                                ielem = data%mesh%bc_patch_group(isurface)%patch(ipatch)%ielement_l_%at(ibc_face)
                                iface = data%mesh%bc_patch_group(isurface)%patch(ipatch)%iface_%at(ibc_face)

                                ! Get coordinate value at point
                                if ( data%mesh%domain(idom)%elems(ielem)%coordinate_system == CYLINDRICAL ) then
                                    r     = real(data%mesh%domain(idom)%faces(ielem,iface)%interp_coords_def(:,1),rdouble)
                                    theta = real(data%mesh%domain(idom)%faces(ielem,iface)%interp_coords_def(:,2),rdouble)
                                    z     = real(data%mesh%domain(idom)%faces(ielem,iface)%interp_coords_def(:,3),rdouble)
                                    if (icoord == 1) val = r*cos(theta)
                                    if (icoord == 2) val = r*sin(theta)
                                    if (icoord == 3) val = z

                                else
                                    if (allocated(val)) deallocate(val)
                                    allocate(val(size(data%mesh%domain(idom)%faces(ielem,iface)%interp_coords_def,1)))
                                    val = real(data%mesh%domain(idom)%faces(ielem,iface)%interp_coords_def(:,icoord),rdouble)
                                end if

                                tecstat = tecZoneVarWriteDoubleValues(handle, zone_index, icoord, ipartition, int(size(val),c_int64_t), val)
                                if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_domains: Error in call to tecZoneVarWriteDoubleValues")

                            end do ! coords
                        end do !ibc_face
                    end do !ipatch



                    ! For each actual face, create a sub-sampling of faces to resolve solution variation
                    do ipatch = 1,data%mesh%bc_patch_group(isurface)%npatches()
                        do ibc_face = 1,data%mesh%bc_patch_group(isurface)%patch(ipatch)%nfaces()

                            ! Get face index from the current element
                            idom  = data%mesh%bc_patch_group(isurface)%patch(ipatch)%idomain_l()
                            ielem = data%mesh%bc_patch_group(isurface)%patch(ipatch)%ielement_l_%at(ibc_face)
                            iface = data%mesh%bc_patch_group(isurface)%patch(ipatch)%iface_%at(ibc_face)

                            ! Update location
                            elem_info = data%mesh%get_element_info(idom,ielem)

                            call worker%set_element(elem_info)
                            call worker%set_face(iface)

                            ! Update the element cache
                            call cache_handler%update(worker,data%eqnset, data%bc_state_group, components    = 'interior faces',    &
                                                                                               face          = iface,               &
                                                                                               differentiate = NO_DIFF,             &
                                                                                               lift          = .false.)

                            ! Retrieve name of current field, retrieve interpolation, write interpolation to file
                            do ifield = 1,data%eqnset(eqn_ID)%prop%nio_fields()
                                var_string = data%eqnset(eqn_ID)%prop%get_io_field_name(ifield)
                                var = worker%get_field(var_string, 'value', 'face interior')

                                tecstat = tecZoneVarWriteDoubleValues(handle, zone_index, 3+ifield, ipartition, int(size(var),c_int64_t), real(var(:)%x_ad_,rdouble))
                                if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_domains: Error in call to tecZoneVarWriteDoubleValues")
                            end do ! ifield

                        end do !ibc_face
                    end do !ipatch




                    !
                    ! Write element connectivity
                    !
                    npts_element     = (OUTPUT_RES+1)*(OUTPUT_RES+1)
                    nsub_per_element = (OUTPUT_RES*OUTPUT_RES)
                    nsub_elements    = nsub_per_element * nfaces

                    if (allocated(connectivity)) deallocate(connectivity)
                    !allocate(connectivity(4,nsub_elements), stat=ierr)
                    allocate(connectivity(8,nsub_elements), stat=ierr)
                    if (ierr /= 0) call AllocationError

                    nelem_eta  = OUTPUT_RES
                    nelem_xi   = OUTPUT_RES


                    ! For each patch/face
                    current_face = 1
                    do ipatch = 1,data%mesh%bc_patch_group(isurface)%npatches()
                        do ibc_face = 1,data%mesh%bc_patch_group(isurface)%patch(ipatch)%nfaces()

                            ! Within each element, sample into sub-elements to account for variation within element
                            do ielem_eta = 1,nelem_eta
                                do ielem_xi = 1,nelem_xi

                                    ielem_start  = ielem_xi  +  (ielem_eta-1)*(nelem_xi)
                                    ielem_global = ielem_start   +   (current_face-1)*(nelem_eta*nelem_xi)
                                    ielem_offset = (npts_element)*(current_face-1)
                                    istart       = ielem_xi  +  (ielem_eta-1)*(nelem_xi+1)

                                    connectivity(1,ielem_global) = istart + ielem_offset                     
                                    connectivity(2,ielem_global) = istart + ielem_offset + 1 
                                    connectivity(3,ielem_global) = istart + ielem_offset + 1 + (OUTPUT_RES+1)
                                    connectivity(4,ielem_global) = istart + ielem_offset     + (OUTPUT_RES+1)     
                                    ! Could duplicate connectivity here if we wanted to use FEBRICK for parallel IO
                                    connectivity(5,ielem_global) = istart + ielem_offset                     
                                    connectivity(6,ielem_global) = istart + ielem_offset + 1 
                                    connectivity(7,ielem_global) = istart + ielem_offset + 1 + (OUTPUT_RES+1)
                                    connectivity(8,ielem_global) = istart + ielem_offset     + (OUTPUT_RES+1)     

                                end do
                            end do
                        
                            current_face = current_face + 1

                        end do !ibc_face
                    end do ! ipatch

                    tecstat = tecZoneNodeMapWrite32(handle,zone_index,ipartition,1,int(size(connectivity),c_int64_t),reshape(connectivity,[size(connectivity)]))
                    if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_domains: Error in call to tecZoneNodeMapWrite32")
                
                end if ! nfaces /= 0

            end do ! isurface
        end do ! itime

    end subroutine write_tecio_surfaces
    !***********************************************************************************









    !>  Initialize TecIO file to write to.
    !!
    !!  This opens a new tecplot binary file for writing and initializes the type, 
    !!  and number of data fields to accept
    !!
    !!  Zeroes strandID counter for zone entities.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  title       Title of the dataset.
    !!  @param[in]  variables   Comma-separated list of variables that will be written.
    !!  @param[in]  filename    Name of the file to be written.
    !!  @param[in]  filetype    Indicating grid or solution file.[0 = full(grid+solution)]
    !!
    !-----------------------------------------------------------------------------------------
    function init_tecio_file(title,variables,file_name,file_type) result(handle)
        character(*)    :: title
        character(*)    :: variables
        character(*)    :: file_name
        integer(TEC)    :: file_type

        integer(c_int32_t)  :: tecstat
        integer(TEC)        :: file_format  = 1      ! 0 = .plt         1 = .szplt
        integer(TEC)        :: data_type    = 1      ! 0=double  1=single  2=32bit int   3=16bit int
        
        type(c_ptr) :: handle
        type(c_ptr) :: empty_handle = C_NULL_PTR


        tecstat = tecFileWriterOpen(trim(file_name)//C_NULL_CHAR,   &
                                    trim(title)//C_NULL_CHAR,       &
                                    trim(variables)//C_NULL_CHAR,   &
                                    file_format,                    &
                                    file_type,                      &
                                    data_type,                      &
                                    empty_handle,                   &
                                    handle)
        if (tecstat /= 0) call chidg_signal(FATAL,"init_tecio_file: Error in TecIO file initialization.")


        tecstat = tecMPIInitialize(handle,ChiDG_COMM%mpi_val,GLOBAL_MASTER)
        if (tecstat /= 0) call chidg_signal(FATAL,"init_tecio_file: Error in TecIO MPI initialization.")


        ! Reset strandID count. Gets incremented every time a new zone is added.
        strandID = 0

    end function init_tecio_file
    !*****************************************************************************************










    !>
    !!
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function init_tecio_volume_zones(handle,data,eqn_ID,itime) result(ndomains_g)
        type(c_ptr),        intent(in)  :: handle
        type(chidg_data_t), intent(in)  :: data
        integer(ik),        intent(in)  :: eqn_ID
        integer(ik),        intent(in)  :: itime

        integer(c_int32_t)  :: zonetype         = 5    ! 5 = FEBRICK
        integer(c_int32_t)  :: fnmode           = 0
        integer(c_int32_t)  :: sharconnfrom     = 0
        integer(c_int64_t)  :: nfconns          = 0

        integer(TEC),   allocatable :: datatypes(:)
        integer(TEC),   allocatable :: locations(:)
        integer(TEC),   allocatable :: sharevars(:)
        integer(TEC),   allocatable :: passivevars(:)

        ! These should be declared with (ik) since we are using MPI_INTEGER4
        ! to communicate. We can convert afterwards to c_int if necessary.
        integer(ik) :: nnodes_g, nnodes_g_reduced
        integer(ik) :: nelements_g, nelements_g_reduced
        integer(ik) :: domain_master, domain_master_reduced
        integer(ik) :: ndomain_procs

        ! Partition Data
        integer(c_int32_t)              :: zoneindex, tecstat, ipartition, npartitions
        integer(c_int64_t)              :: nghost_nodes, nghost_elements
        integer(c_int32_t), allocatable :: ghost_nodes(:), neighbor_nodes(:), ghost_elements(:)
        integer(c_int32_t), allocatable :: neighbor_partitions(:), partition_ranks(:)

        integer(ik)             :: ierr, idomain_g_local_max, idomain_g_global_max, idom,   &
                                   idomain_g, nvars, ndomains_g
        integer(ik), allocatable :: domain_procs(:)
        character(100)          :: domain_name
        real(rk) :: solution_time


        !
        ! Set up some auxiliary TecIO information specific to the way we use the library
        !
        nvars = 3 + data%eqnset(eqn_ID)%prop%nio_fields()
        allocate(datatypes(nvars), sharevars(nvars), locations(nvars), passivevars(nvars), stat=ierr)
        if (ierr /= 0) call AllocationError
        datatypes   = 1
        locations   = 1
        sharevars   = 0
        passivevars = 0


        !
        ! Get maximum idomain_g
        !

        ! Get proc-local maximum idomain_g
        idomain_g_local_max = 0
        do idom = 1,data%mesh%ndomains()
            if (data%mesh%domain(idom)%idomain_g > idomain_g_local_max) idomain_g_local_max = data%mesh%domain(idom)%idomain_g
        end do

        ! Synchronize with all ranks
        call MPI_AllReduce(idomain_g_local_max,idomain_g_global_max,1,MPI_INTEGER4,MPI_MAX,ChiDG_COMM,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'mod_tecio%init_tecio_volume_zones: Error reducing max global domain index.')

        !
        ! SET FUNCTION RETURN VALUE
        !
        ndomains_g = idomain_g_global_max


        !
        ! Need to know nnodes_g, nelements_g for all domains
        !
        do idomain_g = 1,ndomains_g


            ! Handle time index
            strandID = strandID + 1

            ! Does current rank have idomain_g
            nnodes_g      = 0
            nelements_g   = 0
            domain_master = NO_PROC
            ndomain_procs = 0
            do idom = 1,data%mesh%ndomains()
                if (data%mesh%domain(idom)%idomain_g == idomain_g) then
                    nnodes_g    = (OUTPUT_RES+1)*(OUTPUT_RES+1)*(OUTPUT_RES+1) * data%mesh%domain(idom)%get_nelements_global()
                    nelements_g = (OUTPUT_RES*OUTPUT_RES*OUTPUT_RES)           * data%mesh%domain(idom)%get_nelements_global()
                    domain_master = data%mesh%domain(idom)%procs(1)
                    domain_name   = 'Domain '//trim(data%mesh%domain(idom)%name)
                    domain_procs  = data%mesh%domain(idom)%procs
                    ndomain_procs = size(domain_procs)
                    exit
                end if
            end do !idom

            
            ! Synchronize across all ranks. 
            ! NOTE: be sure not to mix different c_int datatypes with MPI_INTEGER4
            call MPI_AllReduce(nnodes_g,nnodes_g_reduced,1,MPI_INTEGER4,MPI_MAX,ChiDG_COMM,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'mod_tecio%init_tecio_volume_zones: Error reducing nnodes_g.')
            call MPI_AllReduce(nelements_g,nelements_g_reduced,1,MPI_INTEGER4,MPI_MAX,ChiDG_COMM,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'mod_tecio%init_tecio_volume_zones: Error reducing nelements_g.')
            call MPI_AllReduce(domain_master,domain_master_reduced,1,MPI_INTEGER4,MPI_MAX,ChiDG_COMM,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'mod_tecio%init_tecio_volume_zones: Error reducing domain_master.')



            ! Broadcast metadata from domain_master rank
            call MPI_BCast(domain_name,len(domain_name),MPI_CHARACTER,domain_master_reduced,ChiDG_COMM,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'mod_tecio%init_tecio_volume_zones: Error broadcasting domain_name.')


            call MPI_BCast(ndomain_procs,1,MPI_INTEGER4,domain_master_reduced,ChiDG_COMM,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'mod_tecio%init_tecio_volume_zones: Error broadcasting ndomain_procs.')


            if (IRANK /= domain_master_reduced) then
                if (allocated(domain_procs)) deallocate(domain_procs)
                allocate(domain_procs(ndomain_procs), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if
            
            call MPI_BCast(domain_procs,ndomain_procs,MPI_INTEGER4,domain_master_reduced,ChiDG_COMM,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'mod_tecio%init_tecio_volume_zones: Error broadcasting domain_procs.')


            ! Create zone
            tecstat = tecZoneCreateFE(handle,                               &
                                      trim(domain_name)//char(0),           &
                                      int(zonetype,c_int32_t),              &
                                      int(nnodes_g_reduced,c_int64_t),      &
                                      int(nelements_g_reduced,c_int64_t),   &
                                      datatypes,                            &
                                      sharevars,                            &
                                      locations,                            &
                                      passivevars,                          &
                                      sharconnfrom,                         &
                                      nfconns,                              &
                                      fnmode,                               &
                                      zoneindex)
            if(tecstat /= 0) call chidg_signal(FATAL,"mod_tecio%init_tecio_volume_zones: Error in TecIO zone initialization.")


            ! Associate MPI Ranks to partitions
            npartitions     = size(domain_procs)
            partition_ranks = int(domain_procs,c_int32_t)

            tecstat = tecZoneMapPartitionsToMPIRanks(handle, zoneindex, npartitions, partition_ranks)
            if(tecstat /= 0) call chidg_signal(FATAL,"mod_tecio%init_tecio_volume_zone: Error in TecIO zone partition to MPI map.")


            ! Set Unsteady StrandID informcation
            solution_time = data%time_manager%times(itime)
            tecstat = tecZoneSetUnsteadyOptions(handle,zoneindex,real(solution_time,rdouble), StrandID)

        end do !idomain_g

    end function init_tecio_volume_zones
    !****************************************************************************************







    !> Create a partition that will be written to by the current MPI rank.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/23/2019
    !!
    !-----------------------------------------------------------------------------------------
    function init_tecio_volume_partition(handle,domain,zone_index) result(ipartition)
        type(c_ptr),    intent(in)  :: handle
        type(domain_t), intent(in)  :: domain
        integer(ik),    intent(in)  :: zone_index

        integer(c_int64_t)  :: nnodes_l
        integer(c_int64_t)  :: nelements_l

        ! Partition Data
        integer(c_int32_t)              :: ipartition, npartitions, iproc
        integer(c_int64_t)              :: nghost_nodes, nghost_elements
        integer(c_int32_t), allocatable :: ghost_nodes(:), neighbor_nodes(:), ghost_elements(:)
        integer(c_int32_t), allocatable :: neighbor_partitions(:), partition_ranks(:)

        integer(c_int32_t)      :: tecstat
        integer(ik)             :: ierr

        ! Get proc-local description for creating partition
        nnodes_l    = (OUTPUT_RES+1)*(OUTPUT_RES+1)*(OUTPUT_RES+1) * domain%get_nelements_local()
        nelements_l = (OUTPUT_RES*OUTPUT_RES*OUTPUT_RES)           * domain%get_nelements_local()


        ! Get partition index associated with current rank
        ipartition = 0
        do iproc = 1,size(domain%procs)
            if (domain%procs(iproc) == IRANK) then
                ipartition = iproc
                exit
            end if
        end do
        if (ipartition < 1 .or. ipartition > size(domain%procs)) call chidg_signal(FATAL,'init_tecio_volume_partition: partition index not found in domain process list.')


        ! Set ghost node/element data particular to how we use TecIO
        nghost_nodes    = 0
        nghost_elements = 0
        allocate(ghost_nodes(0), neighbor_nodes(0), ghost_elements(0), neighbor_partitions(0), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Initial tests indicate PartitionCreate should not be called if there is only a single process
        npartitions = size(domain%procs)
        if (npartitions > 1) then
            tecstat = tecFEPartitionCreate32(handle, int(zone_index,c_int32_t), ipartition, nnodes_l, nelements_l, &
                                             nghost_nodes, ghost_nodes, neighbor_partitions, neighbor_nodes, nghost_elements, ghost_elements)
            if(tecstat /= 0) call chidg_signal(FATAL,"init_tecio_volume_partition: Error in TecIO zone partition creation.")
        else
            ! Return ipartition=0 if zone not partitioned
            ipartition = 0
        end if


    end function init_tecio_volume_partition
    !****************************************************************************************

















    !>  This begins a new surface zone in the current opened file. Must be called
    !!  after init_tecplot_file, because it needs an open binary file.
    !!  If multiple files are open, you can switch between them with the
    !!  TECFIL142 call, as long as the they can be identified by integer values
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/6/2017    Updated to new TecIO API for .szplt
    !!
    !!  @param[in]  zonetitle   Name for the zone being initialized.
    !!  @param[in]  domain      domain_t containing the domain description to be initialized.
    !!  @param[in]  timeindex   Integer index of time strand.
    !!
    !-----------------------------------------------------------------------------------------
    function init_tecio_surface_zones(handle,data,eqn_ID,itime) result(global_names)
        type(c_ptr),        intent(in)  :: handle
        type(chidg_data_t), intent(in)  :: data
        integer(ik),        intent(in)  :: eqn_ID
        integer(ik),        intent(in)  :: itime

        integer(ik)         :: ierr
        integer(TEC)        :: zoneindex
        integer(TEC)        :: tecstat
        !integer(TEC)        :: zonetype         = 3    ! 3 = FEQUADRILATERAL
        integer(TEC)        :: zonetype         = 5    ! 5 = FEBRICK
        integer(TEC)        :: fnmode           = 0
        integer(TEC)        :: sharconnfrom     = 0
        integer(c_int64_t)  :: nfconns          = 0


        ! Partition Data
        integer(c_int32_t)              :: ipartition, npartitions
        integer(c_int64_t)              :: nghost_nodes, nghost_elements
        integer(c_int32_t), allocatable :: ghost_nodes(:), neighbor_nodes(:), ghost_elements(:)
        integer(c_int32_t), allocatable :: neighbor_partitions(:), partition_ranks(:)

        ! These should be declared with (ik) since we are using MPI_INTEGER4
        ! to communicate. We can convert afterwards to c_int if necessary.
        integer(ik) :: nnodes_g, nnodes_g_reduced
        integer(ik) :: nelements_g, nelements_g_reduced
        integer(ik) :: group_ID, bc_ID, nnames, iname, name_size, isurface_g, &
                       nsurfaces_g, nvars, nsurface_procs, surface_master, iproc

        integer(TEC),   allocatable :: datatypes(:)
        integer(TEC),   allocatable :: locations(:)
        integer(TEC),   allocatable :: sharevars(:)
        integer(TEC),   allocatable :: passivevars(:)

        type(svector_t)             :: group_names, global_names
        type(string_t)              :: group_name_string
        character(:),   allocatable :: group_name
        character(1024)             :: group_name_buffer
        logical                     :: ranks_with_surface(NRANK), ranks_with_surface_reduced(NRANK)

        type(mpi_comm)      :: bc_comm
        type(ivector_t)     :: surface_procs
        real(rk)            :: solution_time

        !
        ! Set up some auxiliary TecIO information specific to the way we use the library
        !
        nvars = 3 + data%eqnset(eqn_ID)%prop%nio_fields()
        allocate(datatypes(nvars), sharevars(nvars), locations(nvars), passivevars(nvars), stat=ierr)
        if (ierr /= 0) call AllocationError
        datatypes   = 1
        locations   = 1
        sharevars   = 0
        passivevars = 0


        !
        ! Assemble ordered list of surfaces
        !
        do iproc = 0,NRANK-1

            if (iproc == IRANK) then

                ! First assemble all local surfaces 
                do group_ID = 1,data%mesh%nbc_patch_groups()
                    if (data%mesh%bc_patch_group(group_ID)%name /= 'empty' .and. &
                        data%mesh%bc_patch_group(group_ID)%nfaces() /= 0) then
                        call group_names%push_back_unique(string_t(trim(data%mesh%bc_patch_group(group_ID)%name)))
                    end if 
                end do !isurface

                ! Send to master
                if (iproc /= GLOBAL_MASTER) then

                    nnames = group_names%size()
                    call MPI_Send(nnames,1,MPI_INTEGER4,GLOBAL_MASTER,0,ChiDG_COMM,ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,'init_tecio_surface_zones: error sending nnames to master rank.')

                    do iname = 1,group_names%size()
                        group_name_string = group_names%at(iname)
                        group_name_buffer = group_name_string%str
                        name_size = len(trim(group_name_string%str))

                        call MPI_Send(name_size,1,MPI_INTEGER4,GLOBAL_MASTER,0,ChiDG_COMM,ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,'init_tecio_surface_zones: error sending name size to master rank.')
                        call MPI_Send(group_name_buffer,1024,MPI_CHARACTER,GLOBAL_MASTER,0,ChiDG_COMM,ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,'init_tecio_surface_zones: error sending group name to master rank.')
                    end do
                end if

            end if

            ! Only master rank receives and only if another process 
            ! is sending in this loop iteration
            if ( (IRANK == GLOBAL_MASTER) .and. (iproc /= GLOBAL_MASTER) ) then

                ! Get number of incoming names
                call MPI_Recv(nnames,1,MPI_INTEGER4,iproc,0,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)
                if (ierr /= 0) call chidg_signal(FATAL,'init_tecio_surface_zones: error receiving name size.')

                ! Receive each surface name and only store unique
                do iname = 1,nnames
                    call MPI_Recv(name_size,1,MPI_INTEGER4,iproc,0,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,'init_tecio_surface_zones: error receiving name size.')

                    if (allocated(group_name)) deallocate(group_name)
                    allocate(character(len=name_size) :: group_name,stat=ierr)
                    if (ierr /= 0) call AllocationError

                    !call MPI_Recv(group_name,name_size,MPI_CHARACTER,iproc,0,ChiDG_COMM,ierr)
                    call MPI_Recv(group_name_buffer,1024,MPI_CHARACTER,iproc,0,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)
                    if (ierr /= 0) call chidg_signal(FATAL,'init_tecio_surface_zones: error receiving name size.')

                    call group_names%push_back_unique(string_t(trim(group_name_buffer)))
                end do !nnames
            end if

            call MPI_Barrier(ChiDG_COMM,ierr)

        end do



        !
        ! Distribute group_names to all ranks
        !
        if (IRANK == GLOBAL_MASTER) nnames = group_names%size()

        call MPI_BCast(nnames,1,MPI_INTEGER4,GLOBAL_MASTER,ChiDG_COMM,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'init_tecio_surface_zones: error broadcasting number of global surface names.')


        do iname = 1,nnames


            if (IRANK == GLOBAL_MASTER) then
                group_name_string = group_names%at(iname)
                group_name = group_name_string%str
                name_size = len(trim(group_name))
            end if


            call MPI_BCast(name_size,1,MPI_INTEGER4,GLOBAL_MASTER,ChiDG_COMM,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'init_tecio_surface_zones: error broadcasting name size.')


            if (IRANK /= GLOBAL_MASTER) then
                if (allocated(group_name)) deallocate(group_name)
                allocate(character(len=name_size) :: group_name, stat=ierr)
                if (ierr /= 0) call AllocationError
            end if


            call MPI_BCast(group_name,name_size,MPI_CHARACTER,GLOBAL_MASTER,ChiDG_COMM,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'init_tecio_surface_zones: error broadcasting name.')


            if (IRANK /= GLOBAL_MASTER) call global_names%push_back_unique(string_t(trim(group_name)))


        end do


        ! The master rank has collected all the group names from everyone, so
        ! global_names is group_names
        if (IRANK == GLOBAL_MASTER) global_names = group_names


        !
        ! Need to know nnodes_g, nelements_g for all domains
        !
        do isurface_g = 1,global_names%size()

            ! Handle time index
            strandID = strandID + 1

            ! Clear variables for current surface
            nnodes_g           = 0
            nelements_g        = 0
            surface_master     = NO_PROC
            ranks_with_surface = .false.

            ! Does current rank have isurface_g
            group_name_string = global_names%at(isurface_g)
            group_name = group_name_string%str
            group_ID = data%mesh%get_bc_patch_group_id(trim(group_name))
            if ( group_ID /= NO_ID ) then
                if (data%mesh%bc_patch_group(group_ID)%nfaces() /= 0) then
                    ranks_with_surface(IRANK+1) = .true.
                end if
            end if

            ! Reduce registered participating ranks so everyone knows
            ranks_with_surface_reduced = .false.
            call MPI_AllReduce(ranks_with_surface,ranks_with_surface_reduced,NRANK,MPI_LOGICAL,MPI_LOR,ChiDG_COMM,ierr)
                
            ! Accumulate participating ranks and count
            call surface_procs%clear()
            do iproc = 1,NRANK
                if (ranks_with_surface_reduced(iproc)) call surface_procs%push_back(iproc-1)
            end do
            nsurface_procs = surface_procs%size()
            surface_master = surface_procs%at(1)

            ! Try and find surface on local proc
            if (group_ID /= NO_ID) then
                if (data%mesh%bc_patch_group(group_ID)%nfaces() /= 0) then
                    bc_ID   = data%get_bc_state_group_id(trim(group_name))
                    bc_COMM = data%bc_state_group(bc_ID)%bc_COMM
                    nnodes_g    = (OUTPUT_RES+1)*(OUTPUT_RES+1) * data%bc_state_group(bc_ID)%nfaces_g
                    nelements_g = (OUTPUT_RES*OUTPUT_RES)       * data%bc_state_group(bc_ID)%nfaces_g
                end if
            end if

            ! Synchronize across all ranks. 
            ! NOTE: be sure not to mix different c_int datatypes with MPI_INTEGER4
            call MPI_AllReduce(nnodes_g,nnodes_g_reduced,1,MPI_INTEGER4,MPI_MAX,ChiDG_COMM,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'mod_tecio%init_tecio_surface_zones: Error reducing nnodes_g.')
            call MPI_AllReduce(nelements_g,nelements_g_reduced,1,MPI_INTEGER4,MPI_MAX,ChiDG_COMM,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'mod_tecio%init_tecio_surface_zones: Error reducing nelements_g.')

            ! Create zone
            tecstat = tecZoneCreateFE(handle,                               &
                                      trim(group_name)//char(0),            &
                                      int(zonetype,c_int32_t),              &
                                      int(nnodes_g_reduced,c_int64_t),      &
                                      int(nelements_g_reduced,c_int64_t),   &
                                      datatypes,                            &
                                      sharevars,                            &
                                      locations,                            &
                                      passivevars,                          &
                                      sharconnfrom,                         &
                                      nfconns,                              &
                                      fnmode,                               &
                                      zoneindex)
            if(tecstat /= 0) call chidg_signal(FATAL,"mod_tecio%init_surface_zones: Error in TecIO zone initialization.")

            ! Associate MPI Ranks to partitions
            npartitions     = surface_procs%size()
            partition_ranks = int(surface_procs%data(),c_int32_t)

            tecstat = tecZoneMapPartitionsToMPIRanks(handle, zoneindex, npartitions, partition_ranks)
            if(tecstat /= 0) call chidg_signal(FATAL,"mod_tecio%init_surface_zone: Error in TecIO zone partition to MPI map.")

            ! Set Unsteady StrandID informcation
            solution_time = data%time_manager%times(itime)
            tecstat = tecZoneSetUnsteadyOptions(handle,zoneindex,real(solution_time,rdouble), StrandID)

        end do !isurface_g

    end function init_tecio_surface_zones
    !****************************************************************************************





    !> Create a partition that will be written to by the current MPI rank.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/23/2019
    !!
    !-----------------------------------------------------------------------------------------
    function init_tecio_surface_partition(handle,bc_patch_group,bc_state_group,zone_index) result(ipartition)
        type(c_ptr),            intent(in)  :: handle
        type(bc_patch_group_t), intent(in)  :: bc_patch_group
        type(bc_state_group_t), intent(in)  :: bc_state_group
        integer(ik),            intent(in)  :: zone_index

        integer(c_int64_t)  :: nnodes_l
        integer(c_int64_t)  :: nelements_l

        ! Partition Data
        integer(c_int32_t)              :: ipartition, npartitions, iproc
        integer(c_int64_t)              :: nghost_nodes, nghost_elements
        integer(c_int32_t), allocatable :: ghost_nodes(:), neighbor_nodes(:), ghost_elements(:)
        integer(c_int32_t), allocatable :: neighbor_partitions(:), partition_ranks(:)

        integer(c_int32_t)      :: tecstat
        integer(ik)             :: ierr

        ! Get proc-local description for creating partition
        nnodes_l    = (OUTPUT_RES+1)*(OUTPUT_RES+1) * bc_patch_group%get_nfaces_local()
        nelements_l = (OUTPUT_RES*OUTPUT_RES)       * bc_patch_group%get_nfaces_local()


        ! Get partition index associated with current rank
        ipartition = 0
        do iproc = 1,size(bc_state_group%bc_procs)
            if (bc_state_group%bc_procs(iproc) == IRANK) then
                ipartition = iproc
                exit
            end if
        end do
        if (ipartition < 1 .or. ipartition > size(bc_state_group%bc_procs)) call chidg_signal(FATAL,'init_tecio_surface_partition: partition index not found in bc_patch_group process list.')


        ! Set ghost node/element data particular to how we use TecIO
        nghost_nodes    = 0
        nghost_elements = 0
        allocate(ghost_nodes(0), neighbor_nodes(0), ghost_elements(0), neighbor_partitions(0), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Initial tests indicate PartitionCreate should not be called if there is only a single process
        npartitions = size(bc_state_group%bc_procs)
        if (npartitions > 1) then
            tecstat = tecFEPartitionCreate32(handle, int(zone_index,c_int32_t), ipartition, nnodes_l, nelements_l, &
                                             nghost_nodes, ghost_nodes, neighbor_partitions, neighbor_nodes, nghost_elements, ghost_elements)
            if(tecstat /= 0) call chidg_signal(FATAL,"init_tecio_surface_partition: Error in TecIO zone partition creation.")
        else
            ! Return ipartition=0 if zone not partitioned
            ipartition = 0
        end if


    end function init_tecio_surface_partition
    !****************************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine close_tecio_file(handle)
        type(c_ptr),    intent(inout)  :: handle

        integer(c_int32_t) :: tecstat

        tecstat = tecFileWriterClose(handle)
        if (tecstat /= 0) call chidg_signal(FATAL,"close_tecio: Error in TecIO file end.")

    end subroutine close_tecio_file
    !*****************************************************************************************






end module mod_tecio
