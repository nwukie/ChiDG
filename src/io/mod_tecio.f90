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
!!  write_tecio
!!      : init_tecio_file
!!      : write_tecio_domains
!!          : init_tecio_volume_zone
!!      : write_tecio_surfaces
!!          : init_tecio_surface_zone
!!      : close_tecio_file
!!
!!
!!  @author Nathan A. Wukie
!!
!--------------------------------------------------------------------------------------
module mod_tecio
#include <messenger.h>
    use mod_kinds,              only: rk,ik,rdouble,TEC
    use mod_constants,          only: ONE, HALF, TWO, OUTPUT_RES, XI_MIN, XI_MAX, &
                                      ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, NO_ID, CYLINDRICAL

    use type_chidg_data,        only: chidg_data_t
    use type_domain,            only: domain_t
    use type_bc_patch_group,    only: bc_patch_group_t

    use type_chidg_worker,      only: chidg_worker_t
    use type_chidg_cache,       only: chidg_cache_t
    use type_cache_handler,     only: cache_handler_t
    use type_element_info,      only: element_info_t
    use type_timer,             only: timer_t
    use DNAD_D
    use iso_c_binding
    implicit none

#include "tecio.f90"


    !
    ! This gets incremented every time a new zone is initialized. 
    ! This get zeroed every time write_tecio is called.
    !
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
    subroutine write_tecio(data,filename,write_domains,write_surfaces)
        type(chidg_data_t),     intent(inout)           :: data
        character(*),           intent(in)              :: filename
        logical,                intent(in)              :: write_domains
        logical,                intent(in)              :: write_surfaces


        integer(ik)     :: ierr, ieq, eqn_ID
        type(c_ptr)     :: handle
        character(100)  :: varstring
        type(timer_t)   :: timer

        call timer%start()

        !
        ! Assemble variables string.
        !
        !   Default: Grid coordinates
        !
        ieq = 1
        varstring = "X,Y,Z"
        eqn_ID = data%mesh%domain(1)%eqn_ID
        do while (ieq <= data%eqnset(eqn_ID)%prop%nio_fields())
            varstring = trim(varstring)//","//trim(data%eqnset(eqn_ID)%prop%get_io_field_name(ieq))
            ieq = ieq + 1
        end do

        
        !
        ! Open and initialize TecIO file
        !
        handle = init_tecio_file('solnfile',trim(varstring),filename,0)


        !
        ! Write volume data from 'mesh%domains'
        !
        if (write_domains) call write_tecio_domains(handle,data)


        !
        ! Write surface data from 'mesh%domains'
        !
        if (write_surfaces) call write_tecio_surfaces(handle,data)


        !
        ! Close the current TecIO file context
        !
        call close_tecio_file(handle)

        call timer%stop()
        call timer%report('Time:')

    end subroutine write_tecio
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
    subroutine write_tecio_domains(handle,data)
        type(c_ptr),            intent(in)      :: handle
        type(chidg_data_t),     intent(inout)   :: data


        integer(ik)        :: nelem_xi, nelem_eta, nelem_zeta,  &
                              ielem_xi, ielem_eta, ielem_zeta,  &
                              ielem, ielem_global, ielem_offset,&
                              npts_element, nsub_per_element,   &
                              nsub_elements, ierr,              &
                              nelem, istart, ielem_start,       &
                              idom, itime, icoord, inode, ifield

        integer(TEC)                                :: zone_index
        integer(TEC)                                :: tecstat
        integer(TEC),   allocatable                 :: connectivity(:,:)
        real(rdouble),  allocatable, dimension(:)   :: val, r, theta, z


        integer(ik)                 :: eqn_ID, numvars
        type(AD_D),     allocatable :: var(:)
        character(:),   allocatable :: zone_string, var_string

        type(chidg_worker_t)        :: worker
        type(chidg_cache_t)         :: cache
        type(cache_handler_t)       :: cache_handler
        type(element_info_t)        :: elem_info



        call write_line("   TECIO: Writing domains...")


        !
        ! Initialize Chidg Worker references
        !
        call worker%init(data%mesh, data%eqnset(:)%prop, data%sdata, cache)


        eqn_ID = data%mesh%domain(1)%eqn_ID
        numvars = 3 + data%eqnset(eqn_ID)%prop%nio_fields()
        !io_ID = data%mesh%domain(1)%io_ID
        !numvars = data%io(io_ID)%prop%nio_fields()


        !
        ! Loop time/domains/elements/fields
        !
        do itime = 1,data%sdata%q_out%get_ntime()
            worker%itime = itime
            do idom = 1,data%mesh%ndomains()


                !
                ! Initialize new zone in the TecIO file for the current domain
                !
                zone_string = 'Domain '//data%mesh%domain(idom)%name
                zone_index = init_tecio_volume_zone(handle,zone_string,data%mesh%domain(idom),data%time_manager%times(itime),numvars)

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
                            val = real(data%mesh%domain(idom)%elems(ielem)%interp_coords_def(:,icoord),rdouble)
                        end if

                        tecstat = tecZoneVarWriteDoubleValues(handle, zone_index, icoord, 0, int(size(val),c_int64_t), val)
                        if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_domains: Error in call to tecZoneVarWriteDoubleValues")


                    end do ! coords
                end do !ielem





                ! For each variable in equation set, compute value pointwise and save
                eqn_ID = data%mesh%domain(idom)%eqn_ID
                do ielem = 1,nelem

                    ! Update location
                    elem_info%idomain_g  = data%mesh%domain(idom)%elems(ielem)%idomain_g
                    elem_info%idomain_l  = data%mesh%domain(idom)%elems(ielem)%idomain_l
                    elem_info%ielement_g = data%mesh%domain(idom)%elems(ielem)%ielement_g
                    elem_info%ielement_l = data%mesh%domain(idom)%elems(ielem)%ielement_l
                    call worker%set_element(elem_info)

                    ! Update the element cache
                    call cache_handler%update(worker,data%eqnset, data%bc_state_group, components    = 'element',   &
                                                                                       face          = NO_ID,       &
                                                                                       differentiate = .false.,     &
                                                                                       lift          = .false.)


                    ! Retrieve name of current field, retrieve interpolation, write interpolation to file
                    do ifield = 1,data%eqnset(eqn_ID)%prop%nio_fields()
                        var_string = data%eqnset(eqn_ID)%prop%get_io_field_name(ifield)
                        var = worker%get_field(var_string, 'value', 'element')
                        tecstat = tecZoneVarWriteDoubleValues(handle, zone_index, 3+ifield, 0, int(size(var),c_int64_t), real(var(:)%x_ad_,rdouble))
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

                                istart = ielem_xi  +  (ielem_eta-1)*(nelem_xi+1)  +  (ielem_zeta-1)*((nelem_xi+1)*(nelem_eta+1))

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

                tecstat = tecZoneNodeMapWrite32(handle,zone_index,0,1,int(size(connectivity),c_int64_t),reshape(connectivity,[size(connectivity)]))
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
    subroutine write_tecio_surfaces(handle,data)
        type(c_ptr),        intent(in)      :: handle
        type(chidg_data_t), intent(inout)   :: data


        integer(ik)        :: nelem_xi, nelem_eta, nelem_zeta,  &
                              ielem_xi, ielem_eta, ielem_zeta,  &
                              ielem, ielem_global, ielem_offset,&
                              npts_element, nsub_per_element,   &
                              nsub_elements, ierr,              &
                              nelem, istart, ielem_start,       &
                              idom, isurface, itime, icoord,    &
                              nfaces, ibc_face, current_face,   &
                              iface, ipatch, inode, ifield

        type(AD_D),     allocatable                 :: var(:)
        integer(TEC),   allocatable                 :: connectivity(:,:)
        real(rdouble),  allocatable, dimension(:)   :: val, r, theta, z

    
        integer(TEC)                :: zone_index
        integer(TEC)                :: tecstat, numvars


        integer(ik)                 :: eqn_ID
        character(:),   allocatable :: zone_string, var_string

        type(chidg_worker_t)        :: worker
        type(chidg_cache_t)         :: cache
        type(cache_handler_t)       :: cache_handler
        type(element_info_t)        :: elem_info




        call write_line("   TECIO: Writing surfaces...")

        !
        ! Initialize Chidg Worker references
        !
        call worker%init(data%mesh, data%eqnset(:)%prop, data%sdata, cache)



        eqn_ID = data%mesh%domain(1)%eqn_ID
        numvars = 3 + data%eqnset(eqn_ID)%prop%nio_fields()

        !
        ! Loop time instances
        !
        do itime = 1,data%sdata%q_out%get_ntime()


            worker%itime = itime


            !
            ! Loop surfaces
            !
            do isurface = 1,data%mesh%nbc_patch_groups()

                !
                ! Get number of surface elements
                !
                nfaces = data%mesh%bc_patch_group(isurface)%nfaces()


                !
                ! Initialize new zone in the TecIO file for the current domain
                !
                zone_string = data%mesh%bc_patch_group(isurface)%name
                zone_index = init_tecio_surface_zone(handle,zone_string,data%mesh%bc_patch_group(isurface),data%time_manager%times(itime),numvars)



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
                                val = real(data%mesh%domain(idom)%faces(ielem,iface)%interp_coords_def(:,icoord),rdouble)
                            end if

                            tecstat = tecZoneVarWriteDoubleValues(handle, zone_index, icoord, 0, int(size(val),c_int64_t), val)
                            if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_domains: Error in call to tecZoneVarWriteDoubleValues")

                        end do ! coords
                    end do !ibc_face
                end do !ipatch





                ! For each variable in equation set, compute value pointwise and save
                eqn_ID = data%mesh%domain(idom)%eqn_ID

                ! For each actual face, create a sub-sampling of faces to resolve solution variation
                do ipatch = 1,data%mesh%bc_patch_group(isurface)%npatches()
                    do ibc_face = 1,data%mesh%bc_patch_group(isurface)%patch(ipatch)%nfaces()

                        !
                        ! Get face index from the current element
                        !
                        idom  = data%mesh%bc_patch_group(isurface)%patch(ipatch)%idomain_l()
                        ielem = data%mesh%bc_patch_group(isurface)%patch(ipatch)%ielement_l_%at(ibc_face)
                        iface = data%mesh%bc_patch_group(isurface)%patch(ipatch)%iface_%at(ibc_face)


                        ! Update location
                        elem_info%idomain_g  = data%mesh%domain(idom)%elems(ielem)%idomain_g
                        elem_info%idomain_l  = data%mesh%domain(idom)%elems(ielem)%idomain_l
                        elem_info%ielement_g = data%mesh%domain(idom)%elems(ielem)%ielement_g
                        elem_info%ielement_l = data%mesh%domain(idom)%elems(ielem)%ielement_l
                        call worker%set_element(elem_info)
                        call worker%set_face(iface)

                        ! Update the element cache
                        call cache_handler%update(worker,data%eqnset, data%bc_state_group, components    = 'interior faces',   &
                                                                                           face          = iface,   &
                                                                                           differentiate = .false.,  &
                                                                                           lift          = .false.)

                        ! Retrieve name of current field, retrieve interpolation, write interpolation to file
                        do ifield = 1,data%eqnset(eqn_ID)%prop%nio_fields()
                            var_string = data%eqnset(eqn_ID)%prop%get_io_field_name(ifield)
                            var = worker%get_field(var_string, 'value', 'face interior')

                            tecstat = tecZoneVarWriteDoubleValues(handle, zone_index, 3+ifield, 0, int(size(var),c_int64_t), real(var(:)%x_ad_,rdouble))
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
                allocate(connectivity(4,nsub_elements), stat=ierr)
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

                            end do
                        end do
                    
                        current_face = current_face + 1

                    end do !ibc_face
                end do ! ipatch

                tecstat = tecZoneNodeMapWrite32(handle,zone_index,0,1,int(size(connectivity),c_int64_t),reshape(connectivity,[size(connectivity)]))
                if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_domains: Error in call to tecZoneNodeMapWrite32")
                

            end do ! idom

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

        integer(4)          :: tecstat
        character           :: NULLCHAR     = char(0)
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


        ! Reset strandID count. Gets incremented every time a new zone is added.
        strandID = 0

    end function init_tecio_file
    !*****************************************************************************************










    !>  This begins a new volume zone in the current opened file. Must be called
    !!  after init_tecplot_file, because it needs an open binary file.
    !!  If multiple files are open, you can switch between them with the
    !!  TECFIL142 call, as long as the they can be identified by integer values
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  zonetitle   Name for the zone being initialized.
    !!  @param[in]  domain      domain_t containing the domain description to be initialized.
    !!  @param[in]  timeindex   Integer index of time strand.
    !!
    !-----------------------------------------------------------------------------------------
    function init_tecio_volume_zone(handle,zonetitle,domain,solutiontime,numvars) result(zoneindex)
        type(c_ptr),    intent(in)  :: handle
        character(*),   intent(in)  :: zonetitle
        type(domain_t), intent(in)  :: domain
        real(rk),       intent(in)  :: solutiontime
        integer(ik),    intent(in)  :: numvars

        integer(TEC)        :: zonetype         = 5    ! 5 = FEBRICK
        integer(TEC)        :: fnmode           = 0
        integer(TEC)        :: sharconnfrom     = 0
        integer(c_int64_t)  :: nfconns          = 0
        integer(c_int64_t)  :: numnodes
        integer(c_int64_t)  :: numelements

        integer(TEC),   allocatable :: datatypes(:)
        integer(TEC),   allocatable :: locations(:)
        integer(TEC),   allocatable :: sharevars(:)
        integer(TEC),   allocatable :: passivevars(:)


        integer(TEC)            :: zoneindex
        integer(4)              :: tecstat, ierr
        integer(TEC),   pointer :: NullPtr(:) => null()    ! Null pointer array


        !
        ! Handle time index
        !
        strandID = strandID + 1


        !
        ! Handle domain discretization
        !
        numnodes    = (OUTPUT_RES+1)*(OUTPUT_RES+1)*(OUTPUT_RES+1) * domain%nelem
        numelements = (OUTPUT_RES*OUTPUT_RES*OUTPUT_RES) * domain%nelem



        allocate(datatypes(numvars), sharevars(numvars), locations(numvars), passivevars(numvars), stat=ierr)
        if (ierr /= 0) call AllocationError

        datatypes   = 1
        locations   = 1
        sharevars   = 0
        passivevars = 0

        ! Create zone
        tecstat = tecZoneCreateFE(handle,                   &
                                  trim(zonetitle)//char(0), &
                                  zonetype,                 &
                                  numnodes,                 &
                                  numelements,              &
                                  datatypes,                &
                                  sharevars,                &
                                  locations,                &
                                  passivevars,              &
                                  sharconnfrom,             &
                                  nfconns,                  &
                                  fnmode,                   &
                                  zoneindex)

        if(tecstat /= 0) call chidg_signal(FATAL,"init_tecio_volume_zone: Error in TecIO zone initialization.")


        ! Set time 
        tecstat = tecZoneSetUnsteadyOptions(handle,zoneindex,real(solutiontime,rdouble), strandid)
        if(tecstat /= 0) call chidg_signal(FATAL,"init_tecio_volume_zone: Error in TecIO zone initialization.")


    end function init_tecio_volume_zone
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
    function init_tecio_surface_zone(handle,zonetitle,patch_group,solutiontime,numvars) result(zoneindex)
        type(c_ptr),            intent(in)  :: handle
        character(*),           intent(in)  :: zonetitle
        type(bc_patch_group_t), intent(in)  :: patch_group
        real(rk),               intent(in)  :: solutiontime
        integer(ik),            intent(in)  :: numvars

        integer(ik)         :: ierr
        integer(TEC)        :: zoneindex
        integer(TEC)        :: tecstat
        integer(TEC)        :: zonetype         = 3    ! 3 = FEQUADRILATERAL
        integer(TEC)        :: fnmode           = 0
        integer(TEC)        :: sharconnfrom     = 0
        integer(c_int64_t)  :: nfconns          = 0
        integer(c_int64_t)  :: numnodes
        integer(c_int64_t)  :: numelements

        integer(TEC),   allocatable :: datatypes(:)
        integer(TEC),   allocatable :: locations(:)
        integer(TEC),   allocatable :: sharevars(:)
        integer(TEC),   allocatable :: passivevars(:)


        !
        ! Handle time index
        !
        strandID = strandID + 1


        !
        ! Accumulate total number of points to be written on the surface
        !
        numnodes    = (OUTPUT_RES+1)*(OUTPUT_RES+1) * patch_group%nfaces()
        numelements = (OUTPUT_RES*OUTPUT_RES)       * patch_group%nfaces()


        allocate(datatypes(numvars), sharevars(numvars), locations(numvars), passivevars(numvars), stat=ierr)
        if (ierr /= 0) call AllocationError

        datatypes   = 1
        locations   = 1
        sharevars   = 0
        passivevars = 0

        ! Create zone
        tecstat = tecZoneCreateFE(handle,                   &
                                  trim(zonetitle)//char(0), &
                                  zonetype,                 &
                                  numnodes,                 &
                                  numelements,              &
                                  datatypes,                &
                                  sharevars,                &
                                  locations,                &
                                  passivevars,              &
                                  sharconnfrom,             &
                                  nfconns,                  &
                                  fnmode,                   &
                                  zoneindex)

        if(tecstat /= 0) call chidg_signal(FATAL,"init_tecio_surface_zone: Error in TecIO zone initialization.")


        ! Set time 
        tecstat = tecZoneSetUnsteadyOptions(handle,zoneindex,real(solutiontime,rdouble), strandid)
        if(tecstat /= 0) call chidg_signal(FATAL,"init_tecio_surface_zone: Error in TecIO zone initialization.")


    end function init_tecio_surface_zone
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

        integer(TEC) :: tecstat

        tecstat = tecFileWriterClose(handle)
        if (tecstat /= 0) call chidg_signal(FATAL,"close_tecio: Error in TecIO file end.")

    end subroutine close_tecio_file
    !*****************************************************************************************


















end module mod_tecio
