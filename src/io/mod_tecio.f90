!--------------------------------------------------------------------------------------
!!
!!  API for writing tecplot files using the TecIO library.
!!
!!  Right now, these capabilities exist strictly to support post-processing
!!  and visualization of solution data. Modal data is interpolated to 
!!  discrete points and elements/faces are sub-divided to resolve polynomial
!!  variation within these entities.
!!
!!  High-level procedures:
!!  ----------------------
!!  write_tecio
!!
!!
!!  Mid-level procedures:
!!  ---------------------
!!  write_tecio_domains
!!  write_tecio_surfaces
!!
!!
!!  Low-level procedures:
!!  ---------------------
!!  init_tecio_file
!!  init_tecio_volume_zone
!!  init_tecio_surface_zone
!!  close_tecio_file
!!
!!
!--------------------------------------------------------------------------------------
module mod_tecio
#include <messenger.h>
    use mod_kinds,              only: rk,ik,rdouble,TEC
    use mod_constants,          only: ONE, HALF, TWO, OUTPUT_RES, XI_MIN, XI_MAX, &
                                      ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX

    use type_chidg_data,        only: chidg_data_t
    use type_domain,            only: domain_t
    use type_bc_patch_group,    only: bc_patch_group_t
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
    !!  @param[inout]   data            chidg_data instance containing initialized mesh/solution data.
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


        integer(ik)     :: ierr, ieq
        character(100)  :: varstring


        !
        ! Assemble variables string.
        !
        !   Default: Grid coordinates
        !
        varstring = "X,Y,Z"


        !
        ! TODO: Generalized TECIO for different equation set in each domain.
        !
        ieq = 1
        do while (ieq <= data%eqnset(1)%prop%nprimary_fields())
            varstring = trim(varstring)//","//trim(data%eqnset(1)%prop%get_primary_field_name(ieq))
            ieq = ieq + 1
        end do

        
        !
        ! Open and initialize TecIO file
        !
        call init_tecio_file('solnfile',trim(varstring),filename,0)


        !
        ! Write volume data from 'mesh%domains'
        !
        if (write_domains) call write_tecio_domains(data)


        !
        ! Write surface data from 'mesh%domains'
        !
        if (write_surfaces) call write_tecio_surfaces(data)


        !
        ! Close the current TecIO file context
        !
        call close_tecio_file()


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
    subroutine write_tecio_domains(data)
        type(chidg_data_t),     intent(inout)   :: data


        integer(ik)        :: nelem_xi, nelem_eta, nelem_zeta,  &
                              ielem_xi, ielem_eta, ielem_zeta,  &
                              npts_xi,  npts_eta,  npts_zeta,   &
                              ipt_xi,   ipt_eta,   ipt_zeta,    &
                              xilim,    etalim,    zetalim,     &
                              ielem, ielem_global, ielem_offset,&
                              npts_element, nsub_per_element,   &
                              nsub_elements, npts, ierr,        &
                              nelem, istart, ielem_start,       &
                              ivar, idom, itime, icoord

        real(rdouble)      :: val(1), r, theta, z
        real(TEC)          :: valeq(1)
        equivalence           (valeq(1), val(1))

    
        integer(4)                  :: tecstat
        integer(4),     allocatable :: connectivity(:,:)


        real(rk)                    :: xi,eta,zeta,p, sumsqr, d_normalization, grad1_d, grad2_d, grad3_d
        integer(ik)                 :: eqn_ID
        character(:),   allocatable :: zonestring



        ! using (output_res+1) so that the skip number used in tecplot to
        ! correctly display the element surfaces is the same as the number
        ! specified in the input file
        npts = OUTPUT_RES+1



        !
        ! Loop time instances
        !
        do itime = 1,data%sdata%q_out%get_ntime()


            !
            ! Loop domains
            !
            do idom = 1,data%mesh%ndomains()

                nelem    = data%mesh%domain(idom)%nelem

                !
                ! Initialize new zone in the TecIO file for the current domain
                !
                zonestring = 'Domain '//data%mesh%domain(idom)%name
                call init_tecio_volume_zone(zonestring,data%mesh%domain(idom),itime)


                xilim   = npts
                etalim  = npts
                zetalim = npts

                ! For each coordinate, compute it's value pointwise and save
                do icoord = 1,3

                    ! For each actual element, create a sub-sampling of elements to resolve solution variation
                    do ielem = 1,nelem
                        

                        ! Write sampling for current element
                        do ipt_zeta = 1,zetalim
                            zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
                            do ipt_eta = 1,etalim
                                eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
                                do ipt_xi = 1,xilim
                                    xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                    ! Get coordinate value at point
                                    if ( data%mesh%domain(idom)%elems(ielem)%coordinate_system == 'Cylindrical' ) then

                                        r     = real(data%mesh%domain(idom)%elems(ielem)%grid_point(1,xi,eta,zeta),rdouble)
                                        theta = real(data%mesh%domain(idom)%elems(ielem)%grid_point(2,xi,eta,zeta),rdouble)
                                        z     = real(data%mesh%domain(idom)%elems(ielem)%grid_point(3,xi,eta,zeta),rdouble)

                                        if (icoord == 1) then
                                            val = r*cos(theta)
                                        else if (icoord == 2) then
                                            val = r*sin(theta)
                                        else if (icoord == 3) then
                                            val = z
                                        end if

                                    else

                                        val = real(data%mesh%domain(idom)%elems(ielem)%grid_point(icoord,xi,eta,zeta),rdouble)

                                    end if



                                    tecstat = TECDAT142(1,valeq,1)
                                    if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_domains: Error in call to TECDAT142")

                                end do ! ipt_xi
                            end do ! ipt_eta
                        end do ! ipt_zeta



                    end do !ielem


                end do ! coords





                ! For each variable in equation set, compute value pointwise and save
                eqn_ID = data%mesh%domain(idom)%eqn_ID
                do ivar = 1,data%eqnset(eqn_ID)%prop%nprimary_fields()

                    ! For each actual element, create a sub-sampling of elements to resolve solution variation
                    do ielem = 1,nelem
                        
                        ! Write sampling for current element
                        do ipt_zeta = 1,zetalim
                            zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
                            do ipt_eta = 1,etalim
                                eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
                                do ipt_xi = 1,xilim
                                    xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO


                                    !
                                    ! Get solution value at point
                                    !   
                                    val = real(data%mesh%domain(idom)%elems(ielem)%solution_point(data%sdata%q_out%dom(idom)%vecs(ielem),ivar,itime,xi,eta,zeta),rdouble)
                                    tecstat = TECDAT142(1,valeq,1)
                                    if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_domains: Error in call to TECDAT142")
                                        

                                end do
                            end do
                        end do


                    end do

                end do ! ivar










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


                            end do
                        end do
                    end do
                    

                end do ! ielem

                tecstat = TECNOD142(connectivity)
                if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_domains: Error in call to TECNOD142")

                

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
    subroutine write_tecio_surfaces(data)
        type(chidg_data_t), intent(inout)   :: data


        integer(ik)        :: nelem_xi, nelem_eta, nelem_zeta,  &
                              ielem_xi, ielem_eta, ielem_zeta,  &
                              npts_xi,  npts_eta,  npts_zeta,   &
                              ipt_xi,   ipt_eta,   ipt_zeta,    &
                              xilim,    etalim,    zetalim,     &
                              ielem, ielem_global, ielem_offset,&
                              npts_element, nsub_per_element,   &
                              nsub_elements, npts, ierr,        &
                              nelem, istart, ielem_start,       &
                              ivar, idom, isurface, itime, icoord, nfaces, ibc_face, current_face, iface, ipatch

        real(rdouble)      :: val(1), r, theta, z
        real(TEC)          :: valeq(1)
        equivalence           (valeq(1), val(1))

    
        integer(4)                  :: tecstat
        integer(4),     allocatable :: connectivity(:,:)


        real(rk)                    :: xi,eta,zeta,p, sumsqr, d_normalization, grad1_d, grad2_d, grad3_d
        integer(ik)                 :: eqn_ID
        character(:),   allocatable :: zonestring



        ! using (output_res+1) so that the skip number used in tecplot to
        ! correctly display the element surfaces is the same as the number
        ! specified in the input file
        npts = OUTPUT_RES+1



        !
        ! Loop time instances
        !
        do itime = 1,data%sdata%q_out%get_ntime()


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
                zonestring = data%mesh%bc_patch_group(isurface)%name
                call init_tecio_surface_zone(zonestring,data%mesh%bc_patch_group(isurface),itime)



                ! For each coordinate, compute it's value pointwise and save
                do icoord = 1,3

                    ! For each face in each patch, create a sub-sampling of faces to resolve solution variation
                    do ipatch = 1,data%mesh%bc_patch_group(isurface)%npatches()
                        do ibc_face = 1,data%mesh%bc_patch_group(isurface)%patch(ipatch)%nfaces()
                            
                            !
                            ! Get face index from the current element
                            !
                            idom  = data%mesh%bc_patch_group(isurface)%patch(ipatch)%idomain_l_%at(ibc_face)
                            ielem = data%mesh%bc_patch_group(isurface)%patch(ipatch)%ielement_l_%at(ibc_face)
                            iface = data%mesh%bc_patch_group(isurface)%patch(ipatch)%iface_%at(ibc_face)


                            !
                            ! Set parameters based on face
                            ! 
                            select case(iface)
                                case(XI_MIN)
                                    xi      = -ONE
                                    xilim   = 1
                                    etalim  = npts
                                    zetalim = npts

                                case(XI_MAX)
                                    xi      = ONE
                                    xilim   = 1
                                    etalim  = npts
                                    zetalim = npts

                                case(ETA_MIN)
                                    eta     = -ONE
                                    xilim   = npts
                                    etalim  = 1
                                    zetalim = npts
                                    

                                case(ETA_MAX)
                                    eta     = ONE
                                    xilim   = npts
                                    etalim  = 1
                                    zetalim = npts

                                case(ZETA_MIN)
                                    zeta    = -ONE
                                    xilim   = npts
                                    etalim  = npts
                                    zetalim = 1

                                case(ZETA_MAX)
                                    zeta    = ONE
                                    xilim   = npts
                                    etalim  = npts
                                    zetalim = 1

                                case default
                                    call chidg_signal(FATAL,"write_tecio_surfaces: Invalid face index.")

                            end select


                            !
                            ! Write sub-sampling for current element
                            !   Note: don't vary xi/eta/zeta if that is the face we are writing
                            !
                            do ipt_zeta = 1,zetalim
                                if ((iface /= ZETA_MIN) .and. (iface /= ZETA_MAX)) then
                                    zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
                                end if

                                do ipt_eta = 1,etalim
                                    if ((iface /= ETA_MIN) .and. (iface /= ETA_MAX)) then
                                        eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
                                    end if

                                    do ipt_xi = 1,xilim
                                        if ((iface /= XI_MIN) .and. (iface /= XI_MAX)) then
                                            xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
                                        end if

                                        ! Get coordinate value at point
                                        if ( data%mesh%domain(idom)%elems(ielem)%coordinate_system == 'Cylindrical' ) then

                                            r     = real(data%mesh%domain(idom)%elems(ielem)%grid_point(1,xi,eta,zeta),rdouble)
                                            theta = real(data%mesh%domain(idom)%elems(ielem)%grid_point(2,xi,eta,zeta),rdouble)
                                            z     = real(data%mesh%domain(idom)%elems(ielem)%grid_point(3,xi,eta,zeta),rdouble)

                                            if (icoord == 1) then
                                                val = r*cos(theta)
                                            else if (icoord == 2) then
                                                val = r*sin(theta)
                                            else if (icoord == 3) then
                                                val = z
                                            end if

                                        else

                                            val = real(data%mesh%domain(idom)%elems(ielem)%grid_point(icoord,xi,eta,zeta),rdouble)

                                        end if



                                        tecstat = TECDAT142(1,valeq,1)
                                        if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_surfaces: Error in call to TECDAT142")

                                    end do ! ipt_xi
                                end do ! ipt_eta
                            end do ! ipt_zeta



                        end do !ibc_face
                    end do !ipatch


                end do ! coords





                ! For each variable in equation set, compute value pointwise and save
                eqn_ID = data%mesh%domain(idom)%eqn_ID
                do ivar = 1,data%eqnset(eqn_ID)%prop%nprimary_fields()

                    ! For each actual face, create a sub-sampling of faces to resolve solution variation
                    do ipatch = 1,data%mesh%bc_patch_group(isurface)%npatches()
                        do ibc_face = 1,data%mesh%bc_patch_group(isurface)%patch(ipatch)%nfaces()
                            
                            !
                            ! Get face index from the current element
                            !
                            idom  = data%mesh%bc_patch_group(isurface)%patch(ipatch)%idomain_l_%at(ibc_face)
                            ielem = data%mesh%bc_patch_group(isurface)%patch(ipatch)%ielement_l_%at(ibc_face)
                            iface = data%mesh%bc_patch_group(isurface)%patch(ipatch)%iface_%at(ibc_face)


                            !
                            ! Set parameters based on face
                            ! 
                            select case(iface)
                                case(XI_MIN)
                                    xi      = -ONE
                                    xilim   = 1
                                    etalim  = npts
                                    zetalim = npts

                                case(XI_MAX)
                                    xi      = ONE
                                    xilim   = 1
                                    etalim  = npts
                                    zetalim = npts

                                case(ETA_MIN)
                                    eta     = -ONE
                                    xilim   = npts
                                    etalim  = 1
                                    zetalim = npts
                                    

                                case(ETA_MAX)
                                    eta     = ONE
                                    xilim   = npts
                                    etalim  = 1
                                    zetalim = npts

                                case(ZETA_MIN)
                                    zeta    = -ONE
                                    xilim   = npts
                                    etalim  = npts
                                    zetalim = 1

                                case(ZETA_MAX)
                                    zeta    = ONE
                                    xilim   = npts
                                    etalim  = npts
                                    zetalim = 1

                                case default
                                    call chidg_signal(FATAL,"write_tecio_surfaces: Invalid face index.")

                            end select




                            !
                            ! Write sub-sampling for current element
                            !   Note: don't vary xi/eta/zeta if that is the face we are writing
                            !
                            do ipt_zeta = 1,zetalim
                                if ((iface /= ZETA_MIN) .and. (iface /= ZETA_MAX)) then
                                    zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
                                end if

                                do ipt_eta = 1,etalim
                                    if ((iface /= ETA_MIN) .and. (iface /= ETA_MAX)) then
                                        eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
                                    end if

                                    do ipt_xi = 1,xilim
                                        if ((iface /= XI_MIN) .and. (iface /= XI_MAX)) then
                                            xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
                                        end if

                                        !
                                        ! Get solution value at point
                                        !   
                                        val = real(data%mesh%domain(idom)%elems(ielem)%solution_point(data%sdata%q_out%dom(idom)%vecs(ielem),ivar,itime,xi,eta,zeta),rdouble)
                                        tecstat = TECDAT142(1,valeq,1)
                                        if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_surfaces: Error in call to TECDAT142")
                                            

                                    end do
                                end do
                            end do


                        end do !ibc_face
                    end do !ipatch

                end do ! ivar










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

                tecstat = TECNOD142(connectivity)
                if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_domains: Error in call to TECNOD142")

                

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
    subroutine init_tecio_file(title,variables,filename,filetype)
        character(*)    :: title
        character(*)    :: variables
        character(*)    :: filename
        integer(TEC)    :: filetype

        integer(4)      :: tecstat
        character       :: NULLCHAR = char(0)
        integer(TEC)    :: fileformat = 0       ! 0 = .plt         1 = subzone loadable .szplt
        integer(TEC)    :: isdouble   = 1       ! 0 = single prec  1 = double prec
        integer(TEC)    :: debug      = 0       ! 0 = debug off

        tecstat = TECINI142(trim(title)//NULLCHAR,      &
                            trim(variables)//NULLCHAR,  &
                            trim(filename)//NULLCHAR,   &
                            '.'//NULLCHAR,              &
                            fileformat,                 &
                            filetype,                   &
                            debug,                      &
                            isdouble)

        if (tecstat /= 0) call chidg_signal(FATAL,"init_tecio_file: Error in TecIO file initialization.")

        ! Reset strandID count. Gets incremented every time a new zone is added.
        strandID = 0

    end subroutine init_tecio_file
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
    subroutine init_tecio_volume_zone(zonetitle,domain,timeindex)
            character(*),   intent(in)  :: zonetitle
            type(domain_t), intent(in)  :: domain
            integer(ik),    intent(in)  :: timeindex

            integer(TEC)   :: zonetype                  = 5    ! 5 = FEBRICK
            integer(TEC)   :: numpts
            integer(TEC)   :: numelements
            integer(TEC)   :: numfaces                  = 0    ! not used
            integer(TEC)   :: icellmax                  = 0    ! not used
            integer(TEC)   :: jcellmax                  = 0    ! not used
            integer(TEC)   :: kcellmax                  = 0    ! not used
            real(rk)       :: solutiontime              = 0._rk
!            integer(TEC)   :: strandid                  = 0    ! strandID is now a module variable
            integer(TEC)   :: parentzone                = 0
            integer(TEC)   :: isblock                   = 1
            integer(TEC)   :: nfconns                   = 0
            integer(TEC)   :: fnmode                    = 0
            integer(TEC)   :: totalnumfacenodes         = 1
            integer(TEC)   :: totalnumbndryfaces        = 1
            integer(TEC)   :: totalnumbndryconnections  = 1
            integer(TEC)   :: passivevars(3)            = 0    ! null = all vars active
            integer(TEC)   :: vallocation(3)            = 1    ! null = all vars node-centered
            integer(TEC)   :: sharvarfrom(3)            = 0    ! null = zones share no data
            integer(TEC)   :: sharconnfrom              = 0

            integer(4)              :: tecstat
            integer(TEC),   pointer :: NullPtr(:) => null()    ! Null pointer array


            !
            ! Handle time index
            !
            strandID = strandID + 1
            solutiontime = real(timeindex,rk)


            !
            ! Handle domain discretization
            !
            numpts      = (OUTPUT_RES+1)*(OUTPUT_RES+1)*(OUTPUT_RES+1) * domain%nelem
            numelements = (OUTPUT_RES*OUTPUT_RES*OUTPUT_RES) * domain%nelem



            tecstat = TECZNE142(trim(zonetitle)//char(0),   &
                                zonetype,                   &
                                numpts,                     &
                                numelements,                &
                                numfaces,                   &
                                icellmax,                   &
                                jcellmax,                   &
                                kcellmax,                   &
                                real(solutiontime,rdouble), &
                                strandid,                   &
                                parentzone,                 &
                                isblock,                    &
                                nfconns,                    &
                                fnmode,                     &
                                totalnumfacenodes,          &
                                totalnumbndryfaces,         &
                                totalnumbndryconnections,   &
                                NullPtr,                    &
                                NullPtr,                    &
                                NullPtr,                    &
                                sharconnfrom)

            if(tecstat /= 0) call chidg_signal(FATAL,"init_tecio_volume_zone: Error in TecIO zone initialization.")

    end subroutine init_tecio_volume_zone
    !****************************************************************************************








    !>  This begins a new surface zone in the current opened file. Must be called
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
    subroutine init_tecio_surface_zone(zonetitle,patch_group,timeindex)
            character(*),           intent(in)  :: zonetitle
            type(bc_patch_group_t), intent(in)  :: patch_group
            integer(ik),            intent(in)  :: timeindex

            integer(TEC)   :: zonetype                  = 3    ! 3 = FEQUADRILATERAL
            integer(TEC)   :: numpts
            integer(TEC)   :: numelements
            integer(TEC)   :: numfaces                  = 0    ! not used
            integer(TEC)   :: icellmax                  = 0    ! not used
            integer(TEC)   :: jcellmax                  = 0    ! not used
            integer(TEC)   :: kcellmax                  = 0    ! not used
            real(rk)       :: solutiontime              = 0._rk
!            integer(TEC)   :: strandid                  = 0    ! strandID is now a module variable
            integer(TEC)   :: parentzone                = 0
            integer(TEC)   :: isblock                   = 1
            integer(TEC)   :: nfconns                   = 0
            integer(TEC)   :: fnmode                    = 0
            integer(TEC)   :: totalnumfacenodes         = 1
            integer(TEC)   :: totalnumbndryfaces        = 1
            integer(TEC)   :: totalnumbndryconnections  = 1
            integer(TEC)   :: passivevars(3)            = 0    ! null = all vars active
            integer(TEC)   :: vallocation(3)            = 1    ! null = all vars node-centered
            integer(TEC)   :: sharvarfrom(3)            = 0    ! null = zones share no data
            integer(TEC)   :: sharconnfrom              = 0

            integer(4)              :: tecstat
            integer(TEC),   pointer :: NullPtr(:) => null()    ! Null pointer array


            !
            ! Handle time index
            !
            strandID = strandID + 1
            solutiontime = real(timeindex,rk)


            !
            ! Accumulate total number of points to be written on the surface
            !
            numpts      = (OUTPUT_RES+1)*(OUTPUT_RES+1) * patch_group%nfaces()
            numelements = (OUTPUT_RES*OUTPUT_RES)       * patch_group%nfaces()



            tecstat = TECZNE142(trim(zonetitle)//char(0),   &
                                zonetype,                   &
                                numpts,                     &
                                numelements,                &
                                numfaces,                   &
                                icellmax,                   &
                                jcellmax,                   &
                                kcellmax,                   &
                                real(solutiontime,rdouble), &
                                strandid,                   &
                                parentzone,                 &
                                isblock,                    &
                                nfconns,                    &
                                fnmode,                     &
                                totalnumfacenodes,          &
                                totalnumbndryfaces,         &
                                totalnumbndryconnections,   &
                                NullPtr,                    &
                                NullPtr,                    &
                                NullPtr,                    &
                                sharconnfrom)

            if(tecstat /= 0) call chidg_signal(FATAL,"init_tecio_surface_zone: Error in TecIO zone initialization.")

    end subroutine init_tecio_surface_zone
    !****************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine close_tecio_file()
        integer(kind=TEC) :: tecstat

        tecstat = TECEND142()
        if (tecstat /= 0) call chidg_signal(FATAL,"close_tecio: Error in TecIO file end.")

    end subroutine close_tecio_file
    !*****************************************************************************************


















end module mod_tecio
