module mod_tecio
#include <messenger.h>
    use mod_kinds,              only: rk,ik,rdouble,TEC
    use mod_constants,          only: ONE, HALF, TWO, OUTPUT_RES, XI_MIN, XI_MAX, &
                                      ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use mod_tecio_interface,    only: init_tecio_file, init_tecio_volume_zone, &
                                      init_tecio_surface_zone, finalize_tecio

    use type_chidg_data,        only: chidg_data_t
    implicit none

#include "tecio.f90"

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/7/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine write_tecio(data,filename)
        type(chidg_data_t),     intent(inout)           :: data
        character(*),           intent(in)              :: filename


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
        call write_tecio_domains(data)


        !
        ! Write surface data from 'mesh%domains'
        !
        call write_tecio_surfaces(data)


        !
        ! Close the current TecIO file context
        !
        call finalize_tecio()


    end subroutine write_tecio
    !***********************************************************************************










    !>  Write tecio domains.
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

                nelem = data%mesh%domain(idom)%nelem

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






























end module mod_tecio
