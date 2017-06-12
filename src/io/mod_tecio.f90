module mod_tecio
#include <messenger.h>
    use mod_kinds,              only: rk,ik,rdouble,TEC
    use mod_constants,          only: ONE, HALF, TWO, OUTPUT_RES
    use mod_tecio_interface,    only: init_tecio_file, init_tecio_zone, &
                                      finalize_tecio

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
    subroutine write_tecio_variables(data,filename)
        type(chidg_data_t),     intent(inout)           :: data
        character(*),           intent(in)              :: filename


        integer(ik)        :: nelem_xi, nelem_eta, nelem_zeta,  &
                              ielem_xi, ielem_eta, ielem_zeta,  &
                              npts_xi,  npts_eta,  npts_zeta,   &
                              ipt_xi,   ipt_eta,   ipt_zeta,    &
                              xilim,    etalim,    zetalim,     &
                              ielem, ielem_global, ielem_offset,&
                              npts_element, nsub_per_element,   &
                              nsub_elements, npts, icorrd, ierr,&
                              nelem, istart, ielem_start, ieq,  &
                              ivar, idom, itime, icoord, time_index

        real(rdouble)      :: val(1), r, theta, z
        real(TEC)          :: valeq(1)
        equivalence           (valeq(1), val(1))

    
        integer(4)                  :: tecstat
        integer(4),     allocatable :: connectivity(:,:)


        real(rk)                    :: xi,eta,zeta,p, sumsqr, d_normalization, grad1_d, grad2_d, grad3_d
        integer(ik)                 :: eqn_ID
        character(100)              :: varstring
        character(:),   allocatable :: zonestring



        ! using (output_res+1) so that the skip number used in tecplot to
        ! correctly display the element surfaces is the same as the number
        ! specified in the input file
        npts = OUTPUT_RES+1



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
        ! Loop time instances
        !
        !   TODO: get 'ntime' here from q_in
        !
        do itime = 1,data%sdata%q_out%get_ntime()



            do idom = 1,data%mesh%ndomains()
                !
                ! Get number of elements
                !
                nelem = data%mesh%domain(idom)%nelem


                !
                ! Initialize new zone in the TecIO file for the current domain
                !
                zonestring = 'Domain '//data%mesh%domain(idom)%name
                call init_tecio_zone(zonestring,data%mesh%domain(idom),itime)


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
                                    if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_variables: Error in call to TECDAT142")

                                end do ! ipt_xi
                            end do ! ipt_eta
                        end do ! ipt_zeta



                    end do !ielem


                end do ! coords





                ! For each variable in equation set, compute value pointwise and save
                !do ivar = 1,data%eqnset(idom)%prop%nprimary_fields()
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
                                    !   TODO: swap time_index here with 'itime' once itime gets plugged in to 'ntime' from q_in
                                    !
                                    time_index = 1
                                    !val = real(data%mesh%domain(idom)%elems(ielem)%solution_point(data%sdata%q_in%dom(idom)%vecs(ielem),ivar,itime,xi,eta,zeta),rdouble)
                                    val = real(data%mesh%domain(idom)%elems(ielem)%solution_point(data%sdata%q_out%dom(idom)%vecs(ielem),ivar,time_index,xi,eta,zeta),rdouble)
                                    tecstat = TECDAT142(1,valeq,1)
                                    if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_variables: Error in call to TECDAT142")
                                        

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
                if (tecstat /= 0) call chidg_signal(FATAL,"write_tecio_variables: Error in call to TECNOD142")

                

            end do ! idom

        end do ! itime




        !
        ! Close the current TecIO file context
        !
        call finalize_tecio()




    end subroutine write_tecio_variables
    !***********************************************************************************





end module mod_tecio
