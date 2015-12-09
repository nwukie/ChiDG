module mod_tecio
    use mod_kinds,              only: rk,ik,TEC
    use mod_constants,          only: ONE, HALF, TWO
    use mod_grid_operators,     only: mesh_point, solution_point
    use mod_tecio_interface,    only: init_tecio_file, init_tecio_zone, finalize_tecio

    use type_element,           only: element_t
    use type_blockvector,       only: blockvector_t
    use type_solverdata,        only: solverdata_t
    use type_chidg_data,        only: chidg_data_t

    use mod_io,                 only: OUTPUT_RES
    implicit none

#include "tecio.f90"

contains



    subroutine write_tecio_variables(data,filename,timeindex)
        type(chidg_data_t),     intent(inout)           :: data
        character(*),           intent(in)              :: filename
        integer(ik),            intent(in)              :: timeindex


        integer(ik)        :: nelem_xi, nelem_eta, nelem_zeta
        integer(ik)        :: ielem_xi, ielem_eta, ielem_zeta
        integer(ik)        :: npts_xi,  npts_eta,  npts_zeta
        integer(ik)        :: ipt_xi,   ipt_eta,   ipt_zeta
        integer(ik)        :: xilim,    etalim,    zetalim
        integer(ik)        :: npts, icoord, ielem
        integer(4)         :: tecstat

        real(rk)           :: val(1)
        real(TEC)          :: valeq(1)
        equivalence           (valeq(1), val(1))
        real(rk)           :: xi,eta,zeta
        character(100)     :: varstring
        integer(ik)        :: ieq, ivar, idom

        !type(element_t),      pointer :: elem(:,:,:)
        !type(blockvector_t),  pointer :: q


        ! using (output_res+1) so that the skip number used in tecplot to
        ! correctly display the element surfaces is the same as the number
        ! specified in the input file
        npts = OUTPUT_RES+1



        !
        ! Assemble variables string
        !
        varstring = "X Y Z"     ! Initialize variables string with mesh coordinates
        ieq = 1
        !& DEBUG - DOMAINS - Assumes same equation set in all domains
        do while (ieq <= data%eqnset(1)%item%neqns)
            !varstring = trim(varstring)//" "//trim(domain%eqnset%eqns(ieq)%name)
            varstring = trim(varstring)//" "//trim(data%eqnset(1)%item%prop%eqns(ieq)%name)
            ieq = ieq + 1
        end do


        !
        ! Open and initialize TecIO file
        !
        call init_tecio_file('solnfile',trim(varstring),filename,0)



        do idom = 1,data%ndomains()
        !associate (mesh => data%mesh(idom), eqnset => data%eqnset(idom)%item, sdata => data%sdata)

            !
            ! Store element indices for current block
            !
            nelem_xi   = data%mesh(idom)%nelem_xi
            nelem_eta  = data%mesh(idom)%nelem_eta
            nelem_zeta = data%mesh(idom)%nelem_zeta


            !
            ! Remap elements array to block matrix
            !
            !elem => mesh%elems_m
            !q    => sdata%q%dom(idom)


            !
            ! Initialize new zone in the TecIO file for the current domain
            !
            call init_tecio_zone('solnzone',data%mesh(idom),1,timeindex)



            xilim   = npts
            etalim  = npts
            zetalim = npts

            ! For each coordinate, compute it's value pointwise and save
            do icoord = 1,3

               ! Loop through elements and get structured points
                do ielem_zeta = 1,nelem_zeta
                    do ipt_zeta = 1,zetalim
                        zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                        do ielem_eta = 1,nelem_eta
                            do ipt_eta = 1,etalim
                                eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                do ielem_xi = 1,nelem_xi
                                    do ipt_xi = 1,xilim
                                        xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                        ! Get coordinate value at point
                                        !val = mesh_point(elem(ielem_xi,ielem_eta,ielem_zeta),icoord,xi,eta,zeta)
                                        !val = mesh_point(data%mesh(idom)%elems_m(ielem_xi,ielem_eta,ielem_zeta),icoord,xi,eta,zeta)
                                        ielem = ielem_xi + (nelem_xi)*(ielem_eta-1) + (nelem_xi * nelem_eta)*(ielem_zeta-1)
                                        val = mesh_point(data%mesh(idom)%elems(ielem),icoord,xi,eta,zeta)
                                        tecstat = TECDAT142(1,valeq,1)

                                    end do
                                end do

                            end do
                        end do

                    end do
                end do

            end do  ! coords





            ! For each variable in equation set, compute value pointwise and save
            do ivar = 1,data%eqnset(idom)%item%neqns

                do ielem_zeta = 1,nelem_zeta
                    do ipt_zeta = 1,zetalim
                        zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                        do ielem_eta = 1,nelem_eta
                            do ipt_eta = 1,etalim
                                eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                do ielem_xi = 1,nelem_xi
                                    do ipt_xi = 1,xilim
                                        xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                        ! Get solution value at point
                                        !val = solution_point(q%lvecs_m(ielem_xi,ielem_eta,ielem_zeta),ivar,xi,eta,zeta)
                                        ielem = ielem_xi + (nelem_xi)*(ielem_eta-1) + (nelem_xi * nelem_eta)*(ielem_zeta-1)
                                        !val = solution_point(q%lvecs(ielem),ivar,xi,eta,zeta)
                                        val = solution_point(data%sdata%q%dom(idom)%lvecs(ielem),ivar,xi,eta,zeta)

                                        tecstat = TECDAT142(1,valeq,1)
                                    
                                    
                                    end do
                                end do

                            end do
                        end do

                    end do
                end do

            end do ! ivar





        !end associate
        end do ! idom




        !
        ! Close the current TecIO file context
        !
        call finalize_tecio()




    end subroutine















































end module mod_tecio
