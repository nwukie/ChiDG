module mod_tecio
    use mod_kinds,              only: rk,ik,TEC
    use mod_constants,          only: IO_RES, ONE, HALF, TWO
    use type_domain,            only: domain_t
    use type_element,           only: element_t
    use type_expansion,         only: expansion_t
    use mod_grid_operators,     only: mesh_point, solution_point
    use mod_tecio_interface,    only: init_tecio_file, init_tecio_zone, finalize_tecio

    implicit none

#include "tecio.f90"

contains



    subroutine write_tecio_variable(domain,filename,timeindex)
        type(domain_t), intent(in), target  :: domain
        character(*),   intent(in)          :: filename
        integer(ik),    intent(in)          :: timeindex


        integer(ik)        :: nelem_xi, nelem_eta, nelem_zeta
        integer(ik)        :: ielem_xi, ielem_eta, ielem_zeta
        integer(ik)        :: npts_xi,  npts_eta,  npts_zeta
        integer(ik)        :: ipt_xi,   ipt_eta,   ipt_zeta
        integer(ik)        :: xilim,    etalim,    zetalim
        integer(ik)        :: npts, icoord
        integer(4)         :: tecstat

        real(rk)           :: val(1)
        real(TEC)          :: valeq(1)
        equivalence           (valeq(1), val(1))
        real(rk)           :: xi,eta,zeta
        character(100)     :: varstring
        integer(ik)        :: ieq, ivar

        type(element_t),    pointer :: elems(:,:,:)
        type(expansion_t),  pointer :: q(:,:,:)


        ! Store element indices for current block
        nelem_xi   = domain%mesh%nelem_xi
        nelem_eta  = domain%mesh%nelem_eta
        nelem_zeta = domain%mesh%nelem_zeta

        ! Remap elements array to block matrix
        elems(1:nelem_xi,1:nelem_eta,1:nelem_zeta) => domain%mesh%elems
        q(1:nelem_xi,1:nelem_eta,1:nelem_zeta)     => domain%q

        ! using (output_res+1) so that the skip number used in tecplot to
        ! correctly display the element surfaces is the same as the number
        ! specified in the input file
        npts = IO_RES+1

        ! Initialize variables string with mesh coordinates
        varstring = "X Y Z"

        ! Assemble variables string
        ieq = 1
        do while (ieq <= domain%eqnset%neqns)
            varstring = trim(varstring)//" "//trim(domain%eqnset%eqns(ieq)%name)
            ieq = ieq + 1
        end do



        ! Initialize TECIO binary for solution file
        call init_tecio_file('solnfile',trim(varstring),filename,0)
        call init_tecio_zone('solnzone',domain%mesh,0,timeindex)


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
                                    val = mesh_point(elems(ielem_xi,ielem_eta,ielem_zeta),icoord,xi,eta,zeta)
                                    tecstat = TECDAT142(1,valeq,1)

                                end do
                            end do

                        end do
                    end do

                end do
            end do

        end do  ! coords


        ! For each variable in equation set, compute value pointwise and save
        do ivar = 1,domain%eqnset%neqns

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
                                    val = solution_point(q(ielem_xi,ielem_eta,ielem_zeta),ivar,xi,eta,zeta)
                                    tecstat = TECDAT142(1,valeq,1)

                                end do
                            end do

                        end do
                    end do

                end do
            end do

        end do ! ivar


        call finalize_tecio()


    end subroutine



















!    subroutine write_tecio_mesh(domain,writetype,filename,timeindex)
!        type(domain_t), intent(in) :: domain
!        character(*),   intent(in) :: writetype
!        character(*),   intent(in) :: filename
!        integer(ik),    intent(in) :: timeindex
!
!
!        integer(ik)        :: nelem_xi, nelem_eta, nelem_zeta
!        integer(ik)        :: ielem_xi, ielem_eta, ielem_zeta
!        integer(ik)        :: npts_xi,  npts_eta,  npts_zeta
!        integer(ik)        :: ipt_xi,   ipt_eta,   ipt_zeta
!        integer(ik)        :: xilim,    etalim,    zetalim
!        integer(ik)        :: npts
!        integer(4)         :: tecstat
!
!        real(rk)           :: xval(1),yval(1),zval(1),val(1)
!        real(TEC)          :: xeq(1),yeq(1),zeq(1),valeq(1)
!        equivalence           (xeq(1),   xval(1))
!        equivalence           (yeq(1),   yval(1))
!        equivalence           (zeq(1),   zval(1))
!        equivalence           (valeq(1), val(1))
!        real(rk)           :: xi,eta,zeta
!
!        ! Store element indices for current block
!        nelem_xi   = domain%mesh%nelem_xi
!        nelem_eta  = domain%mesh%nelem_eta
!        nelem_zeta = domain%mesh%nelem_zeta
!
!        ! using (output_res+1) so that the skip number used in tecplot to
!        ! correctly display the element surfaces is the same as the number
!        ! specified in the input file
!        npts = IO_RES+1
!
!
!        !=====================================================
!        !
!        !   Write mesh routine
!        !
!        !=====================================================
!        ! Initialize TECIO binary for mesh file
!        call init_tecio_file('meshfile','X Y Z',filename,1)
!        call init_tecio_zone('meshzone',domain%mesh,0,1)
!
!
!
!        ! Write x-coordinate
!        do ielem_zeta = 1,nelem_zeta
!            if (ielem_zeta == nelem_zeta) then
!                zetalim = npts
!            else
!                zetalim = npts-1
!            end if
!            do ipt_zeta = 1,zetalim
!                zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
!
!                do ielem_eta = 1,nelem_eta
!                    if (ielem_eta == nelem_eta) then
!                        etalim = npts
!                    else
!                        etalim = npts-1
!                    end if
!                    do ipt_eta = 1,etalim
!                        eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
!
!                        do ielem_xi = 1,nelem_xi
!                            if (ielem_xi == nelem_xi) then
!                                xilim = npts
!                            else
!                                xilim = npts-1
!                            end if
!
!                            do ipt_xi = 1,xilim
!                                xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
!
!                                ! Get x-mesh point
!                                xval = self%mesh%elements(ielem_xi,ielem_eta,ielem_zeta)%x(xi,eta,zeta)
!                                tecstat = TECDAT142(1,xeq,1)
!
!                            end do
!                        end do
!
!                    end do
!                end do
!
!            end do
!        end do
!
!        ! Write y-coordinate
!        do ielem_zeta = 1,nelem_zeta
!            if (ielem_zeta == nelem_zeta) then
!                zetalim = npts
!            else
!                zetalim = npts-1
!            end if
!            do ipt_zeta = 1,zetalim
!                zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
!
!                do ielem_eta = 1,nelem_eta
!                    if (ielem_eta == nelem_eta) then
!                        etalim = npts
!                    else
!                        etalim = npts-1
!                    end if
!
!                    do ipt_eta = 1,etalim
!                        eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
!
!                        do ielem_xi = 1,nelem_xi
!                            if (ielem_xi == nelem_xi) then
!                                xilim = npts
!                            else
!                                xilim = npts-1
!                            end if
!
!                            do ipt_xi = 1,xilim
!                                xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
!
!                                ! Get y-mesh point
!                                yval = self%mesh%elements(ielem_xi,ielem_eta,ielem_zeta)%y(xi,eta,zeta)
!                                tecstat = TECDAT142(1,yeq,1)
!
!                            end do
!                        end do
!
!                    end do
!                end do
!
!            end do
!        end do
!
!
!        ! Write z-coordinate
!        do ielem_zeta = 1,nelem_zeta
!            if (ielem_zeta == nelem_zeta) then
!                zetalim = npts
!            else
!                zetalim = npts-1
!            end if
!
!            do ipt_zeta = 1,zetalim
!                zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
!
!                do ielem_eta = 1,nelem_eta
!                    if (ielem_eta == nelem_eta) then
!                        etalim = npts
!                    else
!                        etalim = npts-1
!                    end if
!
!                    do ipt_eta = 1,etalim
!                        eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
!
!                        do ielem_xi = 1,nelem_xi
!                            if (ielem_xi == nelem_xi) then
!                                xilim = npts
!                            else
!                                xilim = npts-1
!                            end if
!
!                            do ipt_xi = 1,xilim
!                                xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO
!
!                                ! Get z-mesh point
!                                zval = self%mesh%elements(ielem_xi,ielem_eta,ielem_zeta)%z(xi,eta,zeta)
!                                tecstat = TECDAT142(1,zeq,1)
!
!                            end do
!                        end do
!
!                    end do
!                end do
!
!            end do
!        end do
!
!        call finalize_tecio()
!
!
!    end subroutine
!

































end module mod_tecio
