module mod_spatial
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: NFACES, DIAG
    use type_domain,    only: domain_t
    use type_timer,     only: timer_t

    implicit none


contains

    subroutine update_space(domain,timing,info)
        implicit none
        type(domain_t), intent(inout)   :: domain
        real(rk),       optional        :: timing
        integer(ik),    optional        :: info

        type(timer_t)               :: timer
        integer(ik)                 :: iblk, ielem, iface, nelem, i, ibc, iflux, nflux
        logical                     :: skip = .false.


        !
        ! Start timer on spatial discretization update
        !
        call timer%start()

        associate ( mesh => domain%mesh, sdata => domain%sdata, eqnset => domain%eqnset, prop => domain%eqnset%prop)

            nelem = mesh%nelem
            !------------------------------------------------------------------------------------------
            !                                      Interior Scheme
            !------------------------------------------------------------------------------------------

            !> Loop through given element and neighbors and compute the corresponding linearization
            !> XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, DIAG

            do iblk = 1,7

                !> Loop through elements in the domain
                do ielem = 1,nelem


                    !> If the block direction is DIAG, then we want to compute all valid faces with neighbors.
                    !! If the block direction is not DIAG, then we only want to compute faces in the block direction
                    !! if it has a neighbor element.
                    if (iblk /= DIAG) then
                        !> Check if there is an element to linearize against in the iblk direction. If not, cycle
                        if (domain%mesh%faces(ielem,iblk)%ineighbor == 0) then
                            skip = .true.
                        else
                            skip = .false.
                        end if
                    else
                        skip = .false.  !> Don't skip DIAG
                    end if




                    !-----------------------------------------------------------------------------------------
                    if ( .not. skip) then
                        ! For the current element, compute the contributions from boundary integrals
                        do iface = 1,NFACES
                            !
                            ! Only call the following routines for interior faces -- ftype == 0
                            ! Furthermore, only call the routines if we are computing derivatives for the neighbor of
                            ! iface or for the current element(DIAG). This saves a lot of unnecessary compute_boundary calls.
                            !
                            if (domain%mesh%faces(ielem,iface)%ftype == 0 .and. &
                                (iblk == iface .or. iblk == DIAG)) then



                                !
                                ! Call all boundary advective flux components
                                !
                                if (allocated(eqnset%boundary_advective_flux)) then
                                    nflux = size(eqnset%boundary_advective_flux)
                                    do iflux = 1,nflux
                                        call eqnset%boundary_advective_flux(iflux)%flux%compute(mesh,sdata,ielem,iface,iblk,prop)
                                    end do
                                end if



!                                !
!                                ! Call all boundary diffusive flux components
!                                !
!                                if (allocated(eqnset%boundary_diffusive_flux)) then
!                                    nflux = size(eqnset%boundary_diffusive_flux)
!                                    do iflux = 1,nflux
!                                        call eqnset%boundary_diffusive_flux(iflux)%flux%compute(mesh,sdata,ielem,iface,iblk,prop)
!                                    end do
!                                end if



                            end if
                        end do 



                        !
                        ! Call all volume advective flux components
                        !
                        if (allocated(eqnset%volume_advective_flux)) then
                            nflux = size(eqnset%volume_advective_flux)
                            do iflux = 1,nflux
                                call eqnset%volume_advective_flux(iflux)%flux%compute(mesh,sdata,ielem,iblk,prop)
                            end do
                        end if


!                        !
!                        ! Call all volume diffusive flux components
!                        !
!                        if (allocated(eqnset%volume_diffusive_flux)) then
!                            nflux = size(eqnset%volume_diffusive_flux)
!                            do iflux = 1,nflux
!                                call eqnset%volume_diffusive_flux(iflux)%flux%compute(mesh,sdata,ielem,iblk,prop)
!                            end do
!                        end if


                    end if
                    !-----------------------------------------------------------------------------------------

                end do !elem

            end do !block






            !------------------------------------------------------------------------------------------
            !                                      Boundary Scheme
            !------------------------------------------------------------------------------------------
            !> For boundary conditions, the linearization only depends on Q-, which is the solution vector
            !! for the interior element. So, we only need to compute derivatives for the interior element (DIAG)
            iblk = 7    !> DIAG
            call domain%bcset%apply(mesh,sdata,iblk,prop)






            call timer%stop()
            call timer%report('Spatial Discretization Time')
            if (present(timing)) then
                timing = timer%elapsed()
            end if

        end associate

    end subroutine













    !> Visit each element and compute generic function
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !---------------------------------------------------------------------------------------------
    !subroutine loop_space(domain,fcn)
    !    type(domain_t),     intent(in)  :: domain
    !    type(function_t),   intent(in)  :: fcn
    !
    !
    !
    !end subroutine






























end module mod_spatial
