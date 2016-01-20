module mod_spatial
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: NFACES, DIAG, CHIMERA, INTERIOR, XI_MAX
    use type_chidg_data,    only: chidg_data_t
    use type_timer,         only: timer_t

    use type_face_indices,  only: face_indices_t
    use type_flux_indices,  only: flux_indices_t

    use mod_DNAD_tools

    implicit none


contains

    subroutine update_space(data,timing,info)
        implicit none
        type(chidg_data_t), intent(inout)   :: data
        real(rk),           optional        :: timing
        integer(ik),        optional        :: info

        type(timer_t)               :: timer
        integer(ik)                 :: nelem, nflux, ndonors
        integer(ik)                 :: idom, ielem, iface, iblk, idonor, iflux, i, ibc, ChiID
        logical                     :: interior_face         = .false.
        logical                     :: chimera_face          = .false.
        logical                     :: compute_face          = .false.

        logical                     :: flux_needs_computed   = .false.
        logical                     :: flux_needs_linearized = .false.


        type(face_indices_t)        :: face_indices
        type(flux_indices_t)        :: flux_indices


        integer(ik) :: irow, ientry
        real(rk)    :: res


        !
        ! Start timer on spatial discretization update
        !
        call timer%start()

        

            !------------------------------------------------------------------------------------------
            !                                      Interior Scheme
            !------------------------------------------------------------------------------------------

            !
            ! Set flux linearization logicals to false
            !
            data%sdata%flux_computed    = .false.
            data%sdata%flux_linearized  = .false.








            !
            ! Loop through given element and neighbors and compute the corresponding linearization
            !
            ! XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, DIAG
            !
            do iblk = 1,7           ! (1-6 = linearization of neighbor blocks, 7 = linearization of Q- block)



                !
                ! Loop through local domains
                !
                do idom = 1,data%ndomains()
                    associate ( mesh => data%mesh(idom), sdata => data%sdata, eqnset => data%eqnset(idom)%item, prop => data%eqnset(idom)%item%prop)

                    nelem = mesh%nelem


                    !
                    ! Loop through elements in the current domain
                    !
                    do ielem = 1,nelem


                            !
                            ! Faces loop. For the current element, compute the contributions from boundary integrals
                            !
                            do iface = 1,NFACES

                                !
                                ! Define face indices
                                !
                                face_indices%idomain  = idom
                                face_indices%ielement = ielem
                                face_indices%iface    = iface



                                associate ( face => mesh%faces(ielem,iface) )

                                !
                                ! Only call the following routines for interior faces -- ftype == 0
                                ! Furthermore, only call the routines if we are computing derivatives for the neighbor of
                                ! iface or for the current element(DIAG). This saves a lot of unnecessary compute_boundary calls.
                                !
                                interior_face = ( mesh%faces(ielem,iface)%ftype == INTERIOR )
                                chimera_face  = ( mesh%faces(ielem,iface)%ftype == CHIMERA )

                                compute_face = (interior_face .or. chimera_face) .and. ( (iblk == iface) .or. (iblk == DIAG) )



                                if ( compute_face ) then

                                    !
                                    ! If Chimera face, how many donor elements are there that need the linearization computed
                                    !
                                    if ( chimera_face ) then
                                        if ( iblk /= DIAG) then ! only need to compute multiple times when we need the linearization of the chimera neighbors
                                            ChiID  = mesh%faces(ielem,iface)%ChiID
                                            ndonors = mesh%chimera%recv%data(ChiID)%ndonors
                                        else                    ! If we are linearizing the interior receiver element, only need to compute one.
                                            ndonors = 1
                                        end if
                                    else
                                        ndonors = 1
                                    end if


                                    !
                                    ! Test ndonors > 0
                                    !
                                    if (ndonors == 0) call chidg_signal(FATAL,'update_residual: no available donors for boundary calculation')




                                    !
                                    ! Call all boundary advective flux components
                                    !
                                    if (allocated(eqnset%boundary_advective_flux)) then

                                        nflux = size(eqnset%boundary_advective_flux)


                                        do iflux = 1,nflux



                                            if ( data%sdata%flux_computed(idom,ielem,iface,iflux,1) .eqv. .false. ) then
                                                flux_needs_computed = .true.
                                            else
                                                flux_needs_computed = .false.
                                            end if
                                            if ( data%sdata%flux_linearized(idom,ielem,iface,iflux,1,iblk) .eqv. .false. ) then
                                                flux_needs_linearized = .true.
                                            else
                                                flux_needs_linearized = .false.
                                            end if



                                            if ( flux_needs_computed .or. flux_needs_linearized ) then

                                                !
                                                ! Compute boundary flux once for each donor. For interior faces ndonors == 1. For Chimera faces ndonors is potentially > 1.
                                                !
                                                do idonor = 1,ndonors


                                                    face_indices%seed   = compute_seed(data%mesh,idom,ielem,iface,idonor,iblk)

                                                    flux_indices%iflux  = iflux
                                                    flux_indices%idonor = idonor
                                                    flux_indices%iblk   = iblk

                                                    call  eqnset%boundary_advective_flux(iflux)%flux%compute(data%mesh,data%sdata,prop,face_indices,flux_indices)
                                                end do

                                            end if



                                        end do ! iflux

                                    end if ! boundary_advective_flux loop




                                end if ! compute_face




                                end associate

                            end do  ! faces loop
                            



                            !
                            ! Call all volume advective flux components
                            !
                            ! only need to compute volume_advective_flux for linearization of interior block
                            if (iblk == DIAG) then
                                if (allocated(eqnset%volume_advective_flux)) then
                                    nflux = size(eqnset%volume_advective_flux)
                                    do iflux = 1,nflux

                                        call eqnset%volume_advective_flux(iflux)%flux%compute(data%mesh,data%sdata,prop,idom,ielem,iblk)

                                    end do
                                end if
                            end if




                    end do  ! ielem

                    end associate
                end do  ! idom



            end do  ! iblk







            !------------------------------------------------------------------------------------------
            !                                      Boundary Scheme
            !------------------------------------------------------------------------------------------
            !
            ! Boundary conditions
            !
            ! For boundary conditions, the linearization only depends on Q-, which is the solution vector
            ! for the interior element. So, we only need to compute derivatives for the interior element (DIAG)

            iblk = 7    !> DIAG
            do idom = 1,data%ndomains()
                call data%bcset(idom)%apply(data%mesh,data%sdata,data%eqnset(idom)%item%prop,idom,iblk)
            end do






            call timer%stop()
            call timer%report('Spatial Discretization Time')
            if (present(timing)) then
                timing = timer%elapsed()
            end if


    end subroutine update_space
    !******************************************************************************************************************




























end module mod_spatial
