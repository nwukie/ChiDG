module mod_spatial
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: NFACES, DIAG
    use type_chidgData, only: chidgData_t
    use type_timer,     only: timer_t

    implicit none


contains

    subroutine update_space(data,timing,info)
        implicit none
        type(chidgData_t), intent(inout)   :: data
        real(rk),          optional        :: timing
        integer(ik),       optional        :: info

        type(timer_t)               :: timer
        integer(ik)                 :: nelem, nflux, ndonors
        integer(ik)                 :: idom, ielem, iface, iblk, idonor, iflux, i, ibc
        logical                     :: skip = .false.
        logical                     :: interior_face         = .false.
        logical                     :: chimera_face          = .false.
        logical                     :: compute_face_interior = .false.
        logical                     :: compute_face_chimera  = .false.


        !
        ! Start timer on spatial discretization update
        !
        call timer%start()


            !------------------------------------------------------------------------------------------
            !                                      Interior Scheme
            !------------------------------------------------------------------------------------------

            !
            ! Loop through given element and neighbors and compute the corresponding linearization
            !
            ! XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, DIAG
            !
            do iblk = 0,7           ! (0 = linearization of chimera blocks, 1-6 = linearization of interior neighbor blocks, 7 = linearization of Q- block

                do idom = 1,size(data%domains)
                    associate ( domain => data%domains(idom), mesh => data%domains(idom)%mesh, sdata => data%sdata, eqnset => data%domains(idom)%eqnset, prop => data%domains(idom)%eqnset%prop, chimera => data%domains(idom)%mesh%chimera )

                    nelem = data%domains(idom)%mesh%nelem

                    ! Loop through elements in the domain
                    do ielem = 1,nelem


                        !
                        ! If the block direction is DIAG, then we want to compute all valid faces with neighbors.
                        ! If the block direction is not DIAG, then we only want to compute faces in the block direction
                        ! if it has a neighbor element.
                        !
                        if (iblk /= DIAG) then
                            ! Check if there is an element to linearize against in the iblk direction. If not, cycle
                            if (domain%mesh%faces(ielem,iblk)%ineighbor == 0) then
                                skip = .true.
                            else
                                skip = .false.
                            end if
                        else
                            skip = .false.  ! Don't skip DIAG
                        end if




                        !-----------------------------------------------------------------------------------------
                        if ( .not. skip) then
                            !
                            ! Faces loop. For the current element, compute the contributions from boundary integrals
                            !
                            do iface = 1,NFACES

                                associate ( face => data%domains(idom)%mesh%faces(ielem,iface) )

                                !
                                ! Only call the following routines for interior faces -- ftype == 0
                                ! Furthermore, only call the routines if we are computing derivatives for the neighbor of
                                ! iface or for the current element(DIAG). This saves a lot of unnecessary compute_boundary calls.
                                !
                                interior_face = ( domain%mesh%faces(ielem,iface)%ftype == 0 )
                                chimera_face  = ( domain%mesh%faces(ielem,iface)%ftype == 2 )
                                compute_face_interior = ( interior_face .and. (iblk == iface .or. iblk == DIAG) )
                                compute_face_chimera  = ( chimera_face  .and. (iblk == 0     .or. iblk == DIAG) )

                                if ( compute_face_interior .or. compute_face_chimera ) then



                                    !
                                    ! Call all boundary advective flux components
                                    !
                                    if (allocated(eqnset%boundary_advective_flux)) then
                                        nflux = size(eqnset%boundary_advective_flux)
                                        do iflux = 1,nflux



                                            !
                                            ! If Chimera face, how many donor elements are there that need the linearization computed
                                            !
                                            !if ( chimera_face ) then
                                            !    ichimera_recv = face%ichimera_recv
                                            !    ndonors       = chimera%recv(ichimera_recv)
                                            !else
                                                ndonors = 1
                                            !end if


                                            !
                                            ! Test ndonors > 0
                                            !
                                            if (ndonors == 0) call signal(FATAL,'update_residual: no available donors for boundary calculation')


                                            !
                                            ! Compute boundary flux once for each donor. For interior faces ndonors == 1. For Chimera faces ndonors is potentially > 1.
                                            !
                                            do idonor = 1,ndonors
                                                call eqnset%boundary_advective_flux(iflux)%flux%compute(data%domains(:)%mesh,sdata,prop,idom,ielem,iface,iblk,idonor)
                                            end do




                                        end do
                                    end if

                                end if



                                end associate

                            end do  ! faces loop
                            



                            !
                            ! Call all volume advective flux components
                            !
                            if (allocated(eqnset%volume_advective_flux)) then
                                nflux = size(eqnset%volume_advective_flux)
                                do iflux = 1,nflux
                                    call eqnset%volume_advective_flux(iflux)%flux%compute(data%domains(:)%mesh,sdata,prop,idom,ielem,iblk)
                                end do
                            end if




                        end if
                        !-----------------------------------------------------------------------------------------

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
            !> For boundary conditions, the linearization only depends on Q-, which is the solution vector
            !! for the interior element. So, we only need to compute derivatives for the interior element (DIAG)
            iblk = 7    !> DIAG
            do idom = 1,size(data%domains)
                call data%domains(idom)%bcset%apply(data%domains(:)%mesh,data%sdata,idom,iblk,data%domains(idom)%eqnset%prop)
            end do








            call timer%stop()
            call timer%report('Spatial Discretization Time')
            if (present(timing)) then
                timing = timer%elapsed()
            end if


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
