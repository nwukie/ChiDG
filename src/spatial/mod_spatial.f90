module mod_spatial
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: NFACES, DIAG, CHIMERA, INTERIOR, XI_MAX, &
                                  BOUNDARY_ADVECTIVE_FLUX
    use type_chidg_data,    only: chidg_data_t
    use type_timer,         only: timer_t

    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t

    use mod_DNAD_tools

    implicit none


contains


    !> Spatial loop through domains, elements, and faces. Functions get called for each element/face.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------------------------------------
    subroutine update_space(data,timing,info)
        implicit none
        type(chidg_data_t), intent(inout)   :: data
        real(rk),           optional        :: timing
        integer(ik),        optional        :: info

        type(timer_t)               :: timer
        integer(ik)                 :: nelem, nfcn, ndonors
        integer(ik)                 :: idom, ielem, iface, iblk, idonor, ifcn, i, ibc, ChiID
        logical                     :: interior_face         = .false.
        logical                     :: chimera_face          = .false.
        logical                     :: compute_face          = .false.

        logical                     :: compute_function      = .false.
        logical                     :: linearize_function    = .false.


        type(face_info_t)           :: face_info
        type(function_info_t)       :: function_info


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
            ! Clear function_status data. This tracks if a function has already been called. So, in this way
            ! we can compute a function on a face and apply it to both elements. The function is just registered
            ! as computed for both. So we need to reset all of that data here. This is only tracked for the interior scheme.
            ! Boundary condition evaluations and Chimera faces are not tracked.
            !
            call data%sdata%function_status%clear()








            !
            ! Loop through given element and neighbors and compute the corresponding linearization
            !
            ! XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, DIAG
            !
            do iblk = 1,7           ! 1-6 = linearization of neighbor blocks, 7 = linearization of Q- block(self)



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
                                face_info%idomain  = idom
                                face_info%ielement = ielem
                                face_info%iface    = iface



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
                                    ! TODO: Relocate this to a subroutine else-where
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
                                        nfcn = size(eqnset%boundary_advective_flux)
                                        do ifcn = 1,nfcn

                                            function_info%type   = BOUNDARY_ADVECTIVE_FLUX
                                            function_info%ifcn   = ifcn
                                            function_info%iblk   = iblk

                                            compute_function     = data%sdata%function_status%compute_function(   face_info, function_info )
                                            linearize_function   = data%sdata%function_status%linearize_function( face_info, function_info )




                                            if ( compute_function .or. linearize_function ) then
                                                !
                                                ! Compute boundary flux once for each donor. For interior faces ndonors == 1. For Chimera faces ndonors is potentially > 1.
                                                !
                                                do idonor = 1,ndonors
                                                    face_info%seed       = compute_seed(data%mesh,idom,ielem,iface,idonor,iblk)
                                                    function_info%idonor = idonor

                                                    call eqnset%boundary_advective_flux(ifcn)%flux%compute(data%mesh,data%sdata,prop,face_info,function_info)
                                                end do
                                            end if


                                        end do ! ifcn

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
                                    nfcn = size(eqnset%volume_advective_flux)
                                    do ifcn = 1,nfcn

                                        call eqnset%volume_advective_flux(ifcn)%flux%compute(data%mesh,data%sdata,prop,idom,ielem,iblk)

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
            ! Apply boundary conditions for each domain.
            !
            do idom = 1,data%ndomains()

                call data%bcset(idom)%apply(data%mesh,data%sdata,data%eqnset(idom)%item%prop,idom)

            end do ! idom






            call timer%stop()
            call timer%report('Spatial Discretization Time')
            if (present(timing)) then
                timing = timer%elapsed()
            end if


    end subroutine update_space
    !******************************************************************************************************************




























end module mod_spatial
