module mod_check_jacobian
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: NFACES, DIAG, ZERO

    use type_domain,        only: domain_t
    use type_solverdata,    only: solverdata_t
    use type_blockvector
    use type_densematrix,   only: densematrix_t

    implicit none


contains





    ! Compute both the Finite Difference representation and DNAD representation of
    ! the volume advective flux jacobian.
    !
    !   @author Nathan A. Wukie
    !
    !   @param[inout]   domain      Domain structure containing grid, eqnset, etc
    !   @param[in]      ielem       Element index for computing the jacobian
    !   @param[in]      iblk        Block index of the linearization
    !   @param[inout]   blk_dnad    Densematrix for storing the DNAD jacobian
    !   @param[inout]   blk_fd      Densematrix for storing the Finite Difference jacobian
    !
    !----------------------------------------------------------------------------------
    subroutine check_jacobian_volume_advective_flux(domains,sdata,ielem,iblk,blk_dnad,blk_fd)
        type(domain_t),      intent(inout)       :: domains(:)
        type(solverdata_t),  intent(inout)       :: sdata
        integer(ik),         intent(in)          :: ielem, iblk
        type(densematrix_t), intent(inout)       :: blk_dnad, blk_fd

        type(blockvector_t), allocatable    :: rhs_r, vec_fd
        real(rk)                            :: qhold, eps
        integer(ik)                         :: nelem, i, iterm, icol, nterms, ivar, iflux, nflux, idom
        integer(ik)                         :: ielem_p              ! element in which the solution is being perturbed for the FD calculation.


        ! ASSUME ONLY ONE DOMAIN IS PRESENT
        idom = 1

        associate ( domain => domains(1), mesh => domains(:)%mesh, q => sdata%q%dom(1), rhs => sdata%rhs%dom(1), lhs => sdata%lhs%dom(1), prop => domains(1)%eqnset%prop )


            nelem = mesh(1)%nelem
            nterms = mesh(1)%nterms_s
            !------------------------------------------------------------------------------------------
            !                                      Interior Scheme
            !------------------------------------------------------------------------------------------
            !
            ! Clear all working data
            !
            call sdata%lhs%clear()
            call sdata%rhs%clear()
            rhs_r  = rhs      ! Alloate supporting storage
            vec_fd = rhs

            !
            ! For the current element, compute the contributions from volume integrals
            !
            !call domain%eqnset%compute_volume_flux(mesh,sdata,ielem,iblk)
            if (allocated(domain%eqnset%volume_advective_flux)) then
                nflux = size(domain%eqnset%volume_advective_flux)
                do iflux = 1,nflux
                    call domain%eqnset%volume_advective_flux(iflux)%flux%compute(mesh,sdata,prop,idom,ielem,iblk)
                end do

            else
                print*, domain%eqnset%name
                call signal(WARN,"No volume advective flux was found")
            end if


            !
            ! Store linearization from DNAD. This is what we want to check against the FD calculation.
            !
            !blk_dnad = sdata%lhs%dom(idom)%lblks(ielem,iblk)
            blk_dnad = lhs%lblks(ielem,iblk)
            blk_fd = blk_dnad                       ! Sourced allocation
            blk_fd%mat = ZERO                       ! Zero storage


            !
            ! Store temporary rhs
            !
            rhs_r%lvecs(ielem) = rhs%lvecs(ielem)


            ! Reset sdata storage
            call lhs%clear()
            call rhs%clear()


            !
            ! Select in which element, the solution is being perturbed
            !
            if (iblk == DIAG) then
                ielem_p = ielem
            else
                ielem_p = domain%mesh%faces(ielem,iblk)%ineighbor
            end if




            eps   = 1.e-8_rk
            ! Loop through terms, perturb term, compute rhs, compute finite difference jacobian, return term.
            do ivar = 1,domain%eqnset%neqns
                do iterm = 1,domain%mesh%nterms_s

                    !
                    ! Perturb the iterm-th term in the solution expansion for variable ivar in element ielem.
                    !
                    qhold = q%lvecs(ielem_p)%getterm(ivar,iterm)
                    call q%lvecs(ielem_p)%setterm(ivar,iterm,qhold + eps)


                    !
                    ! For the current element, compute the contributions from volume integrals
                    !
                    !call domain%eqnset%compute_volume_flux(mesh,sdata,ielem,iblk)
                    if (allocated(domain%eqnset%volume_advective_flux)) then
                        nflux = size(domain%eqnset%volume_advective_flux)
                        do iflux = 1,nflux
                            call domain%eqnset%volume_advective_flux(iflux)%flux%compute(mesh,sdata,prop,idom,ielem,iblk)
                        end do

                    else
                        print*, domain%eqnset%name
                        call signal(WARN,"No volume advective flux was found")
                    end if




                    !
                    ! Return perturbed value to normal state
                    !
                    call q%lvecs(ielem_p)%setterm(ivar,iterm,qhold)


                    !
                    ! Compute finite difference jacobian
                    !
                    vec_fd = (rhs - rhs_r)/eps


                    !
                    ! Store to column of blk_fd
                    !
                    icol = (ivar-1)*nterms + iterm                  !> Compute appropriate column for storing linearization
                    blk_fd%mat(:,icol) = vec_fd%lvecs(ielem)%vec    !> Store finite difference linearization of the residual


                    ! Reset sdata storage
                    call sdata%lhs%clear()
                    call sdata%rhs%clear()

                end do
            end do



        end associate

    end subroutine





    ! compute both the finite difference representation and dnad representation of
    ! the boundary advective flux jacobian.
    !
    !   @author nathan a. wukie
    !
    !   @param[inout]   domain      domain structure containing grid, eqnset, etc
    !   @param[in]      ielem       element index for computing the jacobian
    !   @param[in]      iblk        block index of the linearization
    !   @param[inout]   blk_dnad    densematrix for storing the dnad jacobian
    !   @param[inout]   blk_fd      densematrix for storing the finite difference jacobian
    !
    !----------------------------------------------------------------------------------
    subroutine check_jacobian_boundary_advective_flux(domains,sdata,ielem,iblk,blk_dnad,blk_fd)
        type(domain_t),         intent(inout)   :: domains(:)
        type(solverdata_t),     intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iblk
        type(densematrix_t),    intent(inout)   :: blk_dnad, blk_fd

        type(blockvector_t)  :: rhs_r, vec_fd
        real(rk)    :: qhold, eps
        integer(ik) :: nelem, i, iterm, iface, ivar, icol, nterms, iflux, nflux, idom, idonor
        integer(ik) :: ielem_p              !> ielem_p is the element in which the solution is being perturbed for the finite difference calculation.


        ! ASSUMING THERE IS ONLY ONE DOMAIN
        idom   = 1
        idonor = 1

        !
        ! Select in which element, the solution is being perturbed
        ! Select which face the boundary integral is being computed on based on the direction of iblk.
        !
        if (iblk == DIAG) then
            ielem_p = ielem
            iface   = 1
        else
            ielem_p = domains(1)%mesh%faces(ielem,iblk)%ineighbor
            iface   = iblk
        end if


        !associate ( mesh => domain%mesh, sdata => domain%sdata, prop => domain%eqnset%prop) 
        associate ( domain => domains(1), mesh => domains(:)%mesh, q => sdata%q%dom(1), rhs => sdata%rhs%dom(1), lhs => sdata%lhs%dom(1), prop => domains(1)%eqnset%prop )

            nelem = mesh(1)%nelem
            nterms = mesh(1)%nterms_s
            !------------------------------------------------------------------------------------------
            !                                      Interior Scheme
            !------------------------------------------------------------------------------------------

            !
            ! Zero data storage
            !
            call lhs%clear()      
            call lhs%clear()
            rhs_r  = rhs          ! Allocate supporing storage containers
            vec_fd = rhs


            !
            ! For the current element, compute the contributions from boundary integrals
            !
            !call domain%eqnset%compute_boundary_average_flux(mesh,sdata,ielem,iface,iblk)
            if (allocated(domain%eqnset%boundary_advective_flux)) then
                nflux = size(domain%eqnset%boundary_advective_flux)
                do iflux = 1,nflux
                    call domain%eqnset%boundary_advective_flux(iflux)%flux%compute(mesh,sdata,prop,idom,ielem,iface,iblk,idonor)
                end do

            else
                print*, domain%eqnset%name
                call signal(WARN,"No boundary advective flux was found")
            end if


            !
            ! Store linearization computed from DNAD
            !
            blk_dnad = lhs%lblks(ielem,iblk)
            blk_fd = blk_dnad   ! Sourced allocation
            blk_fd%mat = ZERO   ! Zero storage



            ! Reset sdata storage
            call lhs%clear()
            call rhs%clear()



            !
            ! For the current element, compute the contributions from boundary integrals
            !
            ! Need to use DIAG to get rhs for finite difference calculation. 
            ! This is because RHS is only stored for DIAG in the integrate procedure.
            !call domain%eqnset%compute_boundary_average_flux(mesh,sdata,ielem,iface,DIAG)       
            if (allocated(domain%eqnset%boundary_advective_flux)) then
                nflux = size(domain%eqnset%boundary_advective_flux)
                do iflux = 1,nflux
                    call domain%eqnset%boundary_advective_flux(iflux)%flux%compute(mesh,sdata,prop,idom,ielem,iface,DIAG,idonor)
                end do

            else
                print*, domain%eqnset%name
                call signal(WARN,"No boundary advective flux was found")
            end if









            ! Store temporary rhs
            rhs_r%lvecs(ielem) = rhs%lvecs(ielem)


            ! Reset sdata storage
            call lhs%clear()
            call rhs%clear()




            eps   = 1.e-8_rk
            ! Loop through terms, perturb term, compute rhs, compute finite difference jacobian, return term.
            do ivar = 1,domain%eqnset%neqns
                do iterm = 1,domain%mesh%nterms_s

                    !
                    ! Perturb the iterm-th term in the solution expansion for variable ivar in element ielem.
                    !
                    qhold = q%lvecs(ielem_p)%getterm(ivar,iterm)
                    call q%lvecs(ielem_p)%setterm(ivar,iterm,qhold + eps)


                    !
                    ! For the current element, compute the contributions from volume integrals
                    !
                    ! Need to use DIAG to get rhs for finite difference calculation. 
                    ! This is because RHS is only stored for DIAG in the integrate procedure.
                    !call domain%eqnset%compute_boundary_average_flux(mesh,sdata,ielem,iface,DIAG)    
                    if (allocated(domain%eqnset%boundary_advective_flux)) then
                        nflux = size(domain%eqnset%boundary_advective_flux)
                        do iflux = 1,nflux
                            call domain%eqnset%boundary_advective_flux(iflux)%flux%compute(mesh,sdata,prop,idom,ielem,iface,DIAG,idonor)
                        end do

                    else
                        print*, domain%eqnset%name
                        call signal(WARN,"No boundary advective flux was found")
                    end if






                    !
                    ! Return perturbed value to normal state
                    !
                    call q%lvecs(ielem_p)%setterm(ivar,iterm,qhold)


                    !
                    ! Compute finite difference jacobian
                    !
                    vec_fd = (rhs - rhs_r)/eps


                    !
                    ! Store to column of blk_fd
                    !
                    icol = (ivar-1)*nterms + iterm                  ! Compute appropriate column for storing linearization
                    blk_fd%mat(:,icol) = vec_fd%lvecs(ielem)%vec


                    !
                    ! Reset sdata storage
                    !
                    call lhs%clear()
                    call rhs%clear()

                end do
            end do


        end associate

    end subroutine








end module mod_check_jacobian
