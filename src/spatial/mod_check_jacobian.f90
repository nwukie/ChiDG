module mod_check_jacobian
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: NFACES, DIAG, ZERO
    use type_domain,        only: domain_t
    use type_expansion
    use type_densematrix,   only: densematrix_t

    implicit none


contains

    subroutine check_jacobian_volume_flux(domain,ielem,ivar,iblk,blk_dnad,blk_fd)
        type(domain_t),      intent(inout)       :: domain
        integer(ik),         intent(in)          :: ielem, ivar, iblk
        type(densematrix_t), intent(inout)       :: blk_dnad, blk_fd

        type(expansion_t), allocatable  :: rhs_r(:), vec_fd(:)
        real(rk)    :: qhold, eps
        integer(ik) :: nelem, i, iterm
        integer(ik) :: ielem_p              !> ielem_p is the element in which the solution is being perturbed for the finite difference calculation.



        associate ( mesh => domain%mesh, sdata => domain%sdata)

            nelem = mesh%nelem
            !------------------------------------------------------------------------------------------
            !                                      Interior Scheme
            !------------------------------------------------------------------------------------------

            call sdata%lin%clear()          !> Ensure linearization is zero-ed
            do i=1,size(sdata%rhs)
                sdata%rhs(i)%vec = ZERO     !> Ensure RHS vector is zero-ed
            end do
            rhs_r  = sdata%rhs              !> should result in sourced allocation
            vec_fd = sdata%rhs


            !> For the current element, compute the contributions from volume integrals
            call domain%eqnset%compute_volume_flux(mesh,sdata,ielem,iblk)


            blk_dnad = sdata%lin%lblks(ielem,iblk)  !> Store linearization from DNAD. This is what we want to check against the finite difference calculation.
            blk_fd = blk_dnad                       !> Sourced allocation
            blk_fd%mat = ZERO                       !> Zero storage


            !> Store temporary rhs
            rhs_r(ielem) = sdata%rhs(ielem)


            !> Reset sdata storage
            call sdata%lin%clear()
            do i=1,size(sdata%rhs)
                sdata%rhs(i)%vec = ZERO
            end do


            !> Select in which element, the solution is being perturbed
            if (iblk == DIAG) then
                ielem_p = ielem
            else
                ielem_p = domain%mesh%faces(ielem,iblk)%ineighbor
            end if




            !> Loop through terms, perturb term, compute rhs, compute finite difference jacobian, return term.
            do iterm = 1,domain%mesh%nterms_s

                !> Perturb the iterm-th term in the solution expansion for variable ivar in element ielem.
                eps   = 1.e-8_rk
                qhold = sdata%q(ielem_p)%mat(iterm,ivar)
                sdata%q(ielem_p)%mat(iterm,ivar) = sdata%q(ielem_p)%mat(iterm,ivar) + eps


                !> For the current element, compute the contributions from volume integrals
                call domain%eqnset%compute_volume_flux(mesh,sdata,ielem,iblk)


                sdata%q(ielem_p)%mat(iterm,ivar) = qhold  !> Return perturbed value to normal state

                !> Compute finite difference jacobian
                vec_fd = (sdata%rhs - rhs_r)/eps


                !> Store to column of blk_fd
                blk_fd%mat(:,iterm) = vec_fd(ielem)%vec


                !> Reset sdata storage
                call sdata%lin%clear()
                do i=1,size(sdata%rhs)
                    sdata%rhs(i)%vec = ZERO
                end do


            end do


        end associate

    end subroutine






    subroutine check_jacobian_boundary_average_flux(domain,ielem,ivar,iblk,blk_dnad,blk_fd)
        type(domain_t),      intent(inout)       :: domain
        integer(ik),         intent(in)          :: ielem, ivar, iblk
        type(densematrix_t), intent(inout)       :: blk_dnad, blk_fd

        type(expansion_t), allocatable  :: rhs_r(:), vec_fd(:)
        real(rk)    :: qhold, eps
        integer(ik) :: nelem, i, iterm, iface
        integer(ik) :: ielem_p              !> ielem_p is the element in which the solution is being perturbed for the finite difference calculation.


        !> Select in which element, the solution is being perturbed
        !> Select which face the boundary integral is being computed on based on the direction of iblk.
        if (iblk == DIAG) then
            ielem_p = ielem
            iface   = 1
        else
            ielem_p = domain%mesh%faces(ielem,iblk)%ineighbor
            iface   = iblk
        end if


        associate ( mesh => domain%mesh, sdata => domain%sdata)

            nelem = mesh%nelem
            !------------------------------------------------------------------------------------------
            !                                      Interior Scheme
            !------------------------------------------------------------------------------------------

            call sdata%lin%clear()          !> Ensure linearization is zero-ed
            do i=1,size(sdata%rhs)
                sdata%rhs(i)%vec = ZERO     !> Ensure RHS vector is zero-ed
            end do
            rhs_r  = sdata%rhs              !> should result in sourced allocation
            vec_fd = sdata%rhs


            !> For the current element, compute the contributions from boundary integrals
            call domain%eqnset%compute_boundary_average_flux(mesh,sdata,ielem,iface,iblk)


            !> Store linearization computed from DNAD
            blk_dnad = sdata%lin%lblks(ielem,iblk)
            blk_fd = blk_dnad   !> Sourced allocation
            blk_fd%mat = ZERO   !> Zero storage



            !> Reset sdata storage
            call sdata%lin%clear()
            do i=1,size(sdata%rhs)
                sdata%rhs(i)%vec = ZERO
            end do



            !> For the current element, compute the contributions from boundary integrals
            call domain%eqnset%compute_boundary_average_flux(mesh,sdata,ielem,iface,DIAG)       !> Need to use DIAG to get rhs for finite difference calculation. This is because RHS is only stored for DIAG in the integrate procedure.

            !> Store temporary rhs
            rhs_r(ielem) = sdata%rhs(ielem)


            !> Reset sdata storage
            call sdata%lin%clear()
            do i=1,size(sdata%rhs)
                sdata%rhs(i)%vec = ZERO
            end do




            !> Loop through terms, perturb term, compute rhs, compute finite difference jacobian, return term.
            do iterm = 1,domain%mesh%nterms_s

                !> Perturb the iterm-th term in the solution expansion for variable ivar in element ielem.
                eps   = 1.e-8_rk
                qhold = sdata%q(ielem_p)%mat(iterm,ivar)
                sdata%q(ielem_p)%mat(iterm,ivar) = sdata%q(ielem_p)%mat(iterm,ivar) + eps


                !> For the current element, compute the contributions from volume integrals
                call domain%eqnset%compute_boundary_average_flux(mesh,sdata,ielem,iface,DIAG)    !> Need to use DIAG to get rhs for finite difference calculation. This is because RHS is only stored for DIAG in the integrate procedure.


                sdata%q(ielem_p)%mat(iterm,ivar) = qhold  !> Return perturbed value to normal state

                !> Compute finite difference jacobian
                vec_fd = (sdata%rhs - rhs_r)/eps


                !> Store to column of blk_fd
                blk_fd%mat(:,iterm) = vec_fd(ielem)%vec


                !> Reset sdata storage
                call sdata%lin%clear()
                do i=1,size(sdata%rhs)
                    sdata%rhs(i)%vec = ZERO
                end do


            end do


        end associate

    end subroutine




   subroutine check_jacobian_boundary_upwind_flux(domain,ielem,ivar,iblk,blk_dnad,blk_fd)
        type(domain_t),      intent(inout)       :: domain
        integer(ik),         intent(in)          :: ielem, ivar, iblk
        type(densematrix_t), intent(inout)       :: blk_dnad, blk_fd

        type(expansion_t), allocatable  :: rhs_r(:), vec_fd(:)
        real(rk)    :: qhold, eps
        integer(ik) :: nelem, i, iterm, iface
        integer(ik) :: ielem_p              !> ielem_p is the element in which the solution is being perturbed for the finite difference calculation.



        !> Select in which element, the solution is being perturbed
        !> Select which face the boundary integral is being computed on based on the direction of iblk.
        if (iblk == DIAG) then
            ielem_p = ielem
            iface   = 1
        else
            ielem_p = domain%mesh%faces(ielem,iblk)%ineighbor
            iface   = iblk
        end if



        associate ( mesh => domain%mesh, sdata => domain%sdata)

            nelem = mesh%nelem
            !------------------------------------------------------------------------------------------
            !                                      Interior Scheme
            !------------------------------------------------------------------------------------------

            call sdata%lin%clear()          !> Ensure linearization is zero-ed
            do i=1,size(sdata%rhs)
                sdata%rhs(i)%vec = ZERO     !> Ensure RHS vector is zero-ed
            end do
            rhs_r  = sdata%rhs              !> should result in sourced allocation
            vec_fd = sdata%rhs


            !> For the current element, compute the contributions from volume integrals
            call domain%eqnset%compute_boundary_upwind_flux(mesh,sdata,ielem,iface,iblk)


            !> Store linearization computed from DNAD
            blk_dnad = sdata%lin%lblks(ielem,iblk)
            blk_fd = blk_dnad   !> Sourced allocation
            blk_fd%mat = ZERO   !> Zero storage


            !> Reset sdata storage
            call sdata%lin%clear()
            do i=1,size(sdata%rhs)
                sdata%rhs(i)%vec = ZERO
            end do




            !> For the current element, compute the contributions from volume integrals
            call domain%eqnset%compute_boundary_upwind_flux(mesh,sdata,ielem,iface,DIAG)    !> Need to use DIAG to get RHS for finite difference calculation. This is because RHS is only stored for DIAG in the integrate procedure.


            !> Store temporary rhs
            rhs_r(ielem) = sdata%rhs(ielem)


            !> Reset sdata storage
            call sdata%lin%clear()
            do i=1,size(sdata%rhs)
                sdata%rhs(i)%vec = ZERO
            end do




            !> Loop through terms, perturb term, compute rhs, compute finite difference jacobian, return term.
            do iterm = 1,domain%mesh%nterms_s

                !> Perturb the iterm-th term in the solution expansion for variable ivar in element ielem.
                eps   = 1.e-8_rk
                qhold = sdata%q(ielem_p)%mat(iterm,ivar)
                sdata%q(ielem_p)%mat(iterm,ivar) = sdata%q(ielem_p)%mat(iterm,ivar) + eps


                !> For the current element, compute the contributions from volume integrals
                call domain%eqnset%compute_boundary_average_flux(mesh,sdata,ielem,iface,DIAG)    !> Need to use DIAG to get rhs for finite difference calculation. This is because RHS is only stored for DIAG in the integrate procedure.


                sdata%q(ielem_p)%mat(iterm,ivar) = qhold  !> Return perturbed value to normal state

                !> Compute finite difference jacobian
                vec_fd = (sdata%rhs - rhs_r)/eps


                !> Store to column of blk_fd
                blk_fd%mat(:,iterm) = vec_fd(ielem)%vec


                !> Reset sdata storage
                call sdata%lin%clear()
                do i=1,size(sdata%rhs)
                    sdata%rhs(i)%vec = ZERO
                end do


            end do


        end associate

    end subroutine







end module mod_check_jacobian
