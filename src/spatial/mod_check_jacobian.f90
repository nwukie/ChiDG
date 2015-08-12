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


            !> Store linearization computed from DNAD
            blk_dnad = sdata%lin%lblks(ielem,iblk)
            blk_fd = blk_dnad   !> Sourced allocation
            blk_fd%mat = ZERO   !> Zero storage


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
                qhold = sdata%q(ielem)%mat(iterm,ivar)
                sdata%q(ielem)%mat(iterm,ivar) = sdata%q(ielem)%mat(iterm,ivar) + eps


                !> For the current element, compute the contributions from volume integrals
                call domain%eqnset%compute_volume_flux(mesh,sdata,ielem,iblk)


                sdata%q(ielem)%mat(iterm,ivar) = qhold  !> Return perturbed value to normal state

                !> Compute finite difference jacobian
                vec_fd = (sdata%rhs - rhs_r)/eps


                !> Store to column of blk_fd
                blk_fd%mat(:,iterm) = vec_fd(ielem)%vec
            end do


        end associate

    end subroutine


















end module mod_check_jacobian
