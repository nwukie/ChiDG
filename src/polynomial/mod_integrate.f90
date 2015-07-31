module mod_integrate
    use mod_kinds,  only: rk,ik



    implicit none



contains



    !>  Compute the volume integral of a flux vector
    !!
    !!      - Adds value contribution to the rhs vector
    !!      - Adds the derivative contribution to the linearization matrix
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]      elem    Element being integrated over
    !!  @param[inout]   rhs     Right-hand side vector storage
    !!  @param[inout]   lin     Domain linearization matrix
    !!  @param[in]      iblk    Selected block of the linearization being computed. lin(ielem,iblk), where iblk = (1-7)
    !!  @param[in]      ivar    Index of the variable associated with the flux being integrated
    !!  @param[inout]   flux_x  x-Flux and derivatives at quadrature points
    !!  @param[inout]   flux_y  y-Flux and derivatives at quadrature points
    !!  @param[inout]   flux_z  z-Flux and derivatives at quadrature points
    !--------------------------------------------------------------------------------------------------------
    subroutine integrate_volume_flux(elem,rhs,lin,iblk,ivar,flux_x,flux_y,flux_z)
        type(element_t),        intent(in)      :: elem
        type(expansion_t),      intent(inout)   :: rhs(:)
        type(blockmatrix_t),    intent(inout)   :: lin
        integer(ik),            intent(in)      :: iblk
        integer(ik),            intent(in)      :: ivar
        type(AD_D),             intent(inout)   :: flux_x(:), flux_y(:), flux_z(:)


        integer(ik)                             :: ielem
        type(AD_D), dimension(elem%nterms_s)    :: integral

        ielem = elem%ielem  !> get element index

        ! Multiply each component by quadrature weights and element jacobians
        flux_x = (flux_x) * (self%gq%vol%weights) * (self%jinv)
        flux_y = (flux_y) * (self%gq%vol%weights) * (self%jinv)
        flux_z = (flux_z) * (self%gq%vol%weights) * (self%jinv)

        ! Multiply by column of test function gradients, integrate, and add to RHS
        integral = matmul(transpose(self%dtdx),flux_x)
        rhs(ielem)%mat(:,ivar) = rhs(ielem)%mat(:,ivar) + integral(:)%x_ad_

        ! loop over terms in the residual expansion and store derivatives in the linearization matrix
        call lin%store(integral,ielem,iblk,ivar)





        integral = matmul(transpose(self%dtdy),flux_y)
!        self%rhs%vals(:,varindex) = self%rhs%vals(:,varindex) + integral

        integral = matmul(transpose(self%dtdz),flux_z)
!        self%rhs%vals(:,varindex) = self%rhs%vals(:,varindex) + integral


    end subroutine




end module mod_integrate
