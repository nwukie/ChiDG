module mod_integrate
    use mod_kinds,  only: rk,ik



    implicit none



contains



    subroutine integrate_volume_flux(elem,rhs,A,iblk,ivar,flux_x,flux_y,flux_z)
        type(element_t),        intent(in)      :: elem
        type(expansion_t),      intent(inout)   :: rhs(:)
        type(blockmatrix_t),    intent(inout)   :: A(:,:)        !> Linearization matrix
        integer(ik),            intent(in)      :: ielem
        integer(ik),            intent(in)      :: ivar
        type(AD_D),             intent(inout)   :: flux_x(:), flux_y(:), flux_z(:)


        type(AD_D), dimension(elem%nterms_s) :: integral

        ! Multiply each component by quadrature weights and element jacobians
        flux_x = (flux_x) * (self%gq%vol%weights) * (self%jinv)
        flux_y = (flux_y) * (self%gq%vol%weights) * (self%jinv)
        flux_z = (flux_z) * (self%gq%vol%weights) * (self%jinv)

        ! Multiply by column of test function gradients, integrate, and add to RHS
        integral = matmul(transpose(self%dtdx),flux_x)
        rhs(ielem)%mat(:,ivar) = rhs(ielem)%mat(:,ivar) + integral(:)%x_ad_







        integral = matmul(transpose(self%dtdy),flux_y)
!        self%rhs%vals(:,varindex) = self%rhs%vals(:,varindex) + integral

        integral = matmul(transpose(self%dtdz),flux_z)
!        self%rhs%vals(:,varindex) = self%rhs%vals(:,varindex) + integral


    end subroutine




end module mod_integrate
