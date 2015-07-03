module mod_interpolate
    use mod_kinds,      only: rk,ik
    use DNAD_D
    use type_element,   only: element_t
    use type_face,      only: face_t
    use atype_solver,   only: q

    implicit none

    interface interpolate
        module procedure    interpolate_autodiff, interpolate_standard
    end interface

    contains


    !> Compute variable at quadrature nodes.
    !!
    !!
    !----------------------------------------------------------------
    subroutine interpolate_autodiff(elems,ielem,ivar,var_gq,ielem_seed)
        class(element_t),   intent(in)      :: elems(:)
        integer(ik),        intent(in)      :: ielem
        integer(ik),        intent(in)      :: ivar
        type(AD_D),         intent(inout)   :: var_gq(:)
        integer(ik),        intent(in)      :: ielem_seed

        type(AD_D)  :: qdiff(elems(ielem)%nterms_s)
        integer(ik) :: nderiv, set_deriv, iterm, igq


        !> Get the number of degrees of freedom for the seed element
        !! and set this as the number of partial derivatives to track
        nderiv = elems(ielem_seed)%neqns  *  elems(ielem_seed)%nterms_s

        !> Allocate the derivative array for each autodiff variable
        !! MIGHT NOT NEED THIS IF IT GETS AUTOMATICALLY ALLOCATED ON ASSIGNMENT -- TEST
        do igq = 1,size(var_gq)
            var_gq(igq) = AD_D(nderiv)
        end do

        !> If the current element is being differentiated (ielem == ielem_seed)
        !! then copy the solution modes to local AD variable and seed derivatives
        if (ielem == ielem_seed) then

            !> Allocate derivative arrays for temporary solution variable
            do iterm = 1,elems(ielem)%nterms_s
                qdiff(iterm) = AD_D(nderiv)
            end do

            !> Copy the solution variables from 'q' to 'qdiff'
            qdiff = q(ielem)%var(ivar)

            !> Loop through the terms in qdiff
            do iterm = 1,size(qdiff)
                !> For the given term, seed its appropriate derivative
                set_deriv = (ivar - 1)*elems(ielem)%nterms_s + iterm
                qdiff(iterm)%xp_ad_(set_deriv) = 1.0_rk
            end do

            var_gq = matmul(elems(ielem)%gq%vol%val,qdiff)

        else
            !> If the solution variable derivatives dont need initialized
            !! then just use the q(ielem) values and derivatives get
            !! initialized to zero
            var_gq = matmul(elems(ielem)%gq%vol%val,q(ielem)%var(ivar))
        end if
    end subroutine







    subroutine interpolate_standard(elems)
        class(element_t),   intent(inout)   :: elems(:)

    end subroutine




end module mod_interpolate
