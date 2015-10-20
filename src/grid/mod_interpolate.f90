module mod_interpolate
    use mod_kinds,      only: rk,ik
    use DNAD_D
    use type_mesh,      only: mesh_t
    !use type_element,   only: element_t
    !use type_face,      only: face_t
    !use type_blockvector,   only: blockvector_t
    use type_chidgVector,   only: chidgVector_t

    implicit none




    interface interpolate_element
        module procedure    interpolate_element_autodiff, interpolate_element_standard
    end interface

    interface interpolate_face
        module procedure    interpolate_face_autodiff,    interpolate_face_standard
    end interface





contains


    !> Compute variable at quadrature nodes.
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------
    subroutine interpolate_element_autodiff(mesh,q,idom,ielem,ivar,var_gq,ielem_seed)
        type(mesh_t),        intent(in)      :: mesh(:)
        type(chidgVector_t), intent(inout)   :: q
        integer(ik),         intent(in)      :: idom
        integer(ik),         intent(in)      :: ielem
        integer(ik),         intent(in)      :: ivar
        type(AD_D),          intent(inout)   :: var_gq(:)
        integer(ik),         intent(in)      :: ielem_seed

        type(AD_D)  :: qdiff(mesh(idom)%elems(ielem)%nterms_s)
        integer(ik) :: nderiv, set_deriv, iterm, igq, i


        associate ( elems => mesh(idom)%elems )


            !
            ! Get the number of degrees of freedom for the seed element
            ! and set this as the number of partial derivatives to track
            !
            if (ielem_seed == 0) then
                ! If ielem_seed == 0 then we aren't interested in tracking derivatives
                nderiv = 1
            else
                nderiv = elems(ielem_seed)%neqns  *  elems(ielem_seed)%nterms_s
            end if



            !
            ! Allocate the derivative array for each autodiff variable
            ! MIGHT NOT NEED THIS IF IT GETS AUTOMATICALLY ALLOCATED ON ASSIGNMENT -- TEST
            !
            do igq = 1,size(var_gq)
                var_gq(igq) = AD_D(nderiv)
            end do



            !
            ! If the current element is being differentiated (ielem == ielem_seed)
            ! then copy the solution modes to local AD variable and seed derivatives
            !
            if (ielem == ielem_seed) then

                !
                ! Allocate derivative arrays for temporary solution variable
                !
                do iterm = 1,elems(ielem)%nterms_s
                    qdiff(iterm) = AD_D(nderiv)
                end do


                !
                ! Copy the solution variables from 'q' to 'qdiff'
                !
                qdiff = q%dom(idom)%lvecs(ielem)%getvar(ivar)


                !
                ! Loop through the terms in qdiff
                !
                do iterm = 1,size(qdiff)
                    !
                    ! For the given term, seed its appropriate derivative
                    !
                    set_deriv = (ivar - 1)*elems(ielem)%nterms_s + iterm
                    qdiff(iterm)%xp_ad_(set_deriv) = 1.0_rk
                end do


                var_gq = matmul(elems(ielem)%gq%vol%val,qdiff)

            else
                !
                ! If the solution variable derivatives dont need initialized
                ! then just use the q(ielem) values and derivatives get
                ! initialized to zero
                !
                var_gq = matmul(elems(ielem)%gq%vol%val,q%dom(idom)%lvecs(ielem)%getvar(ivar))
            end if

        end associate

    end subroutine










    !> Compute variable at quadrature nodes.
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------
    subroutine interpolate_face_autodiff(mesh,q,idom,ielem,iface,ivar,var_gq,ielem_seed)
        type(mesh_t),           intent(in)      :: mesh(:)
        type(chidgVector_t),    intent(inout)   :: q
        integer(ik),            intent(in)      :: idom
        integer(ik),            intent(in)      :: ielem
        integer(ik),            intent(in)      :: iface
        integer(ik),            intent(in)      :: ivar
        type(AD_D),             intent(inout)   :: var_gq(:)
        integer(ik),            intent(in)      :: ielem_seed

        type(AD_D)  :: qdiff(mesh(idom)%faces(ielem,iface)%nterms_s)
        integer(ik) :: nderiv, set_deriv, iterm, igq, nterms_s


        associate ( faces => mesh(idom)%faces )


            nterms_s = faces(ielem,iface)%nterms_s

            !
            ! Get the number of degrees of freedom for the seed element
            ! and set this as the number of partial derivatives to track
            !
            if (ielem_seed == 0) then
                ! If ielem_seed == 0 then we aren't interested in tracking derivatives
                nderiv = 1
            else
                nderiv = faces(ielem_seed,1)%neqns  *  faces(ielem_seed,1)%nterms_s     ! using face 1 here, but faces 1-6 point
            end if



            !
            ! Allocate the derivative array for each autodiff variable
            ! MIGHT NOT NEED THIS IF IT GETS AUTOMATICALLY ALLOCATED ON ASSIGNMENT -- TEST
            !
            do igq = 1,size(var_gq)
                allocate(var_gq(igq)%xp_ad_(nderiv))
            end do

            !
            ! If the current element is being differentiated (ielem == ielem_seed)
            ! then copy the solution modes to local AD variable and seed derivatives
            !
            if (ielem == ielem_seed) then

                !
                ! Allocate derivative arrays for temporary solution variable
                !
                do iterm = 1,nterms_s
                    allocate(qdiff(iterm)%xp_ad_(nderiv))
                end do

                !
                ! Copy the solution variables from 'q' to 'qdiff'
                !
                qdiff = q%dom(idom)%lvecs(ielem)%getvar(ivar)


                !
                ! Loop through the terms in qdiff
                !
                do iterm = 1,size(qdiff)
                    ! For the given term, seed its appropriate derivative
                    set_deriv = (ivar - 1)*nterms_s + iterm
                    qdiff(iterm)%xp_ad_(set_deriv) = 1.0_rk
                end do

                var_gq = matmul(faces(ielem,iface)%gq%face%val(:,:,iface),  qdiff)

            else
                !
                ! If the solution variable derivatives dont need initialized
                ! then just use the q(ielem) values and derivatives get
                ! initialized to zero
                !
                var_gq = matmul(faces(ielem,iface)%gq%face%val(:,:,iface),  q%dom(idom)%lvecs(ielem)%getvar(ivar))
            end if

        end associate

    end subroutine













    !> Compute variable at quadrature nodes.
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------
    subroutine interpolate_element_standard(mesh,q,idom,ielem,ivar,var_gq)
        type(mesh_t),           intent(in)      :: mesh(:)
        type(chidgVector_t),    intent(inout)   :: q
        integer(ik),            intent(inout)   :: idom
        integer(ik),            intent(in)      :: ielem
        integer(ik),            intent(in)      :: ivar
        real(rk),               intent(inout)   :: var_gq(:)

        !
        ! Use quadrature instance to compute variable at quadrature nodes.
        ! This takes the form of a matrix multiplication of the quadrature matrix
        ! with the array of modes for the given variable.
        !
        var_gq = matmul(mesh(idom)%elems(ielem)%gq%vol%val, q%dom(idom)%lvecs(ielem)%getvar(ivar))

    end subroutine













    !> Compute variable at quadrature nodes.
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------
    subroutine interpolate_face_standard(mesh,q,idom,ielem,iface,ivar,var_gq)
        type(mesh_t),           intent(in)      :: mesh(:)
        type(chidgVector_t),    intent(inout)   :: q
        integer(ik),            intent(in)      :: idom, ielem, iface, ivar
        real(rk),               intent(inout)   :: var_gq(:)


        !
        ! Use quadrature instance to compute variable at quadrature nodes.
        ! This takes the form of a matrix multiplication of the face quadrature matrix
        ! with the array of modes for the given variable
        !
        var_gq = matmul(mesh(idom)%faces(ielem,iface)%gq%face%val(:,:,iface), q%dom(idom)%lvecs(ielem)%getvar(ivar))



    end subroutine
















end module mod_interpolate







