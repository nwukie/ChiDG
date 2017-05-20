module mod_integrate
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, DIAG, CHIMERA, &
                                  NO_INTERIOR_NEIGHBOR, BOUNDARY, INTERIOR, ZERO
    use mod_chidg_mpi,      only: IRANK

    use type_mesh,          only: mesh_t
    use type_element,       only: element_t
    use type_face,          only: face_t
    use type_element_info,  only: element_info_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t
    use type_solverdata,    only: solverdata_t
    use type_seed,          only: seed_t

    use DNAD_D
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
    !!  @param[in]      ieqn    Index of the variable associated with the flux being integrated
    !!  @param[inout]   flux1   First component of flux at quadrature points
    !!  @param[inout]   flux2   Second component of flux at quadrature points
    !!  @param[inout]   flux3   Third component of flux at quadrature points
    !!
    !--------------------------------------------------------------------------------------------------------
    subroutine integrate_volume_vector_flux(mesh,sdata,elem_info,fcn_info,ieqn,itime,flux1,flux2,flux3)
        type(mesh_t),           intent(in)      :: mesh
        type(solverdata_t),     intent(inout)   :: sdata
        type(element_info_t),   intent(in)      :: elem_info
        type(function_info_t),  intent(in)      :: fcn_info
        integer(ik),            intent(in)      :: ieqn
        integer(ik),            intent(in)      :: itime
        type(AD_D),             intent(inout)   :: flux1(:), flux2(:), flux3(:)


        type(AD_D), dimension(mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%nterms_s)    :: &
            integral, integral1, integral2, integral3


        associate( idom => elem_info%idomain_l, ielem => elem_info%ielement_l, idiff => fcn_info%idiff, &
                   weights     => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%gq%vol%weights, &
                   jinv        => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%jinv,           &
                   grad1_trans => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%grad1_trans,    &
                   grad2_trans => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%grad2_trans,    &
                   grad3_trans => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%grad3_trans )


        !
        ! Multiply each component by quadrature weights and element jacobians
        !
        flux1 = flux1 * weights * jinv
        flux2 = flux2 * weights * jinv
        flux3 = flux3 * weights * jinv


        !
        ! FLUX-1: Multiply by column of test function gradients, integrate
        !
        integral1 = matmul(grad1_trans,flux1)


        !
        ! FLUX-2: Multiply by column of test function gradients, integrate
        !
        integral2 = matmul(grad2_trans,flux2)


        !
        ! FLUX-3: Multiply by column of test function gradients, integrate
        !
        integral3 = matmul(grad3_trans,flux3)



        !
        ! Add componends from each coordiate direction
        !
        integral = integral1 + integral2 + integral3

        
        !
        ! Store integral and derivatives
        !
        call store_volume_integrals(mesh,sdata,elem_info,fcn_info,ieqn,itime,integral)

        end associate

    end subroutine integrate_volume_vector_flux
    !*********************************************************************************************************












    !>  Compute the volume integral of a source function.
    !!
    !!      - Adds value contribution to the rhs vector
    !!      - Adds the derivative contribution to the linearization matrix
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      elem    Element being integrated over
    !!  @param[inout]   rhs     Right-hand side vector storage
    !!  @param[inout]   lin     Domain linearization matrix
    !!  @param[in]      iblk    Selected block of the linearization being computed. lin(ielem,iblk), where iblk = (1-7)
    !!  @param[in]      ieqn    Index of the variable associated with the flux being integrated
    !!  @param[inout]   source  source function to integrate
    !!
    !--------------------------------------------------------------------------------------------------------
    subroutine integrate_volume_scalar_source(mesh,sdata,elem_info,fcn_info,ieqn,itime,source)
        type(mesh_t),           intent(in)      :: mesh
        type(solverdata_t),     intent(inout)   :: sdata
        type(element_info_t),   intent(in)      :: elem_info
        type(function_info_t),  intent(in)      :: fcn_info
        integer(ik),            intent(in)      :: ieqn
        integer(ik),            intent(in)      :: itime
        type(AD_D),             intent(inout)   :: source(:)

        type(AD_D), dimension(mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%nterms_s)    :: integral, integral_x, integral_y, integral_z

        associate( idom => elem_info%idomain_l, ielem => elem_info%ielement_l, idiff => fcn_info%idiff, &
                   weights => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%gq%vol%weights,     &
                   val     => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%gq%vol%val,         &
                   jinv    => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%jinv )

        !
        ! Multiply each component by quadrature weights and element jacobians
        !
        source = source * weights * jinv


        !
        ! Multiply by column of test functions, integrate
        !
        integral = matmul(transpose(val),source)


        !
        ! Store integral and derivatives
        !
        call store_volume_integrals(mesh,sdata,elem_info,fcn_info,ieqn,itime,integral)

        end associate

    end subroutine integrate_volume_scalar_source
    !********************************************************************************************************













    !>  Compute the boundary integral of a flux scalar
    !!
    !!      - Adds value contribution to the rhs vector
    !!      - Adds the derivative contribution to the linearization matrix
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]      face    Face being integrated over
    !!  @param[inout]   rhs     Right-hand side vector storage
    !!  @param[inout]   lin     Domain linearization matrix
    !!  @param[in]      iblk    Selected block of the linearization being computed. lin(ielem,iblk), where iblk = (1-7)
    !!  @param[in]      ieqn    Index of the variable associated with the flux being integrated
    !!  @param[inout]   flux    function to integrate over the boundary
    !!
    !--------------------------------------------------------------------------------------------------------
    subroutine integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,ieqn,itime,integrand)
        type(mesh_t),           intent(in)      :: mesh
        type(solverdata_t),     intent(inout)   :: sdata
        type(face_info_t),      intent(in)      :: face_info
        type(function_info_t),  intent(in)      :: function_info
        integer(ik),            intent(in)      :: ieqn
        integer(ik),            intent(in)      :: itime
        type(AD_D),             intent(inout)   :: integrand(:)
        
        ! Data for applying to self and neighbor
        type(AD_D), allocatable     :: integral(:)
        type(AD_D), allocatable     :: integrand_n(:)
        type(face_info_t)           :: face_n
        type(function_info_t)       :: function_n
        integer(ik)                 :: ineighbor_element_l, ineighbor_face, ineighbor_proc, idiff_n, ierr
        logical                     :: parallel_neighbor, diff_none, diff_interior, diff_exterior



        associate ( idomain_l   => face_info%idomain_l,     &
                    ielement_l  => face_info%ielement_l,    &
                    iface       => face_info%iface,         &
                    ifcn        => function_info%ifcn,      &
                    idonor      => function_info%idepend,   &
                    idiff       => function_info%idiff )



        ! Neighbor indices
        ineighbor_proc      = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_proc
        ineighbor_element_l = mesh%domain(idomain_l)%faces(ielement_l,iface)%get_neighbor_element_l()
        ineighbor_face      = mesh%domain(idomain_l)%faces(ielement_l,iface)%get_neighbor_face()

        parallel_neighbor = ( IRANK /= ineighbor_proc )


        !
        ! Store quadrature flux for neighbor integral
        !
        integrand_n = integrand



        !
        ! Integrate and apply once
        !
        associate ( weights  => mesh%domain(idomain_l)%faces(ielement_l,iface)%gq%face%weights(:,iface),   &
                    jinv     => mesh%domain(idomain_l)%faces(ielement_l,iface)%jinv,                       &
                    val      => mesh%domain(idomain_l)%faces(ielement_l,iface)%gq%face%val(:,:,iface),     &
                    valtrans => mesh%domain(idomain_l)%faces(ielement_l,iface)%gq%face%val_trans(:,:,iface) )


            !
            ! Multiply each component by quadrature weights. The fluxes have already been multiplied by norm
            !
            integrand = integrand * weights
            integral = matmul(valtrans,integrand)


            call store_boundary_integral_residual(     mesh,sdata,face_info,function_info,ieqn,itime,integral)
            call store_boundary_integral_linearization(mesh,sdata,face_info,function_info,ieqn,itime,integral)


        end associate





        !
        ! Integrate and apply second time if there is a neighbor
        !
        if ( ineighbor_element_l /= NO_INTERIOR_NEIGHBOR ) then
            if ( .not. parallel_neighbor ) then

                face_n%idomain_g  = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_g
                face_n%idomain_l  = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_l
                face_n%ielement_g = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_g
                face_n%ielement_l = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_l
                face_n%iface      = mesh%domain(idomain_l)%faces(ielement_l,iface)%get_neighbor_face()

                !
                ! Get linearization block for the neighbor element
                !
                diff_none     = (idiff == 0)
                diff_interior = (idiff == DIAG)
                diff_exterior = ( (idiff == 1) .or. (idiff == 2) .or. &
                                  (idiff == 3) .or. (idiff == 4) .or. &
                                  (idiff == 5) .or. (idiff == 6) )
                if ( diff_exterior ) then
                        idiff_n = DIAG
                else if ( diff_interior ) then
                        idiff_n = ineighbor_face
                else if ( diff_none ) then
                        idiff_n = 0
                else
                    call chidg_signal(FATAL,"store_boundary_integrals: unexpected value")
                end if


                function_n%type    = function_info%type
                function_n%ifcn    = function_info%ifcn
                function_n%idepend = function_info%idepend
                function_n%seed    = function_info%seed
                function_n%idiff   = idiff_n


                associate ( weights_n  => mesh%domain(idomain_l)%faces(ineighbor_element_l,ineighbor_face)%gq%face%weights(:,ineighbor_face),   &
                            jinv_n     => mesh%domain(idomain_l)%faces(ineighbor_element_l,ineighbor_face)%jinv,                                & 
                            val_n      => mesh%domain(idomain_l)%faces(ineighbor_element_l,ineighbor_face)%gq%face%val(:,:,ineighbor_face),     &
                            valtrans_n => mesh%domain(idomain_l)%faces(ineighbor_element_l,ineighbor_face)%gq%face%val_trans(:,:,ineighbor_face) )

                    integrand_n = integrand_n * weights_n

                    !
                    ! Integrate and negate for contribution to neighbor element
                    !
                    integral = -matmul(valtrans_n,integrand_n)

                    call store_boundary_integral_residual(     mesh,sdata,face_n,function_n,ieqn,itime,integral)
                    call store_boundary_integral_linearization(mesh,sdata,face_n,function_n,ieqn,itime,integral)

                end associate

            end if ! .not. parallel_neighbor
        end if ! ielem_n


        end associate

    end subroutine integrate_boundary_scalar_flux
    !********************************************************************************************************










    !> Store volume integral values to RHS vector, and partial derivatives to LIN block matrix
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!  @param[in]      integral    Array of autodiff values containing integrals and partial derivatives for the RHS vector and LIN linearization matrix
    !!  @param[inout]   rhs         Right-hand side vector
    !!  @param[inout]   lin         Block matrix storing the linearization of the spatial scheme
    !!  @param[in]      ielem       Element index for applying to the correct location in RHS and LIN
    !!  @param[in]      ieqn        Variable index
    !!  @param[in]      idiff       Block index for the correct linearization block for the current element
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !---------------------------------------------------------------------------------------------------------
    subroutine store_volume_integrals(mesh,sdata,elem_info,fcn_info,ieqn,itime,integral)
        type(mesh_t),           intent(in)      :: mesh
        type(solverdata_t),     intent(inout)   :: sdata
        type(element_info_t),   intent(in)      :: elem_info
        type(function_info_t),  intent(in)      :: fcn_info
        integer(ik),            intent(in)      :: ieqn
        integer(ik),            intent(in)      :: itime
        type(AD_D),             intent(inout)   :: integral(:)

        integer(ik)         :: i
        logical             :: conforming_face, boundary_face, chimera_face, &
                               diff_none, diff_interior, diff_exterior
        type(face_info_t)   :: face_info
        real(rk)            :: vals(size(integral))

        associate ( idom  => elem_info%idomain_l,   &
                    ielem => elem_info%ielement_l,  &
                    idiff => fcn_info%idiff )


        diff_none     = (idiff == 0)
        diff_interior = (idiff == DIAG)
        diff_exterior = ( (idiff == 1) .or. (idiff == 2) .or. &
                          (idiff == 3) .or. (idiff == 4) .or. &
                          (idiff == 5) .or. (idiff == 6) )


        !
        ! Only store rhs once. 
        !   - If we are differentiating things, only apply function once, when wrt interior
        !       idiff == DIAG
        !   - If we are not differentiating things, just apply function.
        !       idiff == 0
        !
        if ( diff_interior .or. diff_none ) then
            vals = sdata%rhs%dom(idom)%vecs(ielem)%getvar(ieqn,itime) - integral(:)%x_ad_
            call sdata%rhs%dom(idom)%vecs(ielem)%setvar(ieqn,itime,vals)
        end if

        !
        ! Negate derivatives before adding to linearization
        !
        do i = 1,size(integral)
            integral(i)%xp_ad_ = -integral(i)%xp_ad_
        end do



        !
        ! Check if linearization is with respect to an exterior element. 
        ! Only need this for diff_exterior, and idiff is undefined as a face index for 
        ! diff_interior or diff_none.
        !
        if ( diff_exterior ) then
            conforming_face = (mesh%domain(idom)%faces(ielem,idiff)%ftype == INTERIOR)
            boundary_face   = (mesh%domain(idom)%faces(ielem,idiff)%ftype == BOUNDARY)
            chimera_face    = (mesh%domain(idom)%faces(ielem,idiff)%ftype == CHIMERA )
        end if


        !
        ! Initialize a face_info in case it is needed below
        !
        face_info%idomain_g  = elem_info%idomain_g
        face_info%idomain_l  = elem_info%idomain_l
        face_info%ielement_g = elem_info%ielement_g
        face_info%ielement_l = elem_info%ielement_l
        face_info%iface      = idiff


        !
        ! Store linearization
        !
        if ( diff_interior ) then
            call sdata%lhs%store(integral,face_info,fcn_info%seed,ieqn,itime)

        else if ( diff_exterior .and. conforming_face ) then
            call sdata%lhs%store(integral,face_info,fcn_info%seed,ieqn,itime)
        else if ( diff_exterior .and. chimera_face    ) then
            call sdata%lhs%store_chimera(integral,face_info,fcn_info%seed,ieqn,itime)
        else if ( diff_exterior .and. boundary_face   ) then
            call sdata%lhs%store_bc(integral,face_info,fcn_info%seed,ieqn,itime)

        else if ( diff_none ) then
            ! No derivatives to store
        else
            call chidg_signal(FATAL,"store_volume_integrals: Invalid condition for storing integrals. Could be a bad face type or linearization direction")
        end if


        end associate

    end subroutine store_volume_integrals
    !*********************************************************************************************************










    !> Store boundary integral values to RHS vector, and partial derivatives to LIN block matrix
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!  @param[in]      integral    Array of autodiff values containing integrals and partial derivatives for the RHS vector and LIN linearization matrix
    !!  @param[inout]   rhs         Right-hand side vector
    !!  @param[inout]   lin         Block matrix storing the linearization of the spatial scheme
    !!  @param[in]      ielem       Element index for applying to the correct location in RHS and LIN
    !!  @param[in]      ieqn        Variable index
    !!  @param[in]      iblk        Block index for the correct linearization block for the current element
    !!
    !!  @author Mayank Sharma + matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !--------------------------------------------------------------------------------------------------------
    subroutine store_boundary_integral_residual(mesh,sdata,face_info,function_info,ieqn,itime,integral)
        type(mesh_t),           intent(in)      :: mesh
        type(solverdata_t),     intent(inout)   :: sdata
        type(face_info_t),      intent(in)      :: face_info
        type(function_info_t),  intent(in)      :: function_info
        integer(ik),            intent(in)      :: ieqn
        integer(ik),            intent(in)      :: itime
        type(AD_D),             intent(inout)   :: integral(:)

        real(rk)        :: vals(size(integral))

        logical         :: boundary_face, chimera_face, conforming_face, diff_interior, diff_none, diff_exterior
        logical         :: add_flux = .false.


        associate ( idomain_l  => face_info%idomain_l, ielement_l  => face_info%ielement_l, iface => face_info%iface, &
                    ifcn  => function_info%ifcn,       idepend => function_info%idepend,     idiff => function_info%idiff )

            conforming_face = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR)
            boundary_face   = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == BOUNDARY)
            chimera_face    = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA )

            diff_none     = (idiff == 0)
            diff_interior = (idiff == DIAG)
            diff_exterior = ( (idiff == 1) .or. (idiff == 2) .or. &
                              (idiff == 3) .or. (idiff == 4) .or. &
                              (idiff == 5) .or. (idiff == 6) )

            associate ( rhs => sdata%rhs%dom(idomain_l)%vecs, lhs => sdata%lhs)

                !
                ! Only store rhs once. 
                !
                !   For BOUNDARY faces: only store if diff_exterior or diff_none, since the boundary integrals only get computed
                !                       for these cases, not diff_interior actually. This is because the interior element is 
                !                       registered with the boundary condition as a coupled element, and these get computed for
                !                       icompute = [iface] and not icompute = [DIAG]. Only store for idepend == 1, since 
                !                       could be computed multiple times.
                !
                !   For CHIMERA faces: only store if diff_interior or diff_none. Only store for idepend == 1, since could be
                !                      computed multiple times. The residual should be the same for any value of idepend, 
                !                      only the derivatives will change.
                !
                !if (  (boundary_face .and. diff_interior) .or. &
                if (  (boundary_face .and. diff_exterior) .or. &
                      (boundary_face .and. diff_none    ) ) then

                    if (idepend == 1) then
                        vals = rhs(ielement_l)%getvar(ieqn,itime) + integral(:)%x_ad_
                        call rhs(ielement_l)%setvar(ieqn,itime,vals)
                    end if


                else if ( (chimera_face .and. diff_interior) .or. &
                          (chimera_face .and. diff_none    ) ) then

                    if (idepend == 1) then
                        vals = rhs(ielement_l)%getvar(ieqn,itime) + integral(:)%x_ad_
                        call rhs(ielement_l)%setvar(ieqn,itime,vals)
                    end if


                else if ( conforming_face ) then

                    !
                    ! Check if particular flux function has been added already
                    !
                    add_flux = sdata%function_status%compute_function_equation( face_info, function_info, ieqn)


                    ! Store if needed
                    if ( add_flux ) then
                        ! Add to residual and store
                        vals = rhs(ielement_l)%getvar(ieqn,itime) + integral(:)%x_ad_
                        call rhs(ielement_l)%setvar(ieqn,itime,vals)

                        ! Register flux was stored
                        call sdata%function_status%register_function_computed( face_info, function_info, ieqn)
                    end if

                end if


            end associate
        end associate

    end subroutine store_boundary_integral_residual
    !***********************************************************************************************************









    !> Store boundary integral values to RHS vector, and partial derivatives to LIN block matrix
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!  @param[in]      integral    Array of autodiff values containing integrals and partial derivatives for the RHS vector and LIN linearization matrix
    !!  @param[inout]   rhs         Right-hand side vector
    !!  @param[inout]   lin         Block matrix storing the linearization of the spatial scheme
    !!  @param[in]      ielem       Element index for applying to the correct location in RHS and LIN
    !!  @param[in]      ieqn        Variable index
    !!  @param[in]      idiff       Block index for the correct linearization block for the current element
    !!
    !--------------------------------------------------------------------------------------------------------
    subroutine store_boundary_integral_linearization(mesh,sdata,face_info,function_info,ieqn,itime,integral)
        type(mesh_t),           intent(in)      :: mesh
        type(solverdata_t),     intent(inout)   :: sdata
        type(face_info_t),      intent(in)      :: face_info
        type(function_info_t),  intent(in)      :: function_info
        integer(ik),            intent(in)      :: ieqn
        integer(ik),            intent(in)      :: itime
        type(AD_D),             intent(inout)   :: integral(:)

        integer(ik)                 :: i, idomain_l, ielement_l, iface, ChiID
        integer(ik)                 :: idiff, ifcn
        real(rk)                    :: vals(size(integral))

        logical :: conforming_face, boundary_face, chimera_face, &
                   diff_none, diff_interior, diff_exterior, add_linearization

        associate ( idomain_l => face_info%idomain_l, ielement_l => face_info%ielement_l,  iface => face_info%iface, &
                    ifcn      => function_info%ifcn,  idepend    => function_info%idepend, idiff => function_info%idiff, seed => function_info%seed )

        
        conforming_face = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR)
        boundary_face   = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == BOUNDARY)
        chimera_face    = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA )

        diff_none     = ( idiff == 0 )
        diff_interior = ( idiff == DIAG )
        diff_exterior = ( (idiff == 1) .or. (idiff == 2) .or. &
                          (idiff == 3) .or. (idiff == 4) .or. &
                          (idiff == 5) .or. (idiff == 6) )

        associate ( rhs => sdata%rhs%dom(idomain_l)%vecs, lhs => sdata%lhs)

        if (diff_interior .or. diff_exterior) then

            ! Store linearization. Rules for different face types.
            if ( chimera_face ) then

                if (diff_exterior) then
                    ! Store linearization of Chimera boundary donor elements.
                    call lhs%store_chimera(integral,face_info,seed,ieqn,itime)
                else
                    ! Store linearization of Chimera boundary receiver element. Since this could be computed multiple times,
                    ! we just store it once.
                    if (idepend == 1) then
                        call lhs%store(integral,face_info,seed,ieqn,itime)
                    end if
                end if



            else if ( boundary_face ) then
                call lhs%store_bc(integral,face_info,seed,ieqn,itime)



            else if ( conforming_face ) then


                add_linearization = sdata%function_status%linearize_function_equation( face_info, function_info, ieqn)

                ! Store linearization if not already stored
                if ( add_linearization ) then
                    ! Store linearization
                    call lhs%store(integral,face_info,seed,ieqn,itime)

                    ! Register flux as linearized
                    call sdata%function_status%register_function_linearized( face_info, function_info, ieqn)
                end if


            end if ! ftype

        end if ! diff_...
        end associate

        end associate 

    end subroutine store_boundary_integral_linearization
    !***********************************************************************************************************









end module mod_integrate
