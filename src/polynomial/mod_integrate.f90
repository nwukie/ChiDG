module mod_integrate
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, DIAG, CHIMERA, &
                                      NO_INTERIOR_NEIGHBOR, BOUNDARY, INTERIOR, ZERO, NO_DATA, NO_ID, &
                                      dQ_DIFF, dX_DIFF, dBC_DIFF, NO_DIFF, dD_DIFF
    use mod_chidg_mpi,          only: IRANK

    use type_mesh,              only: mesh_t
    use type_element,           only: element_t
    use type_face,              only: face_t
    use type_element_info,      only: element_info_t, element_info
    use type_function_info,     only: function_info_t
    use type_solverdata,        only: solverdata_t
    use type_seed,              only: seed_t
    use mod_differentiate,      only: differentiate_jinv, differentiate_element_interpolator, &
                                      differentiate_face_interior_interpolator,               &
                                      differentiate_normal
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
    !----------------------------------------------------------------------------------------
    subroutine integrate_volume_vector_flux(mesh,sdata,elem_info,fcn_info,ifield,itime,flux1,flux2,flux3)
        type(mesh_t),           intent(in)      :: mesh
        type(solverdata_t),     intent(inout)   :: sdata
        type(element_info_t),   intent(in)      :: elem_info
        type(function_info_t),  intent(in)      :: fcn_info
        integer(ik),            intent(in)      :: ifield
        integer(ik),            intent(in)      :: itime
        type(AD_D),             intent(inout)   :: flux1(:), flux2(:), flux3(:)


        type(AD_D), dimension(mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%nterms_s)    :: &
            integral, integral1, integral2, integral3
        type(AD_D), allocatable :: diff_interp_1(:,:), diff_interp_2(:,:), diff_interp_3(:,:)


        associate( idom    => elem_info%idomain_l, ielem => elem_info%ielement_l, idiff => fcn_info%idiff,             &
                   weights => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%basis_s%weights_element(),   &
                   jinv    => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%jinv,          &
                   grad1_trans => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%grad1_trans,         &
                   grad2_trans => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%grad2_trans,         &
                   grad3_trans => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%grad3_trans )

        ! If dQ differentiation or no differentiation
        if (fcn_info%dtype == dQ_DIFF  .or. &
            fcn_info%dtype == dBC_DIFF .or. &
            fcn_info%dtype == NO_DIFF  .or. &
            fcn_info%dtype == dD_DIFF) then

            ! Multiply each component by quadrature weights and element jacobians
            flux1 = flux1 * weights * jinv
            flux2 = flux2 * weights * jinv
            flux3 = flux3 * weights * jinv

            ! Multiply by column of test function gradients, integrate
            integral1 = matmul(grad1_trans,flux1)
            integral2 = matmul(grad2_trans,flux2)
            integral3 = matmul(grad3_trans,flux3)


        ! If dX differentiation
        else if (fcn_info%dtype == dX_DIFF) then

            ! Multiply each component by quadrature weights and element jacobians
            flux1 = flux1 * weights * differentiate_jinv(mesh,elem_info,0,fcn_info,'element')
            flux2 = flux2 * weights * differentiate_jinv(mesh,elem_info,0,fcn_info,'element')
            flux3 = flux3 * weights * differentiate_jinv(mesh,elem_info,0,fcn_info,'element')

            ! FLUX-1: Multiply by column of test function gradients, integrate
            diff_interp_1 = transpose(differentiate_element_interpolator('grad1',mesh,elem_info,fcn_info,itime))
            integral1 = matmul(diff_interp_1,flux1)
            ! The following statement creates issues on norman
            !integral1 = matmul(transpose(differentiate_element_interpolator('grad1',mesh,elem_info,fcn_info)),flux1)

            ! FLUX-2: Multiply by column of test function gradients, integrate
            diff_interp_2 = transpose(differentiate_element_interpolator('grad2',mesh,elem_info,fcn_info,itime))
            integral2 = matmul(diff_interp_2,flux2)
            !integral2 = matmul(transpose(differentiate_element_interpolator('grad2',mesh,elem_info,fcn_info)),flux2)

            ! FLUX-3: Multiply by column of test function gradients, integrate
            diff_interp_3 = transpose(differentiate_element_interpolator('grad3',mesh,elem_info,fcn_info,itime))
            integral3 = matmul(diff_interp_3,flux3)
            !integral3 = matmul(transpose(differentiate_element_interpolator('grad3',mesh,elem_info,fcn_info)),flux3)

        else

            call chidg_signal_one(FATAL,"integrate_volume_vector_flux: unexpected differentiate type",fcn_info%dtype)

        end if !fcn_info%dtype






        ! Add componends from each coordiate direction
        integral = integral1 + integral2 + integral3

        ! Store integral and derivatives
        call store_volume_integrals(mesh,sdata,elem_info,fcn_info,ifield,itime,integral)

        end associate

    end subroutine integrate_volume_vector_flux
    !*****************************************************************************************












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
    !!  @param[in]      ifield  Index of the variable associated with the flux being integrated
    !!  @param[inout]   source  source function to integrate
    !!
    !----------------------------------------------------------------------------------------
    subroutine integrate_volume_scalar_source(mesh,sdata,elem_info,fcn_info,ifield,itime,source_in)
        type(mesh_t),           intent(in)      :: mesh
        type(solverdata_t),     intent(inout)   :: sdata
        type(element_info_t),   intent(in)      :: elem_info
        type(function_info_t),  intent(in)      :: fcn_info
        integer(ik),            intent(in)      :: ifield
        integer(ik),            intent(in)      :: itime
        type(AD_D),             intent(inout)   :: source_in(:)

        type(AD_D), dimension(mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%nterms_s)    :: integral, integral_x, integral_y, integral_z
        type(AD_D), allocatable :: source(:)

        associate( idom    => elem_info%idomain_l, ielem => elem_info%ielement_l, idiff => fcn_info%idiff,                          &
                   weights => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%basis_s%weights_element(),               &
                   val     => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%basis_s%interpolator_element('Value'),   &
                   jinv    => mesh%domain(elem_info%idomain_l)%elems(elem_info%ielement_l)%jinv )


        ! If dQ differentiation or no differentiation
        if (fcn_info%dtype == dQ_DIFF .or. fcn_info%dtype == dBC_DIFF .or. fcn_info%dtype == NO_DIFF .or. fcn_info%dtype == dD_DIFF) then
            ! Multiply each component by quadrature weights and element jacobians
            source = source_in * weights * jinv

        ! If dX differentiation
        else if (fcn_info%dtype == dX_DIFF) then
            ! Multiply each component by quadrature weights and element jacobians
            source = source * weights * differentiate_jinv(mesh,elem_info,0,fcn_info,'element')

        else
            call chidg_signal_one(FATAL,"integrate_volume_vector_source: unexpected differentiate type",fcn_info%dtype)
        
        end if


        ! Multiply by column of test functions, integrate
        integral = matmul(transpose(val),source)

        ! Store integral and derivatives
        call store_volume_integrals(mesh,sdata,elem_info,fcn_info,ifield,itime,integral)

        end associate

    end subroutine integrate_volume_scalar_source
    !****************************************************************************************







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
    !!  @param[in]      ifield  Index of the variable associated with the flux being integrated
    !!  @param[inout]   flux    function to integrate over the boundary
    !!
    !----------------------------------------------------------------------------------------
    subroutine integrate_boundary_scalar_flux(mesh,sdata,elem_info,function_info,iface,ifield,itime,integrand)
        type(mesh_t),           intent(in)      :: mesh
        type(solverdata_t),     intent(inout)   :: sdata
        type(element_info_t),   intent(in)      :: elem_info
        type(function_info_t),  intent(in)      :: function_info
        integer(ik),            intent(in)      :: iface
        integer(ik),            intent(in)      :: ifield
        integer(ik),            intent(in)      :: itime
        type(AD_D),             intent(inout)   :: integrand(:)
        
        ! Data for applying to self and neighbor
        type(AD_D), allocatable     :: integral(:)
        type(AD_D), allocatable     :: integrand_n(:)
        type(function_info_t)       :: function_n
        type(element_info_t)        :: elem_info_n
        integer(ik)                 :: ineighbor_element_l, ineighbor_face, ineighbor_proc, idiff_n, ierr
        logical                     :: parallel_neighbor, diff_none, diff_interior, diff_exterior


        associate ( idomain_l   => elem_info%idomain_l,     &
                    ielement_l  => elem_info%ielement_l,    &
                    ifcn        => function_info%ifcn,      &
                    idonor      => function_info%idepend,   &
                    idiff       => function_info%idiff )



        ! Neighbor indices
        ineighbor_proc      = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_proc
        ineighbor_element_l = mesh%domain(idomain_l)%faces(ielement_l,iface)%get_neighbor_element_l()

        parallel_neighbor = ( IRANK /= ineighbor_proc )

        ! Store quadrature flux for neighbor integral
        integrand_n = integrand

        ! Integrate and apply once
        associate ( weights  => mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%weights_face(iface),                          &
                    valtrans => transpose(mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%interpolator_face('Value',iface)) )

            ! Multiply each component by quadrature weights. The fluxes have already been multiplied by norm
            integrand = integrand * weights
            integral = matmul(valtrans,integrand)

            call store_boundary_integral_residual(     mesh,sdata,elem_info,function_info,iface,ifield,itime,integral)
            call store_boundary_integral_linearization(mesh,sdata,elem_info,function_info,iface,ifield,itime,integral)

        end associate



        ! Integrate and apply second time if there is a neighbor
        if ( ineighbor_element_l /= NO_INTERIOR_NEIGHBOR ) then
            if ( .not. parallel_neighbor ) then

                elem_info_n = element_info(idomain_g         = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_g,         &
                                           idomain_l         = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_l,         &
                                           ielement_g        = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_g,        &
                                           ielement_l        = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_l,        &
                                           iproc             = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_proc,             &
                                           pelem_ID          = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_pelem_ID,         &
                                           coordinate_system = mesh%domain(idomain_l)%faces(ielement_l,iface)%coordinate_system,          &
                                           eqn_ID            = NO_ID,                                                                     &
                                           nfields           = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_nfields,          &
                                           ntime             = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_ntime,            &
                                           nterms_s          = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_nterms_s,         &
                                           nterms_c          = NO_DATA,                                                                   &
                                           dof_start         = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_dof_start,        &
                                           dof_local_start   = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_dof_local_start,  &
                                           recv_comm         = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_comm,                  &
                                           recv_domain       = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_domain,                &
                                           recv_element      = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_element,               &
                                           recv_dof          = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_dof)


                ! Get neighbor face: we have already decided it is a proc-local, interior face, 
                ! so it must have a local neighbor.
                ineighbor_face = mesh%domain(idomain_l)%faces(ielement_l,iface)%get_neighbor_face()

                    
                ! Get linearization block for the neighbor element
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
                function_n%dtype   = function_info%dtype
                function_n%idiff   = idiff_n


                associate ( weights_n    => mesh%domain(idomain_l)%faces(ineighbor_element_l,ineighbor_face)%basis_s%weights_face(ineighbor_face),               &
                            valtrans_n   => transpose(mesh%domain(idomain_l)%faces(ineighbor_element_l,ineighbor_face)%basis_s%interpolator_face('Value',ineighbor_face)) )

                    integrand_n = integrand_n * weights_n

                    ! Integrate and negate for contribution to neighbor element
                    integral = -matmul(valtrans_n,integrand_n)

                    call store_boundary_integral_residual(     mesh,sdata,elem_info_n,function_n,ineighbor_face,ifield,itime,integral)
                    call store_boundary_integral_linearization(mesh,sdata,elem_info_n,function_n,ineighbor_face,ifield,itime,integral)

                end associate

            end if ! .not. parallel_neighbor
        end if ! ielem_n


        end associate

    end subroutine integrate_boundary_scalar_flux
    !***************************************************************************************














    
!    !>  Compute the boundary integral of a flux scalar using the gradient of the test functions 
!    !!  The gradients are wrt to the physical coordiantes
!    !!  BR2x diffusive scheme
!    !!
!    !!      - Adds value contribution to the rhs vector
!    !!      - Adds the derivative contribution to the linearization matrix
!    !!
!    !!  @author Matteo Ugolotti
!    !!
!    !!  Added specific operations for dX differentiation. The goal here is to avoid to carry out a matmul
!    !!  operation involving two AD factors when it is not necessary (fcn_info%dtype == 'dX' or 'none').
!    !!  In this way we avoid to slow down the matmul operation for standard linearization.
!    !!
!    !!  @author Matteo Ugolotti
!    !!  @date   9/4/2018
!    !!
!    !!  Added BC linearization 
!    !!
!    !!  @author Matteo Ugolotti
!    !!  @date   11/26/2018
!    !!
!    !!  Added distance field linearization 
!    !!
!    !!  @author Matteo Ugolotti
!    !!  @date   5/11/2019
!    !!
!    !--------------------------------------------------------------------------------------------------------
!    subroutine integrate_boundary_scalar_flux_grad(mesh,sdata,face_info,function_info,ieqn,itime,integrand_m,integrand_p)
!        type(mesh_t),           intent(in)      :: mesh
!        type(solverdata_t),     intent(inout)   :: sdata
!        type(face_info_t),      intent(in)      :: face_info
!        type(function_info_t),  intent(in)      :: function_info
!        integer(ik),            intent(in)      :: ieqn
!        integer(ik),            intent(in)      :: itime
!        type(AD_D),             intent(inout)   :: integrand_m(:)
!        type(AD_D),             intent(inout)   :: integrand_p(:)
!        
!        ! Data for applying to self and neighbor
!        type(AD_D), allocatable     :: integral(:), integral_1(:), integral_2(:), integral_3(:),            &
!                                       integral_n(:), integral_1_n(:), integral_2_n(:), integral_3_n(:)
!        type(AD_D), allocatable     :: integrand(:), integrand_1(:), integrand_2(:), integrand_3(:),        &
!                                       integrand_n(:), integrand_1_n(:), integrand_2_n(:), integrand_3_n(:),&
!                                       diff_interp_1(:,:), diff_interp_2(:,:), diff_interp_3(:,:)
!        type(face_info_t)           :: face_n
!        type(function_info_t)       :: function_n
!        integer(ik)                 :: ineighbor_element_l, ineighbor_face, ineighbor_proc, idiff_n, ierr,  &
!                                       ineighbor_domain_l, face_type
!        logical                     :: parallel_neighbor, diff_none, diff_interior, diff_exterior,          &
!                                       local_neighbor
!
!
!        associate ( idomain_l   => face_info%idomain_l,     &
!                    ielement_l  => face_info%ielement_l,    &
!                    iface       => face_info%iface,         &
!                    ifcn        => function_info%ifcn,      &
!                    idonor      => function_info%idepend,   &
!                    idiff       => function_info%idiff )
!
!        
!        ! Retrieve face type (BOUNDAY,INTERIOR,CHIMERA)
!        face_type = mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype
!
!        ! Neighbor indices
!        ineighbor_proc      = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_proc
!        ineighbor_element_l = mesh%domain(idomain_l)%faces(ielement_l,iface)%get_neighbor_element_l()
!        ineighbor_face      = mesh%domain(idomain_l)%faces(ielement_l,iface)%get_neighbor_face()
!
!        parallel_neighbor = ( IRANK /= ineighbor_proc )
!        local_neighbor    = ( IRANK == ineighbor_proc )
!
!
!        !
!        ! Integrate and apply once
!        !
!        
!        !
!        ! Retrieve vectors for interior face
!        !
!        associate ( weights         => mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%weights(iface),   &
!                    val_grad1_trans => transpose(mesh%domain(idomain_l)%faces(ielement_l,iface)%grad1),         &
!                    val_grad2_trans => transpose(mesh%domain(idomain_l)%faces(ielement_l,iface)%grad2),         &
!                    val_grad3_trans => transpose(mesh%domain(idomain_l)%faces(ielement_l,iface)%grad3),         &
!                    norm            => mesh%domain(idomain_l)%faces(ielement_l,iface)%norm )
!
!            !
!            ! if dQ, dBC or none linearization
!            !
!            if (function_info%dtype == dQ_DIFF .or. function_info%dtype == dBC_DIFF .or. function_info%dtype == NO_DIFF .or. function_info%dtype == dD_DIFF) then
!
!                ! Multiply each component by quadrature weights.
!                integrand   = integrand_m * weights
!                integrand_1 = integrand * norm(:,1)
!                integrand_2 = integrand * norm(:,2)
!                integrand_3 = integrand * norm(:,3)
!
!                ! Compute partial integrals (Nathan ok)
!                integral_1 = matmul(val_grad1_trans,integrand_1)
!                integral_2 = matmul(val_grad2_trans,integrand_2)
!                integral_3 = matmul(val_grad3_trans,integrand_3)
!
!            !
!            ! if dX linearization
!            !
!            else if (function_info%dtype == dX_DIFF) then
!
!                ! Multiply each component by quadrature weights.
!                integrand   = integrand_m * weights
!                integrand_1 = integrand * differentiate_normal(mesh,face_info,function_info,1) 
!                integrand_2 = integrand * differentiate_normal(mesh,face_info,function_info,2) 
!                integrand_3 = integrand * differentiate_normal(mesh,face_info,function_info,3) 
!
!                ! Compute partial integrals (Nathan ok)
!                diff_interp_1 = transpose(differentiate_face_interior_interpolator('grad1',mesh,face_info,face_info,function_info))
!                diff_interp_2 = transpose(differentiate_face_interior_interpolator('grad2',mesh,face_info,face_info,function_info))
!                diff_interp_3 = transpose(differentiate_face_interior_interpolator('grad3',mesh,face_info,face_info,function_info))
!                integral_1 = matmul(diff_interp_1,integrand_1)
!                integral_2 = matmul(diff_interp_2,integrand_2)
!                integral_3 = matmul(diff_interp_3,integrand_3)
!                !integral_1 = matmul(transpose(differentiate_face_interior_interpolator('grad1',mesh,face_info,face_info,function_info)),integrand_1)
!                !integral_2 = matmul(transpose(differentiate_face_interior_interpolator('grad2',mesh,face_info,face_info,function_info)),integrand_2)
!                !integral_3 = matmul(transpose(differentiate_face_interior_interpolator('grad3',mesh,face_info,face_info,function_info)),integrand_3)
!
!            else
!
!                call chidg_signal_one(FATAL,"integrate_boundary_scalar_flux_grad: unexpected differentiate type",function_info%dtype)
!            
!            end if
!
!
!            integral = (integral_1 + integral_2 + integral_3)
!
!            call store_boundary_integral_residual(     mesh,sdata,face_info,function_info,ieqn,itime,integral)
!            call store_boundary_integral_linearization(mesh,sdata,face_info,function_info,ieqn,itime,integral)
!
!
!        end associate
!
!
!
!
!        !
!        ! Integrate and apply second time if there is a neighbor
!        !
!        ! Only for interior face (no chimera)
!        ! 
!        ! In case of dX linearization, there are two options:
!        !       1) if ( ineighbor_element_l /= NO_INTERIOR_NEIGHBOR .and. function_info%dtype /= dX_DIFF ) then
!        !       2) if ( ineighbor_element_l /= NO_INTERIOR_NEIGHBOR ) then
!        !
!        ! Both return the correct sensitivities but the derivatives of the flux will be different.
!        !
!        if ( ineighbor_element_l /= NO_INTERIOR_NEIGHBOR ) then
!            if ( .not. parallel_neighbor ) then
!
!                face_n%idomain_g  = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_g
!                face_n%idomain_l  = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_l
!                face_n%ielement_g = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_g
!                face_n%ielement_l = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_l
!                face_n%iface      = mesh%domain(idomain_l)%faces(ielement_l,iface)%get_neighbor_face()
!
!                !
!                ! Get linearization block for the neighbor element
!                !
!                diff_none     = (idiff == 0)
!                diff_interior = (idiff == DIAG)
!                diff_exterior = ( (idiff == 1) .or. (idiff == 2) .or. &
!                                  (idiff == 3) .or. (idiff == 4) .or. &
!                                  (idiff == 5) .or. (idiff == 6) )
!                if ( diff_exterior ) then
!                        idiff_n = DIAG
!                else if ( diff_interior ) then
!                        idiff_n = ineighbor_face
!                else if ( diff_none ) then
!                        idiff_n = 0
!                else
!                    call chidg_signal(FATAL,"store_boundary_integrals: unexpected value")
!                end if
!
!
!                function_n%type    = function_info%type
!                function_n%ifcn    = function_info%ifcn
!                function_n%idepend = function_info%idepend
!                function_n%seed    = function_info%seed
!                function_n%idiff   = idiff_n
!                function_n%dtype   = function_info%dtype
!
!                ! The neighbor is local, retrive the valtrans to compute the jump
!                associate(weights_n         =>  mesh%domain(idomain_l)%faces(ineighbor_element_l,ineighbor_face)%basis_s%weights(ineighbor_face),   &
!                          val_grad1_trans_n => transpose(mesh%domain(idomain_l)%faces(ineighbor_element_l,ineighbor_face)%grad1),                   &
!                          val_grad2_trans_n => transpose(mesh%domain(idomain_l)%faces(ineighbor_element_l,ineighbor_face)%grad2),                   &
!                          val_grad3_trans_n => transpose(mesh%domain(idomain_l)%faces(ineighbor_element_l,ineighbor_face)%grad3),                   &
!                          norm_n            => mesh%domain(idomain_l)%faces(ineighbor_element_l,ineighbor_face)%norm )
!                
!                        
!                    ! Multipy by weigths
!                    integrand_n = integrand_p * weights_n
!                
!                    !
!                    ! if dQ, dBC or none linearization
!                    !
!                    if (function_info%dtype == dQ_DIFF .or. function_info%dtype == dBC_DIFF .or. function_info%dtype == NO_DIFF .or. function_info%dtype == dD_DIFF) then
!                    
!                        
!                        ! Compute three components of the integrand
!                        integrand_1_n = integrand_n * norm_n(:,1)
!                        integrand_2_n = integrand_n * norm_n(:,2)
!                        integrand_3_n = integrand_n * norm_n(:,3)
!                
!                        ! Compute partial integrals (Nathan ok)
!                        integral_1_n = matmul(val_grad1_trans_n,integrand_1_n)
!                        integral_2_n = matmul(val_grad2_trans_n,integrand_2_n)
!                        integral_3_n = matmul(val_grad3_trans_n,integrand_3_n)
!
!
!                    !
!                    ! if dX linearization
!                    !
!                    else if (function_info%dtype == dX_DIFF) then
!
!                        ! Multiply each component by normal of neighbor face.
!                        integrand_1_n = integrand_n * differentiate_normal(mesh,face_n,function_n,1) 
!                        integrand_2_n = integrand_n * differentiate_normal(mesh,face_n,function_n,2) 
!                        integrand_3_n = integrand_n * differentiate_normal(mesh,face_n,function_n,3) 
!
!                        ! Compute partial integrals (Nathan ok)
!                        diff_interp_1 = transpose(differentiate_face_interior_interpolator('grad1',mesh,face_n,face_n,function_n))
!                        diff_interp_2 = transpose(differentiate_face_interior_interpolator('grad2',mesh,face_n,face_n,function_n))
!                        diff_interp_3 = transpose(differentiate_face_interior_interpolator('grad3',mesh,face_n,face_n,function_n))
!                        integral_1_n = matmul(diff_interp_1,integrand_1_n)
!                        integral_2_n = matmul(diff_interp_2,integrand_2_n)
!                        integral_3_n = matmul(diff_interp_3,integrand_3_n)
!
!                    else
!
!                        call chidg_signal_one(FATAL,"integrate_boundary_scalar_flux_grad: unexpected differentiate type",function_info%dtype)
!                    
!                    end if
!
!                        integral_n = (integral_1_n + integral_2_n + integral_3_n)
!
!                        call store_boundary_integral_residual(     mesh,sdata,face_n,function_n,ieqn,itime,integral_n)
!                        call store_boundary_integral_linearization(mesh,sdata,face_n,function_n,ieqn,itime,integral_n)
!                
!                end associate
!
!            end if ! .not. parallel_neighbor
!        end if ! ielem_n
!
!
!        end associate
!
!    end subroutine integrate_boundary_scalar_flux_grad
!    !********************************************************************************************************
!



















    !> Store volume integral values to RHS vector, and partial derivatives to LIN block matrix
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!  @param[in]      integral    Array of autodiff values containing integrals and partial derivatives for the RHS vector and LIN linearization matrix
    !!  @param[inout]   rhs         Right-hand side vector
    !!  @param[inout]   lin         Block matrix storing the linearization of the spatial scheme
    !!  @param[in]      ielem       Element index for applying to the correct location in RHS and LIN
    !!  @param[in]      ifield      Variable index
    !!  @param[in]      idiff       Block index for the correct linearization block for the current element
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !----------------------------------------------------------------------------------------
    subroutine store_volume_integrals(mesh,sdata,elem_info,fcn_info,ifield,itime,integral)
        type(mesh_t),           intent(in)      :: mesh
        type(solverdata_t),     intent(inout)   :: sdata
        type(element_info_t),   intent(in)      :: elem_info
        type(function_info_t),  intent(in)      :: fcn_info
        integer(ik),            intent(in)      :: ifield
        integer(ik),            intent(in)      :: itime
        type(AD_D),             intent(inout)   :: integral(:)

        integer(ik)         :: i
        logical             :: conforming_face, boundary_face, chimera_face, &
                               diff_none, diff_interior, diff_exterior,      &
                               diff_dq, diff_dx, diff_dbc, diff_dd
        real(rk)            :: vals(size(integral)), integral_dbc(size(integral))

        associate ( idom  => elem_info%idomain_l,   &
                    ielem => elem_info%ielement_l,  &
                    idiff => fcn_info%idiff )

        ! Select type of linearization
        diff_dq  = ( fcn_info%dtype == dQ_DIFF .or. fcn_info%dtype == NO_DIFF)
        diff_dx  = ( fcn_info%dtype == dX_DIFF)
        diff_dbc = ( fcn_info%dtype == dBC_DIFF)
        diff_dd  = ( fcn_info%dtype == dD_DIFF)

        ! Select geometric entity
        diff_none     = (idiff == 0)
        diff_interior = (idiff == DIAG)
        diff_exterior = ( (idiff == 1) .or. (idiff == 2) .or. &
                          (idiff == 3) .or. (idiff == 4) .or. &
                          (idiff == 5) .or. (idiff == 6) )


        ! Only store rhs once. 
        !   - If we are differentiating things, only apply function once, when wrt interior
        !       idiff == DIAG
        !   - If we are not differentiating things, just apply function.
        !       idiff == 0
        if ( (diff_interior .or. diff_none) .and. diff_dq ) then
            vals = -integral(:)%x_ad_
            call sdata%rhs%add_field(vals,elem_info,ifield,itime)
        end if

        ! Negate derivatives before adding to linearization
        do i = 1,size(integral)
            integral(i)%xp_ad_ = -integral(i)%xp_ad_
        end do


        ! Check if linearization is with respect to an exterior element. 
        ! Only need this for diff_exterior, and idiff is undefined as a face index for 
        ! diff_interior or diff_none.
        if ( diff_exterior ) then
            conforming_face = (mesh%domain(idom)%faces(ielem,idiff)%ftype == INTERIOR)
            boundary_face   = (mesh%domain(idom)%faces(ielem,idiff)%ftype == BOUNDARY)
            chimera_face    = (mesh%domain(idom)%faces(ielem,idiff)%ftype == CHIMERA )
        end if


        ! Store linearization
        if (diff_dq) then

            if ( diff_interior ) then
                call sdata%lhs%store(integral,elem_info,fcn_info%seed,ifield,itime)
            else if ( diff_exterior .and. conforming_face ) then
                call sdata%lhs%store(integral,elem_info,fcn_info%seed,ifield,itime)
            else if ( diff_exterior .and. chimera_face    ) then
                call sdata%lhs%store_chimera(integral,elem_info,fcn_info%seed,ifield,itime)
            else if ( diff_exterior .and. boundary_face   ) then
                call sdata%lhs%store_bc(integral,elem_info,fcn_info%seed,ifield,itime)
            else if ( diff_none ) then
                ! No derivatives to store
            else
                call chidg_signal(FATAL,"store_volume_integrals: Invalid condition for storing integrals. Could be a bad face type or linearization direction")
            end if

        else if (diff_dx) then

            if ( diff_interior ) then
                call sdata%adjointx%Rx%store(integral,elem_info,fcn_info%seed,ifield,itime)
            else if ( diff_exterior .and. conforming_face ) then
                call sdata%adjointx%Rx%store(integral,elem_info,fcn_info%seed,ifield,itime)
            else if ( diff_exterior .and. chimera_face    ) then
                call sdata%adjointx%Rx%store_chimera(integral,elem_info,fcn_info%seed,ifield,itime)
            else if ( diff_exterior .and. boundary_face   ) then
                call sdata%adjointx%Rx%store_bc(integral,elem_info,fcn_info%seed,ifield,itime)
            else if ( diff_none ) then
                ! No derivatives to store
            else
                call chidg_signal(FATAL,"store_volume_integrals: Invalid condition for storing integrals. Could be a bad face type or linearization direction")
            end if

        else if (diff_dd) then
            if ( diff_interior ) then
                call sdata%adjoint%Rd(1)%store(integral,elem_info,fcn_info%seed,ifield,itime)
            else if ( diff_exterior .and. conforming_face ) then
                call sdata%adjoint%Rd(1)%store(integral,elem_info,fcn_info%seed,ifield,itime)
            else if ( diff_exterior .and. chimera_face    ) then
                call sdata%adjoint%Rd(1)%store_chimera(integral,elem_info,fcn_info%seed,ifield,itime)
            else if ( diff_exterior .and. boundary_face   ) then
                call sdata%adjoint%Rd(1)%store_bc(integral,elem_info,fcn_info%seed,ifield,itime)
            else if ( diff_none ) then
                ! No derivatives to store
            else
                call chidg_signal(FATAL,"store_volume_integrals: Invalid condition for storing integrals. Could be a bad face type or linearization direction")
            end if

        else if (diff_dbc) then

            do i = 1,size(integral)
                 integral_dbc(i) = integral(i)%xp_ad_(1)
            end do
            !vals = sdata%adjointbc%Ra%dom(idom)%vecs(ielem)%getvar(ieqn,itime) + integral_dbc
            !call sdata%adjointbc%Ra%dom(idom)%vecs(ielem)%setvar(ieqn,itime,vals)
            vals = sdata%adjointbc%Ra%get_field(elem_info,ifield,itime) + integral_dbc
            call sdata%adjointbc%Ra%set_field(vals,elem_info,ifield,itime)

        else
            call chidg_signal(FATAL,"store_volume_integrals: type of linearization not defined. Possible linearization dQ, dX, dBC or none.")
        end if





        end associate

    end subroutine store_volume_integrals
    !****************************************************************************************










    !> Store boundary integral values to RHS vector, and partial derivatives to LIN block matrix
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!  @param[in]      integral    Array of autodiff values containing integrals and partial derivatives for the RHS vector and LIN linearization matrix
    !!  @param[inout]   rhs         Right-hand side vector
    !!  @param[inout]   lin         Block matrix storing the linearization of the spatial scheme
    !!  @param[in]      ielem       Element index for applying to the correct location in RHS and LIN
    !!  @param[in]      ifield      Variable index
    !!  @param[in]      iblk        Block index for the correct linearization block for the current element
    !!
    !!  @author Mayank Sharma + matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine store_boundary_integral_residual(mesh,sdata,elem_info,function_info,iface,ifield,itime,integral)
        type(mesh_t),           intent(in)      :: mesh
        type(solverdata_t),     intent(inout)   :: sdata
        type(element_info_t),   intent(in)      :: elem_info
        type(function_info_t),  intent(in)      :: function_info
        integer(ik),            intent(in)      :: iface
        integer(ik),            intent(in)      :: ifield
        integer(ik),            intent(in)      :: itime
        type(AD_D),             intent(inout)   :: integral(:)

        real(rk)        :: vals(size(integral))
        integer(ik)     :: nterms_s,nfields

        logical         :: boundary_face, chimera_face, conforming_face, diff_interior, diff_none, diff_exterior
        logical         :: add_flux = .false.

        if (function_info%dtype == dQ_DIFF .or. function_info%dtype == NO_DIFF) then

            associate ( idomain_l  => elem_info%idomain_l, ielement_l  => elem_info%ielement_l,     &
                        ifcn  => function_info%ifcn,       idepend => function_info%idepend,     idiff => function_info%idiff )

                conforming_face = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR)
                boundary_face   = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == BOUNDARY)
                chimera_face    = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA )

                diff_none     = (idiff == 0)
                diff_interior = (idiff == DIAG)
                diff_exterior = ( (idiff == 1) .or. (idiff == 2) .or. &
                                  (idiff == 3) .or. (idiff == 4) .or. &
                                  (idiff == 5) .or. (idiff == 6) )

                nterms_s = mesh%domain(idomain_l)%faces(ielement_l,iface)%nterms_s
                nfields  = mesh%domain(idomain_l)%faces(ielement_l,iface)%nfields

                 
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
                if ( ((boundary_face .and. diff_exterior) .or. &
                      (boundary_face .and. diff_none    )) .and. &
                      function_info%seed%itime == itime ) then

                    if (idepend == 1) then
                        vals = integral(:)%x_ad_
                        call sdata%rhs%add_field(vals,elem_info,ifield,itime)
                    end if


                else if ( (chimera_face .and. diff_interior) .or. &
                          (chimera_face .and. diff_none    ) ) then

                    if (idepend == 1) then
                        vals = integral(:)%x_ad_
                        call sdata%rhs%add_field(vals,elem_info,ifield,itime)
                    end if


                else if ( conforming_face ) then

                    ! Check if particular flux function has been added already
                    add_flux = sdata%function_status%compute_function_equation( elem_info,function_info,iface,ifield)

                    ! Store if needed
                    if ( add_flux ) then
                        ! Add to residual and store
                        vals = integral(:)%x_ad_
                        call sdata%rhs%add_field(vals,elem_info,ifield,itime)

                        ! Register flux was stored
                        call sdata%function_status%register_function_computed(elem_info,function_info,iface,ifield)
                    end if

                end if

            end associate

        end if !function_info%dtype /= 'dX'

    end subroutine store_boundary_integral_residual
    !******************************************************************************************









    !> Store boundary integral values to RHS vector, and partial derivatives to LIN block matrix
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!  @param[in]      integral    Array of autodiff values containing integrals and partial derivatives for the RHS vector and LIN linearization matrix
    !!  @param[inout]   rhs         Right-hand side vector
    !!  @param[inout]   lin         Block matrix storing the linearization of the spatial scheme
    !!  @param[in]      ielem       Element index for applying to the correct location in RHS and LIN
    !!  @param[in]      ifield      Variable index
    !!  @param[in]      idiff       Block index for the correct linearization block for the current element
    !!
    !---------------------------------------------------------------------------------------
    subroutine store_boundary_integral_linearization(mesh,sdata,elem_info,function_info,iface,ifield,itime,integral)
        type(mesh_t),           intent(in)      :: mesh
        type(solverdata_t),     intent(inout)   :: sdata
        type(element_info_t),   intent(in)      :: elem_info
        type(function_info_t),  intent(in)      :: function_info
        integer(ik),            intent(in)      :: iface
        integer(ik),            intent(in)      :: ifield
        integer(ik),            intent(in)      :: itime
        type(AD_D),             intent(inout)   :: integral(:)

        integer(ik) :: i, idomain_l, ielement_l, ChiID
        integer(ik) :: idiff, ifcn
        real(rk)    :: vals(size(integral)), integral_bc(size(integral))

        logical :: conforming_face, boundary_face, chimera_face, &
                   diff_none, diff_interior, diff_exterior, add_linearization, &
                   diff_dq, diff_dx, diff_dbc, diff_dd

        associate ( idomain_l => elem_info%idomain_l, ielement_l => elem_info%ielement_l,  &
                    ifcn      => function_info%ifcn,  idepend    => function_info%idepend, idiff => function_info%idiff, seed => function_info%seed )

        ! Define the type of linearization and store location based on this
        diff_dq  = ( function_info%dtype == dQ_DIFF .or. function_info%dtype == NO_DIFF)
        diff_dx  = ( function_info%dtype == dX_DIFF)
        diff_dbc = ( function_info%dtype == dBC_DIFF)
        diff_dd  = ( function_info%dtype == dD_DIFF)

        
        conforming_face = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR)
        boundary_face   = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == BOUNDARY)
        chimera_face    = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA )

        diff_none     = ( idiff == 0 )
        diff_interior = ( idiff == DIAG )
        diff_exterior = ( (idiff == 1) .or. (idiff == 2) .or. &
                          (idiff == 3) .or. (idiff == 4) .or. &
                          (idiff == 5) .or. (idiff == 6) )

        associate ( lhs => sdata%lhs,           &
                    !Rd => sdata%adjoint%Rd(1),  &
                    Rx => sdata%adjointx%Rx,    &
                    Ra => sdata%adjointbc%Ra )

        if (diff_interior .or. diff_exterior) then

            ! Store linearization. Rules for different face types.
            if ( chimera_face ) then

                if (diff_exterior) then
                    ! Store linearization of Chimera boundary donor elements.
                    if (diff_dq) call lhs%store_chimera(integral,elem_info,seed,ifield,itime)
                    if (diff_dx) call Rx%store_chimera(integral,elem_info,seed,ifield,itime)
                    if (diff_dd) call sdata%adjoint%Rd(1)%store_chimera(integral,elem_info,seed,ifield,itime)
                else
                    ! Store linearization of Chimera boundary receiver element. Since this could be computed multiple times,
                    ! we just store it once.
                    if (idepend == 1) then
                        if (diff_dq) call lhs%store(integral,elem_info,seed,ifield,itime)
                        if (diff_dx) call Rx%store(integral,elem_info,seed,ifield,itime)
                        if (diff_dd) call sdata%adjoint%Rd(1)%store(integral,elem_info,seed,ifield,itime)
                    end if
                end if


            else if ( boundary_face ) then
                if (diff_dq) call lhs%store_bc(integral,elem_info,seed,ifield,itime)
                if (diff_dx) call Rx%store_bc(integral,elem_info,seed,ifield,itime)
                if (diff_dd) call sdata%adjoint%Rd(1)%store_bc(integral,elem_info,seed,ifield,itime)


            else if ( conforming_face ) then

                add_linearization = sdata%function_status%linearize_function_equation(elem_info,function_info,iface,ifield)

                ! Store linearization if not already stored
                if ( add_linearization ) then
                    if (diff_dq) call lhs%store(integral,elem_info,seed,ifield,itime)
                    if (diff_dx) call Rx%store(integral,elem_info,seed,ifield,itime)
                    if (diff_dd) call sdata%adjoint%Rd(1)%store(integral,elem_info,seed,ifield,itime)
                    ! Register flux as linearized
                    call sdata%function_status%register_function_linearized(elem_info,function_info,iface,ifield)
                end if

            end if ! ftype


            ! Specialized storage for BC linearization
            ! Here, we do not really care if it is chimera, conforming or boundary face.
            ! TODO: this could be moved under 'boundary_face' if-statement, since BC derivatives
            !       comes only from the boundary faces.
            if (diff_dbc) then
                do i = 1,size(integral)
                     integral_bc(i) = integral(i)%xp_ad_(1)
                end do
                !vals = Ra%dom(idomain_l)%vecs(ielement_l)%getvar(ifield,itime) + integral_bc
                !call Ra%dom(idomain_l)%vecs(ielement_l)%setvar(ifield,itime,vals)
                call Ra%add_field(integral_bc,elem_info,ifield,itime)
                
            end if

        end if ! diff_...
        end associate

        end associate 

    end subroutine store_boundary_integral_linearization
    !******************************************************************************************




!    !> Store boundary integral values to RHS vector, and partial derivatives to LIN block matrix
!    !!
!    !!  @author Nathan A. Wukie
!    !!
!    !!
!    !!  @param[in]      integral    Array of autodiff values containing integrals and partial derivatives for the RHS vector and LIN linearization matrix
!    !!  @param[inout]   rhs         Right-hand side vector
!    !!  @param[inout]   lin         Block matrix storing the linearization of the spatial scheme
!    !!  @param[in]      ielem       Element index for applying to the correct location in RHS and LIN
!    !!  @param[in]      ieqn        Variable index
!    !!  @param[in]      idiff       Block index for the correct linearization block for the current element
!    !!
!    !---------------------------------------------------------------------------------------
!    subroutine store_boundary_integral_linearization(mesh,sdata,elem_info,function_info,iface,ieqn,itime,integral)
!        type(mesh_t),           intent(in)      :: mesh
!        type(solverdata_t),     intent(inout)   :: sdata
!        type(element_info_t),   intent(in)      :: elem_info
!        type(function_info_t),  intent(in)      :: function_info
!        integer(ik),            intent(in)      :: iface
!        integer(ik),            intent(in)      :: ieqn
!        integer(ik),            intent(in)      :: itime
!        type(AD_D),             intent(inout)   :: integral(:)
!
!        integer(ik)                 :: i, idomain_l, ielement_l, ChiID
!        integer(ik)                 :: idiff, ifcn
!        real(rk)                    :: vals(size(integral))
!
!        logical :: conforming_face, boundary_face, chimera_face, &
!                   diff_none, diff_interior, diff_exterior, add_linearization
!
!        associate ( idomain_l => elem_info%idomain_l, ielement_l => elem_info%ielement_l,  &
!                    ifcn      => function_info%ifcn,  idepend    => function_info%idepend, idiff => function_info%idiff, seed => function_info%seed )
!
!        
!        conforming_face = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR)
!        boundary_face   = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == BOUNDARY)
!        chimera_face    = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA )
!
!        diff_none     = ( idiff == 0 )
!        diff_interior = ( idiff == DIAG )
!        diff_exterior = ( (idiff == 1) .or. (idiff == 2) .or. &
!                          (idiff == 3) .or. (idiff == 4) .or. &
!                          (idiff == 5) .or. (idiff == 6) )
!
!        associate ( lhs => sdata%lhs)
!
!        if (diff_interior .or. diff_exterior) then
!
!            ! Store linearization. Rules for different face types.
!            if ( chimera_face ) then
!
!                if (diff_exterior) then
!                    ! Store linearization of Chimera boundary donor elements.
!                    call lhs%store_chimera(integral,elem_info,seed,ieqn,itime)
!                else
!                    ! Store linearization of Chimera boundary receiver element. Since this could be computed multiple times,
!                    ! we just store it once.
!                    if (idepend == 1) then
!                        call lhs%store(integral,elem_info,seed,ieqn,itime)
!                    end if
!                end if
!
!
!
!            else if ( boundary_face ) then
!                call lhs%store_bc(integral,elem_info,seed,ieqn,itime)
!
!
!
!            else if ( conforming_face ) then
!
!
!                add_linearization = sdata%function_status%linearize_function_equation( elem_info, function_info, iface, ieqn)
!
!                ! Store linearization if not already stored
!                if ( add_linearization ) then
!                    ! Store linearization
!                    call lhs%store(integral,elem_info,seed,ieqn,itime)
!
!                    ! Register flux as linearized
!                    call sdata%function_status%register_function_linearized( elem_info, function_info, iface, ieqn)
!                end if
!
!
!            end if ! ftype
!
!        end if ! diff_...
!        end associate
!
!        end associate 
!
!    end subroutine store_boundary_integral_linearization
!    !******************************************************************************************






end module mod_integrate
