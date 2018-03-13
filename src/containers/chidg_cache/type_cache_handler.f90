module type_cache_handler
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: NFACES, INTERIOR, CHIMERA, BOUNDARY, DIAG, NO_PROC,   &
                                  ME, NEIGHBOR, HALF, ONE,                              &
                                  XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, NO_ID
    use mod_DNAD_tools,     only: face_compute_seed, element_compute_seed
    use mod_interpolate,    only: interpolate_face_autodiff, interpolate_element_autodiff
    use mod_chidg_mpi,      only: IRANK
    use DNAD_D

    use type_chidg_cache,       only: chidg_cache_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_equation_set,      only: equation_set_t
    use type_bc_state_group,    only: bc_state_group_t

    implicit none



    !>  An object for handling cache operations. Particularly, updating the cache contents.
    !!
    !!  The problem solved here is this. The cache is used in operator_t's to pull data
    !!  computed at quadrature nodes. The cache also requires bc_operators's to precompute
    !!  the boundary condition solution as an external state and also to compute the BR2
    !!  diffusion lifting operators. This introduced a pesky circular dependency.
    !!
    !!  The problem was solved by introducing this cache_handler object. This separates the
    !!  cache behavior from the cache storage. The operator_t's need the cache storage. 
    !!  They don't need to know how the data got there.
    !!
    !!  So this higher-level interface sits outside of the hierarchy that caused the circular
    !!  dependency to handle the cache behavior, such as how it gets updated.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/7/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    type, public :: cache_handler_t



    contains


        procedure   :: update   ! Resize/Update the cache fields


        procedure, private  :: update_auxiliary_fields
        procedure, private  :: update_primary_fields


        procedure, private  :: update_auxiliary_interior
        procedure, private  :: update_auxiliary_exterior
        procedure, private  :: update_auxiliary_element
        procedure, private  :: update_auxiliary_bc

        procedure, private  :: update_primary_interior
        procedure, private  :: update_primary_exterior
        procedure, private  :: update_primary_element
        procedure, private  :: update_primary_bc
        procedure, private  :: update_primary_lift

        procedure, private  :: update_model_interior
        procedure, private  :: update_model_exterior
        procedure, private  :: update_model_element
        procedure, private  :: update_model_bc

        procedure, private  :: update_lift_faces_internal
        procedure, private  :: update_lift_faces_external

    end type cache_handler_t
    !****************************************************************************************





contains


    !>  Resize chidg_cache in worker, update cache components.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/7/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update(self,worker,equation_set,bc_state_group,components,face,differentiate,lift)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),     intent(inout)   :: bc_state_group(:)
        character(*),               intent(in)      :: components
        integer(ik),                intent(in)      :: face
        logical,                    intent(in)      :: differentiate
        logical,                    intent(in)      :: lift

        integer(ik) :: idomain_l, ielement_l, iface, eqn_ID, face_min, face_max
        logical     :: compute_gradients, valid_indices, update_interior_faces, update_exterior_faces, update_element

        type(AD_D), allocatable, dimension(:) :: grad1_mom3, grad2_mom3, grad3_mom3

        !
        ! Check for valid indices
        !
        valid_indices = (worker%element_info%idomain_l /= 0) .and. &
                        (worker%element_info%ielement_l /= 0) .and. &
                        (worker%itime /= 0)

        if (.not. valid_indices) call chidg_signal(FATAL,"cache_handler%update: Bad domain/element/time indices were detected during update.")



        !
        ! Store lift indicator in worker
        !
        worker%contains_lift = lift


        !
        ! Check for valid components
        !
        select case(trim(components))
            case('all')
                update_interior_faces = .true.
                update_exterior_faces = .true.
                update_element        = .true.
            case('element')
                update_interior_faces = .false.
                update_exterior_faces = .false.
                update_element        = .true.
            case('faces')
                update_interior_faces = .true.
                update_exterior_faces = .true.
                update_element        = .false.
            case('interior faces')
                update_interior_faces = .true.
                update_exterior_faces = .false.
                update_element        = .false.
            case('exterior faces')
                update_interior_faces = .false.
                update_exterior_faces = .true.
                update_element        = .false.
            case default
                call chidg_signal_one(FATAL,"cache_handler%update: Bad 'components' argument.",trim(components))
        end select



        !
        ! Set range of faces to update
        !
        if (face == NO_ID) then
            face_min = 1        ! Update all faces
            face_max = NFACES
        else
            face_min = face     ! Only update one face
            face_max = face
        end if



        !
        ! Resize cache
        !
        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        call worker%cache%resize(worker%mesh,worker%prop,idomain_l,ielement_l,differentiate,lift)



        !
        ! Determine if we want to update gradient terms in the cache
        !
        eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
        compute_gradients = (allocated(equation_set(eqn_ID)%volume_diffusive_operator)   .or. &
                             allocated(equation_set(eqn_ID)%boundary_diffusive_operator) )
!        compute_gradients = .true.




        !
        ! Update fields
        !
        call self%update_auxiliary_fields(worker,equation_set,bc_state_group,differentiate)
        call self%update_primary_fields(  worker,equation_set,bc_state_group,differentiate,compute_gradients,update_element, update_interior_faces, update_exterior_faces, face_min, face_max)

        if (update_element) call self%update_model_element(worker,equation_set,bc_state_group,differentiate,model_type='f(Q-)')



        !
        ! Compute f(Q-) models. Interior, Exterior, BC, Element
        !
        do iface = face_min,face_max

            ! Update worker face index
            call worker%set_face(iface)

            ! Update face interior/exterior/bc states.
            if (update_interior_faces) call self%update_model_interior(worker,equation_set,bc_state_group,differentiate,model_type='f(Q-)')
            if (update_exterior_faces) call self%update_model_exterior(worker,equation_set,bc_state_group,differentiate,model_type='f(Q-)')

        end do !iface



        !
        ! Compute f(Q-) models. Interior, Exterior, BC, Element
        !
        do iface = face_min,face_max

            ! Update worker face index
            call worker%set_face(iface)

            if (update_exterior_faces) call self%update_primary_bc(worker,equation_set,bc_state_group,differentiate)
            if (update_exterior_faces) call self%update_model_bc(  worker,equation_set,bc_state_group,differentiate,model_type='f(Q-)')

            if (update_interior_faces) call self%update_model_interior(worker,equation_set,bc_state_group,differentiate,model_type='f(Q-,Q+)')
            if (update_exterior_faces) call self%update_model_exterior(worker,equation_set,bc_state_group,differentiate,model_type='f(Q-,Q+)')
            if (update_exterior_faces) call self%update_model_bc(      worker,equation_set,bc_state_group,differentiate,model_type='f(Q-,Q+)')


        end do !iface











        if (update_element) call self%update_model_element(worker,equation_set,bc_state_group,differentiate,model_type='f(Q-,Q+)')





        !
        ! Compute f(Q-,Q+), f(Grad(Q) models. Interior, Exterior, BC, Element
        !
        !compute_gradients = .false.
        if (compute_gradients) then

!            !
!            ! Store lift indicator in worker
!            !
!            worker%contains_lift = lift


            !
            ! Update lifting operators for second-order pde's
            !
            if (lift) call self%update_primary_lift(worker,equation_set,bc_state_group,differentiate)


            !
            ! Loop through faces and cache 'internal', 'external' interpolated states
            !
            do iface = face_min,face_max

                ! Update worker face index
                call worker%set_face(iface)

                ! Update face interior/exterior/bc states.
                if (update_interior_faces) call self%update_model_interior(worker,equation_set,bc_state_group,differentiate,model_type='f(Grad(Q))')
                if (update_exterior_faces) call self%update_model_exterior(worker,equation_set,bc_state_group,differentiate,model_type='f(Grad(Q))')
                if (update_exterior_faces) call self%update_model_bc(      worker,equation_set,bc_state_group,differentiate,model_type='f(Grad(Q))')

            end do !iface


            
            !
            ! Update model 'element' cache entries
            !
            if (update_element) call self%update_model_element(worker,equation_set,bc_state_group,differentiate,model_type='f(Grad(Q))')

!        else
!            !
!            ! Store lift indicator in worker
!            !
!            worker%contains_lift = .false.

        end if ! compute_gradients



    end subroutine update
    !****************************************************************************************







    !>  Update the cache entries for the primary fields.
    !!
    !!  Activities:
    !!      #1: Loop through faces, update 'face interior', 'face exterior' caches for 
    !!          'value' and 'gradients'
    !!      #2: Update the 'element' cache for 'value' and 'gradients'
    !!      #3: Update the lifting operators for all cache components
    !!          (interior, exterior, element)
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/7/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_primary_fields(self,worker,equation_set,bc_state_group,differentiate,compute_gradients,update_element, update_interior_faces, update_exterior_faces, face_min, face_max)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),     intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate
        logical,                    intent(in)      :: compute_gradients
        logical,                    intent(in)      :: update_element
        logical,                    intent(in)      :: update_interior_faces
        logical,                    intent(in)      :: update_exterior_faces
        integer(ik),                intent(in)      :: face_min
        integer(ik),                intent(in)      :: face_max

        integer(ik)                                 :: idomain_l, ielement_l, iface, &
                                                       idepend, ieqn, idiff
        character(:),   allocatable                 :: field
        type(AD_D),     allocatable, dimension(:)   :: value_gq, grad1_gq, grad2_gq, grad3_gq


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 


        !
        ! Loop through faces and cache 'internal', 'external' interpolated states
        !
        do iface = face_min,face_max

            ! Update worker face index
            call worker%set_face(iface)


            ! Update face interior/exterior/bc states.
            if (update_interior_faces) call self%update_primary_interior(worker,equation_set,bc_state_group,differentiate,compute_gradients)
            if (update_exterior_faces) call self%update_primary_exterior(worker,equation_set,bc_state_group,differentiate,compute_gradients)


        end do !iface


        !
        ! Update 'element' cache
        !
        if (update_element) call self%update_primary_element(worker,equation_set,bc_state_group,differentiate,compute_gradients)


    end subroutine update_primary_fields
    !****************************************************************************************










    !>  Update the cache entries for the auxiliary fields.
    !!
    !!  Activities:
    !!      #1: Loop through faces, update 'face interior', 'face exterior' caches for 
    !!          'value' and 'gradients'
    !!      #2: Update the 'element' cache for 'value' and 'gradients'
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/7/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_auxiliary_fields(self,worker,equation_set,bc_state_group,differentiate)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),     intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate

        integer(ik)                                 :: idomain_l, ielement_l, iface, idepend, &
                                                       ieqn, ifield, iaux_field, idiff
        character(:),   allocatable                 :: field
        type(AD_D),     allocatable, dimension(:)   :: value_gq, grad1_gq, grad2_gq, grad3_gq


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 


        !
        ! Loop through faces and cache internal, external interpolated states
        !
        do iface = 1,NFACES

            ! Update worker face index
            call worker%set_face(iface)


            ! Update face interior/exterior states.
            call self%update_auxiliary_interior(worker,equation_set,bc_state_group,differentiate)
            call self%update_auxiliary_exterior(worker,equation_set,bc_state_group,differentiate)
            call self%update_auxiliary_bc(      worker,equation_set,bc_state_group,differentiate)


        end do !iface



        !
        ! Update cache 'element' data
        !
        call self%update_auxiliary_element(worker,equation_set,bc_state_group,differentiate)



    end subroutine update_auxiliary_fields
    !****************************************************************************************








    !>  Update the primary field 'element' cache entries.
    !!
    !!  Computes the 'value' and 'gradient' entries.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/9/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_primary_element(self,worker,equation_set,bc_state_group,differentiate,compute_gradients)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),     intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate
        logical,                    intent(in)      :: compute_gradients

        integer(ik)                                 :: idepend, ieqn, idomain_l, ielement_l, &
                                                       iface, idiff, eqn_ID
        character(:),   allocatable                 :: field
        real(rk),       allocatable                 :: ale_Dinv(:,:,:)
        real(rk),       allocatable, dimension(:)   :: ale_g, ale_g_grad1, ale_g_grad2, ale_g_grad3
        type(AD_D),     allocatable, dimension(:)   :: value_u, grad1_u, grad2_u, grad3_u, grad1_tmp, grad2_tmp, grad3_tmp


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface


        !
        ! Element primary fields volume 'value' cache. Only depends on interior element
        !
        if (differentiate) then
            idiff = DIAG
        else
            idiff = 0
        end if


        !
        ! Compute Value/Gradients
        !
        idepend = 1
        eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
        worker%function_info%seed    = element_compute_seed(worker%mesh,idomain_l,ielement_l,idepend,idiff)
        worker%function_info%idepend = idepend
        do ieqn = 1,worker%mesh%domain(idomain_l)%neqns
            field = worker%prop(eqn_ID)%get_primary_field_name(ieqn)

            !
            ! Interpolate U
            !

            ! Interpolate modes to nodes on undeformed element
            value_u = interpolate_element_autodiff(worker%mesh,worker%solverdata%q,worker%element_info,worker%function_info,ieqn,worker%itime,'value')

            ! Get ALE transformation data
            ale_g = worker%get_det_jacobian_grid_element('value')

            ! Compute transformation to deformed element
            value_u = (value_u/ale_g)

            ! Store quantities valid on the deformed element
            call worker%cache%set_data(field,'element',value_u,'value',0,worker%function_info%seed)







            !
            ! Interpolate Grad(U)
            !
            if (compute_gradients) then
                !
                ! Interpolate modes to nodes on undeformed element
                !
                grad1_u = interpolate_element_autodiff(worker%mesh,worker%solverdata%q,worker%element_info,worker%function_info,ieqn,worker%itime,'grad1')
                grad2_u = interpolate_element_autodiff(worker%mesh,worker%solverdata%q,worker%element_info,worker%function_info,ieqn,worker%itime,'grad2')
                grad3_u = interpolate_element_autodiff(worker%mesh,worker%solverdata%q,worker%element_info,worker%function_info,ieqn,worker%itime,'grad3')


                !
                ! Get ALE transformation data
                !
                ale_g_grad1 = worker%get_det_jacobian_grid_element('grad1')
                ale_g_grad2 = worker%get_det_jacobian_grid_element('grad2')
                ale_g_grad3 = worker%get_det_jacobian_grid_element('grad3')
                ale_Dinv    = worker%get_inv_jacobian_grid_element()

                
                !
                ! Compute transformation to deformed element
                !
                grad1_tmp = grad1_u-(value_u)*ale_g_grad1
                grad2_tmp = grad2_u-(value_u)*ale_g_grad2
                grad3_tmp = grad3_u-(value_u)*ale_g_grad3

                grad1_u = (ale_Dinv(1,1,:)*grad1_tmp + ale_Dinv(2,1,:)*grad2_tmp + ale_Dinv(3,1,:)*grad3_tmp)/ale_g
                grad2_u = (ale_Dinv(1,2,:)*grad1_tmp + ale_Dinv(2,2,:)*grad2_tmp + ale_Dinv(3,2,:)*grad3_tmp)/ale_g
                grad3_u = (ale_Dinv(1,3,:)*grad1_tmp + ale_Dinv(2,3,:)*grad2_tmp + ale_Dinv(3,3,:)*grad3_tmp)/ale_g


                !
                ! Store quantities valid on the deformed element
                !
                call worker%cache%set_data(field,'element',grad1_u,'gradient',1,worker%function_info%seed)
                call worker%cache%set_data(field,'element',grad2_u,'gradient',2,worker%function_info%seed)
                call worker%cache%set_data(field,'element',grad3_u,'gradient',3,worker%function_info%seed)


            end if

        end do !ieqn



    end subroutine update_primary_element
    !*****************************************************************************************










    !>  Update the primary field 'face interior' cache entries.
    !!
    !!  Computes the 'value' and 'gradient' entries.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/7/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_primary_interior(self,worker,equation_set,bc_state_group,differentiate,compute_gradients)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),     intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate
        logical,                    intent(in)      :: compute_gradients

        integer(ik)                                 :: idepend, ieqn, idomain_l, ielement_l, iface, idiff, eqn_ID
        character(:),   allocatable                 :: field
        real(rk),       allocatable                 :: ale_Dinv(:,:,:)
        real(rk),       allocatable, dimension(:)   :: ale_g, ale_g_grad1, ale_g_grad2, ale_g_grad3
        type(AD_D),     allocatable, dimension(:)   :: value_u, grad1_u, grad2_u, grad3_u, grad1_tmp, grad2_tmp, grad3_tmp


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface

        !
        ! Face interior state. 'values' only depends on interior element.
        !
        idepend = 1


        !
        ! Set differentiation indicator
        !
        if (differentiate) then
            idiff = DIAG
        else
            idiff = 0
        end if


        !
        ! Compute Values
        !
        worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,idiff)
        worker%function_info%idepend = idepend
        worker%function_info%idiff   = idiff
        eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
        do ieqn = 1,worker%mesh%domain(idomain_l)%neqns

            field = worker%prop(eqn_ID)%get_primary_field_name(ieqn)

            
            !
            ! Interpolate modes to nodes on undeformed element
            ! NOTE: we always need to compute the graduent for interior faces for boundary conditions.
            !
            value_u = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ieqn,worker%itime,'value',ME)
            grad1_u = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ieqn,worker%itime,'grad1',ME)
            grad2_u = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ieqn,worker%itime,'grad2',ME)
            grad3_u = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ieqn,worker%itime,'grad3',ME)


            !
            ! Get ALE transformation data
            !
            ale_g       = worker%get_det_jacobian_grid_face('value','face interior')
            ale_g_grad1 = worker%get_det_jacobian_grid_face('grad1','face interior')
            ale_g_grad2 = worker%get_det_jacobian_grid_face('grad2','face interior')
            ale_g_grad3 = worker%get_det_jacobian_grid_face('grad3','face interior')
            ale_Dinv    = worker%get_inv_jacobian_grid_face('face interior')


            !
            ! Compute transformation to deformed element
            !
            value_u   = value_u/ale_g
            grad1_tmp = grad1_u-(value_u)*ale_g_grad1
            grad2_tmp = grad2_u-(value_u)*ale_g_grad2
            grad3_tmp = grad3_u-(value_u)*ale_g_grad3

            grad1_u = (ale_Dinv(1,1,:)*grad1_tmp + ale_Dinv(2,1,:)*grad2_tmp + ale_Dinv(3,1,:)*grad3_tmp)/ale_g
            grad2_u = (ale_Dinv(1,2,:)*grad1_tmp + ale_Dinv(2,2,:)*grad2_tmp + ale_Dinv(3,2,:)*grad3_tmp)/ale_g
            grad3_u = (ale_Dinv(1,3,:)*grad1_tmp + ale_Dinv(2,3,:)*grad2_tmp + ale_Dinv(3,3,:)*grad3_tmp)/ale_g


            !
            ! Store quantities valid on the deformed element
            !
            call worker%cache%set_data(field,'face interior',value_u,'value',   0,worker%function_info%seed,iface)
            call worker%cache%set_data(field,'face interior',grad1_u,'gradient',1,worker%function_info%seed,iface)
            call worker%cache%set_data(field,'face interior',grad2_u,'gradient',2,worker%function_info%seed,iface)
            call worker%cache%set_data(field,'face interior',grad3_u,'gradient',3,worker%function_info%seed,iface)

        end do !ieqn




    end subroutine update_primary_interior
    !*****************************************************************************************










    !>  Update the primary field 'face exterior' cache entries.
    !!
    !!  Computes the 'value' and 'gradient' entries.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/7/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_primary_exterior(self,worker,equation_set,bc_state_group,differentiate,compute_gradients)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),     intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate
        logical,                    intent(in)      :: compute_gradients

        integer(ik)                                 :: idepend, ifield, idomain_l, ielement_l, &
                                                       iface, BC_ID, BC_face, ndepend, idiff, eqn_ID, ChiID
        character(:),   allocatable                 :: field
        real(rk),       allocatable                 :: ale_Dinv(:,:,:)
        real(rk),       allocatable, dimension(:)   :: ale_g, ale_g_grad1, ale_g_grad2, ale_g_grad3
        type(AD_D),     allocatable, dimension(:)   :: value_u, grad1_u, grad2_u, grad3_u, grad1_tmp, grad2_tmp, grad3_tmp


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface


        !
        ! Set differentiation indicator
        !
        if (differentiate) then
            idiff = iface
        else
            idiff = 0
        end if


        ! 
        ! Compute the number of exterior element dependencies for face exterior state
        !
        ndepend = get_ndepend_exterior(worker,equation_set,bc_state_group,differentiate)



        !
        ! Face exterior state. Value
        !
        if ( (worker%face_type() == INTERIOR) .or. &
             (worker%face_type() == CHIMERA ) ) then
            


            if (worker%face_type() == INTERIOR) then
                eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
            else if (worker%face_type() == CHIMERA) then
                ChiID = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%ChiID
                eqn_ID = worker%mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(1)%eqn_ID
            end if


            do ifield = 1,worker%prop(eqn_ID)%nprimary_fields()
                field = worker%prop(eqn_ID)%get_primary_field_name(ifield)
                do idepend = 1,ndepend

                    worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,idiff)
                    worker%function_info%idepend = idepend


                    !
                    ! Interpolate U
                    !

                    ! Interpolate modes to nodes on undeformed element
                    value_u = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ifield,worker%itime,'value',NEIGHBOR)

                    ! Get ALE transformation data
                    ale_g = worker%get_det_jacobian_grid_face('value','face exterior')

                    ! Compute transformation to deformed element
                    value_u = (value_u/ale_g)

                    ! Store quantities valid on the deformed element
                    call worker%cache%set_data(field,'face exterior',value_u,'value',0,worker%function_info%seed,iface)



                    !
                    ! Interpolate Grad(U)
                    !
                    if (compute_gradients) then
                        !
                        ! Interpolate modes to nodes on undeformed element
                        !
                        grad1_u = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ifield,worker%itime,'grad1',NEIGHBOR)
                        grad2_u = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ifield,worker%itime,'grad2',NEIGHBOR)
                        grad3_u = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ifield,worker%itime,'grad3',NEIGHBOR)


                        !
                        ! Get ALE transformation data
                        !
                        ale_g_grad1 = worker%get_det_jacobian_grid_face('grad1','face exterior')
                        ale_g_grad2 = worker%get_det_jacobian_grid_face('grad2','face exterior')
                        ale_g_grad3 = worker%get_det_jacobian_grid_face('grad3','face exterior')
                        ale_Dinv    = worker%get_inv_jacobian_grid_face('face exterior')

                        
                        !
                        ! Compute transformation to deformed element
                        !
                        grad1_tmp = grad1_u-(value_u)*ale_g_grad1
                        grad2_tmp = grad2_u-(value_u)*ale_g_grad2
                        grad3_tmp = grad3_u-(value_u)*ale_g_grad3

                        grad1_u = (ale_Dinv(1,1,:)*grad1_tmp + ale_Dinv(2,1,:)*grad2_tmp + ale_Dinv(3,1,:)*grad3_tmp)/ale_g
                        grad2_u = (ale_Dinv(1,2,:)*grad1_tmp + ale_Dinv(2,2,:)*grad2_tmp + ale_Dinv(3,2,:)*grad3_tmp)/ale_g
                        grad3_u = (ale_Dinv(1,3,:)*grad1_tmp + ale_Dinv(2,3,:)*grad2_tmp + ale_Dinv(3,3,:)*grad3_tmp)/ale_g


                        ! 
                        ! Store quantities valid on the deformed element
                        !
                        call worker%cache%set_data(field,'face exterior',grad1_u,'gradient',1,worker%function_info%seed,iface)
                        call worker%cache%set_data(field,'face exterior',grad2_u,'gradient',2,worker%function_info%seed,iface)
                        call worker%cache%set_data(field,'face exterior',grad3_u,'gradient',3,worker%function_info%seed,iface)
                    end if


                end do !idepend
            end do !ifield

        end if



    end subroutine update_primary_exterior
    !*****************************************************************************************








    !>  Update the primary field BOUNDARY state functions. These are placed in the 
    !!  'face exterior' cache entry.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/7/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_primary_bc(self,worker,equation_set,bc_state_group,differentiate)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),     intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate

        integer(ik)                 :: idepend, ieqn, idomain_l, ielement_l, iface, ndepend, &
                                       istate, bc_ID, group_ID, patch_ID, face_ID, eqn_ID
        character(:),   allocatable :: field


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface


        !
        ! Face bc(exterior) state
        !
        if ( (worker%face_type() == BOUNDARY)  ) then
            
            bc_ID    = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%bc_ID
            group_ID = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%group_ID
            patch_ID = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%patch_ID
            face_ID  = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%face_ID

            ndepend = get_ndepend_exterior(worker,equation_set,bc_state_group,differentiate)

            do istate = 1,size(bc_state_group(bc_ID)%bc_state)
                do idepend = 1,ndepend

                    ! Get coupled bc element to linearize against.
                    if (differentiate) then
                        worker%function_info%seed%idomain_g  = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_g(idepend)
                        worker%function_info%seed%idomain_l  = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_l(idepend)
                        worker%function_info%seed%ielement_g = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_g(idepend)
                        worker%function_info%seed%ielement_l = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_l(idepend)
                        worker%function_info%seed%iproc      = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%proc(idepend)
                    else
                        worker%function_info%seed%idomain_g  = 0
                        worker%function_info%seed%idomain_l  = 0
                        worker%function_info%seed%ielement_g = 0
                        worker%function_info%seed%ielement_l = 0
                        worker%function_info%seed%iproc      = NO_PROC
                    end if

                    eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID

                    call bc_state_group(bc_ID)%bc_state(istate)%state%compute_bc_state(worker,equation_set(eqn_ID)%prop, bc_state_group(bc_ID)%bc_COMM)

                end do !idepend
            end do !istate


        end if



    end subroutine update_primary_bc
    !*****************************************************************************************











    !>  Update the primary field lift functions for diffusion.
    !!
    !!  This only gets computed if there are diffusive operators allocated to the 
    !!  equation set. If not, then there is no need for the lifting operators and they
    !!  are not computed.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/9/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_primary_lift(self,worker,equation_set,bc_state_group,differentiate)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),                 intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate

        integer(ik) :: idomain_l, ielement_l, eqn_ID



        !
        ! Update lifting terms for gradients if diffusive operators are present
        !
        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l
        eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
        if (allocated(equation_set(eqn_ID)%volume_diffusive_operator) .or. &
            allocated(equation_set(eqn_ID)%boundary_diffusive_operator)) then

            call self%update_lift_faces_internal(worker,equation_set,bc_state_group,differentiate)
            call self%update_lift_faces_external(worker,equation_set,bc_state_group,differentiate)

        end if

    end subroutine update_primary_lift
    !*****************************************************************************************











    !>  Update the auxiliary field 'element' cache entries.
    !!
    !!  Computes the 'value' and 'gradient' entries.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/9/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_auxiliary_element(self,worker,equation_set,bc_state_group,differentiate)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),                 intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate

        integer(ik)                                 :: idepend, ieqn, idomain_l, ielement_l, iface, &
                                                       idiff, iaux_field, ifield, eqn_ID
        character(:),   allocatable                 :: field
        type(AD_D),     allocatable, dimension(:)   :: value_gq, grad1_gq, grad2_gq, grad3_gq


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface


        !
        ! Element primary fields volume 'value' cache. Only depends on interior element
        !
        if (differentiate) then
            idiff = DIAG
        else
            idiff = 0
        end if

        idepend = 0 ! no linearization
        eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
        do ifield = 1,worker%prop(eqn_ID)%nauxiliary_fields()

            !
            ! Try to find the auxiliary field in the solverdata_t container; where they are stored.
            !
            field      = worker%prop(eqn_ID)%get_auxiliary_field_name(ifield)
            iaux_field = worker%solverdata%get_auxiliary_field_index(field)

            ! Set seed
            worker%function_info%seed    = element_compute_seed(worker%mesh,idomain_l,ielement_l,idepend,idiff)
            worker%function_info%idepend = idepend
            worker%function_info%idiff   = idiff

            ! Interpolate modes to nodes
            ieqn = 1    !implicitly assuming only 1 equation in the auxiliary field chidgVector
            value_gq = interpolate_element_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%element_info,worker%function_info,ieqn,worker%itime,'value')
            grad1_gq = interpolate_element_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%element_info,worker%function_info,ieqn,worker%itime,'grad1')
            grad2_gq = interpolate_element_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%element_info,worker%function_info,ieqn,worker%itime,'grad2')
            grad3_gq = interpolate_element_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%element_info,worker%function_info,ieqn,worker%itime,'grad3')

            ! Store gq data in cache
            call worker%cache%set_data(field,'element',value_gq,'value',   0,worker%function_info%seed)
            call worker%cache%set_data(field,'element',grad1_gq,'gradient',1,worker%function_info%seed)
            call worker%cache%set_data(field,'element',grad2_gq,'gradient',2,worker%function_info%seed)
            call worker%cache%set_data(field,'element',grad3_gq,'gradient',3,worker%function_info%seed)

        end do !ieqn




    end subroutine update_auxiliary_element
    !*****************************************************************************************









    !>  Update the auxiliary field 'face interior' cache entries.
    !!
    !!  Computes the 'value' and 'gradient' entries.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/7/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_auxiliary_interior(self,worker,equation_set,bc_state_group,differentiate)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),                 intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate

        integer(ik)                                 :: idepend, ifield, idomain_l, ielement_l, iface, &
                                                       iaux_field, iaux, idiff, eqn_ID
        character(:),   allocatable                 :: field
        type(AD_D),     allocatable, dimension(:)   :: value_gq, grad1_gq, grad2_gq, grad3_gq


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface


        !
        ! Set differentiation indicator
        !
        if (differentiate) then
            idiff = DIAG
        else
            idiff = 0
        end if



        !
        ! Face interior state. 
        !
        eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
        idepend = 0 ! no linearization
        do iaux = 1,worker%prop(eqn_ID)%nauxiliary_fields()

            !
            ! Try to find the auxiliary field in the solverdata_t container; where they are stored.
            !
            field      = worker%prop(eqn_ID)%get_auxiliary_field_name(iaux)
            iaux_field = worker%solverdata%get_auxiliary_field_index(field)

            ! Set seed
            worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,idiff)
            worker%function_info%idepend = idepend
            worker%function_info%idiff   = idiff

            ! Interpolate modes to nodes
            ! NOTE: implicitly assuming only 1 field in the auxiliary field chidg_vector
            ifield = 1
            value_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%face_info(),worker%function_info,ifield,worker%itime,'value',ME)
            grad1_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%face_info(),worker%function_info,ifield,worker%itime,'grad1',ME)
            grad2_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%face_info(),worker%function_info,ifield,worker%itime,'grad2',ME)
            grad3_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%face_info(),worker%function_info,ifield,worker%itime,'grad3',ME)

            ! Store gq data in cache
            call worker%cache%set_data(field,'face interior',value_gq,'value',   0,worker%function_info%seed,iface)
            call worker%cache%set_data(field,'face interior',grad1_gq,'gradient',1,worker%function_info%seed,iface)
            call worker%cache%set_data(field,'face interior',grad2_gq,'gradient',2,worker%function_info%seed,iface)
            call worker%cache%set_data(field,'face interior',grad3_gq,'gradient',3,worker%function_info%seed,iface)

        end do !iaux



    end subroutine update_auxiliary_interior
    !*****************************************************************************************














    !>  Update the auxiliary field 'face exterior' cache entries.
    !!
    !!  Computes the 'value' and 'gradient' entries.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/7/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_auxiliary_exterior(self,worker,equation_set,bc_state_group,differentiate)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),                 intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate

        integer(ik)                                 :: idepend, idomain_l, ielement_l, iface, &
                                                       iaux, iaux_field, ifield, idiff, eqn_ID
        character(:),   allocatable                 :: field
        type(AD_D),     allocatable, dimension(:)   :: value_gq, grad1_gq, grad2_gq, grad3_gq


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface


        !
        ! Set differentiation indicator
        !
        if (differentiate) then
            idiff = DIAG
        else
            idiff = 0
        end if


        !
        ! Face exterior state. 
        !
        if ( (worker%face_type() == INTERIOR) .or. (worker%face_type() == CHIMERA) ) then

            eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
            idepend = 0 ! no linearization
            do iaux = 1,worker%prop(eqn_ID)%nauxiliary_fields()

                !
                ! Try to find the auxiliary field in the solverdata_t container; where they are stored.
                !
                field      = worker%prop(eqn_ID)%get_auxiliary_field_name(iaux)
                iaux_field = worker%solverdata%get_auxiliary_field_index(field)

                ! Set seed
                worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,idiff)
                worker%function_info%idepend = idepend
                worker%function_info%idiff   = idiff

                ! Interpolate modes to nodes
                ! WARNING: implicitly assuming only 1 field in the auxiliary field chidg_vector
                ifield = 1
                value_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%face_info(),worker%function_info,ifield,worker%itime,'value',NEIGHBOR)
                grad1_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%face_info(),worker%function_info,ifield,worker%itime,'grad1',NEIGHBOR)
                grad2_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%face_info(),worker%function_info,ifield,worker%itime,'grad2',NEIGHBOR)
                grad3_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%face_info(),worker%function_info,ifield,worker%itime,'grad3',NEIGHBOR)

                ! Store gq data in cache
                call worker%cache%set_data(field,'face exterior',value_gq,'value',   0,worker%function_info%seed,iface)
                call worker%cache%set_data(field,'face exterior',grad1_gq,'gradient',1,worker%function_info%seed,iface)
                call worker%cache%set_data(field,'face exterior',grad2_gq,'gradient',2,worker%function_info%seed,iface)
                call worker%cache%set_data(field,'face exterior',grad3_gq,'gradient',3,worker%function_info%seed,iface)

            end do !iaux

        end if



    end subroutine update_auxiliary_exterior
    !*****************************************************************************************










    !>  Update the auxiliary field bc(face exterior) cache entries.
    !!
    !!  Computes the 'value' and 'gradient' entries.
    !!
    !!  NOTE: This extrapolates information from the 'face interior' and stores in in the
    !!        'face exterior' cache. These are auxiliary fields so they don't exactly have
    !!        a definition outside the domain. An extrapolation is a reasonable assumption.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/7/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_auxiliary_bc(self,worker,equation_set,bc_state_group,differentiate)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),                 intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate

        integer(ik)                                 :: idepend, ifield, idomain_l, ielement_l, iface, &
                                                       iaux_field, iaux, idiff, eqn_ID
        character(:),   allocatable                 :: field
        type(AD_D),     allocatable, dimension(:)   :: value_gq, grad1_gq, grad2_gq, grad3_gq


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface


        !
        ! Set differentiation indicator
        !
        if (differentiate) then
            idiff = DIAG
        else
            idiff = 0
        end if



        !
        ! Face interior state. 
        !
        if ( (worker%face_type() == BOUNDARY) ) then

            eqn_ID  = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
            idepend = 0 ! no linearization
            do iaux = 1,worker%prop(eqn_ID)%nauxiliary_fields()

                !
                ! Try to find the auxiliary field in the solverdata_t container; where they are stored.
                !
                field      = worker%prop(eqn_ID)%get_auxiliary_field_name(iaux)
                iaux_field = worker%solverdata%get_auxiliary_field_index(field)

                ! Set seed
                worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,idiff)
                worker%function_info%idepend = idepend
                worker%function_info%idiff   = idiff

                !
                ! Interpolate modes to nodes
                ifield = 1    !implicitly assuming only 1 equation in the auxiliary field chidgVector
                value_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%face_info(),worker%function_info,ifield,worker%itime,'value',ME)
                grad1_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%face_info(),worker%function_info,ifield,worker%itime,'grad1',ME)
                grad2_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%face_info(),worker%function_info,ifield,worker%itime,'grad2',ME)
                grad3_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%auxiliary_field(iaux_field),worker%face_info(),worker%function_info,ifield,worker%itime,'grad3',ME)

                ! Store gq data in cache
                call worker%cache%set_data(field,'face exterior',value_gq,'value',   0,worker%function_info%seed,iface)
                call worker%cache%set_data(field,'face exterior',grad1_gq,'gradient',1,worker%function_info%seed,iface)
                call worker%cache%set_data(field,'face exterior',grad2_gq,'gradient',2,worker%function_info%seed,iface)
                call worker%cache%set_data(field,'face exterior',grad3_gq,'gradient',3,worker%function_info%seed,iface)

            end do !iaux

        end if



    end subroutine update_auxiliary_bc
    !*****************************************************************************************








    !>  Update the model field 'element' cache entries.
    !!
    !!  Computes the 'value' and 'gradient' entries.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/9/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_model_element(self,worker,equation_set,bc_state_group,differentiate,model_type)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),                 intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate
        character(*),               intent(in)      :: model_type

        logical                     :: diff_none, diff_interior, diff_exterior, compute_model
        integer(ik)                 :: imodel, idomain_l, ielement_l, idepend, idiff, &
                                       ipattern, ndepend, eqn_ID
        integer(ik),    allocatable :: compute_pattern(:)
        character(:),   allocatable :: dependency

        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 


        !
        ! Compute element model field. Potentially differentiated wrt exterior elements.
        !
        eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
        worker%interpolation_source = 'element'
        do imodel = 1,equation_set(eqn_ID)%nmodels()

            !
            ! Get model dependency
            !
            dependency = equation_set(eqn_ID)%models(imodel)%model%get_dependency()

            !
            ! Only execute models specified in incoming model_type
            !
            if (trim(dependency) == trim(model_type)) then

                !
                ! Determine pattern to compute functions. Depends on if we are differentiating 
                ! or not. These will be used to set idiff, indicating the differentiation
                ! direction.
                !
                if (differentiate) then
                    ! compute function, wrt (all exterior)/interior states
                    if (dependency == 'f(Q-)') then
                        compute_pattern = [DIAG]
                    else if ( (dependency == 'f(Q-,Q+)') .or. &
                              (dependency == 'f(Grad(Q))') ) then
                        compute_pattern = [1,2,3,4,5,6,DIAG]
                    else
                        call chidg_signal(FATAL,"cache_handler%update_model_element: Invalid model dependency string.")
                    end if
                else
                    ! compute function, but do not differentiate
                    compute_pattern = [0]
                end if




                !
                ! Execute compute pattern
                !
                do ipattern = 1,size(compute_pattern)

                
                    !
                    ! get differentiation indicator
                    !
                    idiff = compute_pattern(ipattern)

                    diff_none     = (idiff == 0)
                    diff_interior = (idiff == DIAG)
                    diff_exterior = ( (idiff == 1) .or. (idiff == 2) .or. &
                                      (idiff == 3) .or. (idiff == 4) .or. &
                                      (idiff == 5) .or. (idiff == 6) )



                    if (diff_interior .or. diff_none) then
                        compute_model = .true.
                    else if (diff_exterior) then
                        compute_model = ( (worker%mesh%domain(idomain_l)%faces(ielement_l,idiff)%ftype == INTERIOR) .or. &
                                          (worker%mesh%domain(idomain_l)%faces(ielement_l,idiff)%ftype == CHIMERA) )
                    end if



                    if (compute_model) then

                        if (diff_none .or. diff_interior) then
                            ndepend = 1
                        else
                            call worker%set_face(idiff)
                            ndepend = get_ndepend_exterior(worker,equation_set,bc_state_group,differentiate)
                        end if

                        do idepend = 1,ndepend
                            worker%function_info%seed    = element_compute_seed(worker%mesh,idomain_l,ielement_l,idepend,idiff)
                            worker%function_info%idepend = idepend

                            call equation_set(eqn_ID)%models(imodel)%model%compute(worker)
                        end do !idepend
                    end if !compute

                end do !ipattern

            end if ! select model type
        end do !imodel


    end subroutine update_model_element
    !*****************************************************************************************










    !>  Update the model field 'value', 'face interior' cache entries.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/7/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_model_interior(self,worker,equation_set,bc_state_group,differentiate,model_type)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),                 intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate
        character(*),               intent(in)      :: model_type

        logical                     :: exterior_coupling, selected_model
        integer(ik)                 :: idepend, imodel, idomain_l, ielement_l, &
                                       iface, idiff, ndepend, eqn_ID
        integer(ik),    allocatable :: compute_pattern(:)
        character(:),   allocatable :: field, model_dependency, mode
        type(AD_D),     allocatable :: value_gq(:)


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface


        !
        ! Update models for 'face interior'. Differentiated wrt interior.
        !
        idepend = 1
        eqn_ID  = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
        worker%interpolation_source = 'face interior'
        do imodel = 1,equation_set(eqn_ID)%nmodels()

            !
            ! Compute if model dependency matches specified model type in the 
            ! function interface.
            !
            model_dependency = equation_set(eqn_ID)%models(imodel)%model%get_dependency()
            selected_model = (trim(model_type) == trim(model_dependency))

            if (selected_model) then

                !
                ! Set differentiation indicator
                !
                if (differentiate) then
                    idiff = DIAG
                else
                    idiff = 0 
                end if


                worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,idiff)
                worker%function_info%idepend = idepend
                worker%function_info%idiff   = idiff

                call equation_set(eqn_ID)%models(imodel)%model%compute(worker)

            end if !select model

        end do !imodel




        !
        ! Update models for 'face interior'. Differentiated wrt exterior.
        !
        worker%interpolation_source = 'face interior'
        if ( (worker%face_type() == INTERIOR) .or. &
             (worker%face_type() == CHIMERA) ) then

            if (differentiate) then

                do imodel = 1,equation_set(eqn_ID)%nmodels()

                    !
                    ! Get model dependency 
                    !
                    model_dependency = equation_set(eqn_ID)%models(imodel)%model%get_dependency()

                    selected_model    = (trim(model_type) == trim(model_dependency))
                    exterior_coupling = (model_dependency == 'f(Q-,Q+)') .or. (model_dependency == 'f(Grad(Q))')

                    if ( selected_model .and. exterior_coupling ) then

                        !
                        ! Set differentiation indicator
                        !
                        idiff = iface

                        ! 
                        ! Compute the number of exterior element dependencies
                        !
                        ndepend = get_ndepend_exterior(worker,equation_set,bc_state_group,differentiate)

                        !
                        ! Loop through external dependencies and compute model
                        !
                        do idepend = 1,ndepend
                            worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,idiff)
                            worker%function_info%idepend = idepend
                            worker%function_info%idiff   = idiff

                            call equation_set(eqn_ID)%models(imodel)%model%compute(worker)
                        end do !idepend

                    end if ! select model


                end do !imodel

            end if !differentiate
        end if ! INTERIOR or CHIMERA



    end subroutine update_model_interior
    !*****************************************************************************************











    !>  Update the model field 'value', 'face exterior' cache entries.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/7/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine update_model_exterior(self,worker,equation_set,bc_state_group,differentiate,model_type)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),                 intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate
        character(*),               intent(in)      :: model_type

        integer(ik)                 :: idepend, imodel, idomain_l, ielement_l, iface, &
                                       bc_ID, patch_ID, face_ID, ndepend, idiff, eqn_ID, ChiID
        character(:),   allocatable :: field, model_dependency
        logical                     :: selected_model


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface



        !
        ! Face exterior state: interior neighbors and chimera
        !
        !eqn_ID = worker%mesh%domain(idomain_l)%eqn_ID
        if (worker%face_type() == INTERIOR) then
            eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
        else if (worker%face_type() == CHIMERA) then
            ChiID = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%ChiID
            eqn_ID = worker%mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(1)%eqn_ID
        end if




        worker%interpolation_source = 'face exterior'
        if ( (worker%face_type() == INTERIOR) .or. (worker%face_type() == CHIMERA) ) then

            !
            ! Set differentiation indicator. Differentiate 'face exterior' wrt EXTERIOR elements
            !
            if (differentiate) then
                idiff = iface
            else
                idiff = 0
            end if
            
            ! 
            ! Compute the number of exterior element dependencies for face exterior state
            !
            ndepend = get_ndepend_exterior(worker,equation_set,bc_state_group,differentiate)
            do imodel = 1,equation_set(eqn_ID)%nmodels()

                !
                ! Get model dependency 
                !
                model_dependency = equation_set(eqn_ID)%models(imodel)%model%get_dependency()
                selected_model   = (trim(model_type) == trim(model_dependency))

                if (selected_model) then
                    do idepend = 1,ndepend

                        worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,idiff)
                        worker%function_info%idepend = idepend

                        call equation_set(eqn_ID)%models(imodel)%model%compute(worker)

                    end do !idepend
                end if !select model

            end do !imodel


            !
            ! Set differentiation indicator. Differentiate 'face exterior' wrt INTERIOR element
            ! Only need to compute if differentiating
            !
            if (differentiate) then

                idiff = DIAG
            
                ! 
                ! Compute the number of exterior element dependencies for face exterior state
                !
                ndepend = 1
                do imodel = 1,equation_set(eqn_ID)%nmodels()

                    !
                    ! Get model dependency 
                    !
                    model_dependency = equation_set(eqn_ID)%models(imodel)%model%get_dependency()
                    selected_model   = (trim(model_type) == trim(model_dependency))

                    if (selected_model) then
                        do idepend = 1,ndepend

                            worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,idiff)
                            worker%function_info%idepend = idepend

                            call equation_set(eqn_ID)%models(imodel)%model%compute(worker)

                        end do !idepend
                    end if !select model

                end do !imodel

            end if


        end if ! worker%face_type()

    end subroutine update_model_exterior
    !*****************************************************************************************







    !>  Update the model field BOUNDARY state functions. These are placed in the 
    !!  'face exterior' cache entry.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/9/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_model_bc(self,worker,equation_set,bc_state_group,differentiate,model_type)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),     intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate
        character(*),               intent(in)      :: model_type

        integer(ik)                 :: idepend, ieqn, idomain_l, ielement_l, iface, ndepend, &
                                       istate, bc_ID, group_ID, patch_ID, face_ID, imodel, eqn_ID
        character(:),   allocatable :: field, model_dependency
        logical                     :: selected_model


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface



        !
        ! Face exterior state: boundaries
        !
        eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
        worker%interpolation_source = 'face exterior'
        if ( (worker%face_type() == BOUNDARY) ) then

            ! 
            ! Compute the number of exterior element dependencies for face exterior state
            !
            ndepend = get_ndepend_exterior(worker,equation_set,bc_state_group,differentiate)

            bc_ID    = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%bc_ID
            group_ID = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%group_ID
            patch_ID = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%patch_ID
            face_ID  = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%face_ID

            do imodel = 1,equation_set(eqn_ID)%nmodels()

                !
                ! Get model dependency 
                !
                model_dependency = equation_set(eqn_ID)%models(imodel)%model%get_dependency()
                selected_model   = (trim(model_type) == trim(model_dependency))

                if (selected_model) then
                    do idepend = 1,ndepend


                        if (differentiate) then
                            ! Get coupled bc element to differentiate wrt
                            worker%function_info%seed%idomain_g  = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_g(idepend)
                            worker%function_info%seed%idomain_l  = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_l(idepend)
                            worker%function_info%seed%ielement_g = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_g(idepend)
                            worker%function_info%seed%ielement_l = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_l(idepend)
                            worker%function_info%seed%iproc      = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%proc(idepend)

                        else
                            ! Set no differentiation
                            worker%function_info%seed%idomain_g  = 0
                            worker%function_info%seed%idomain_l  = 0
                            worker%function_info%seed%ielement_g = 0
                            worker%function_info%seed%ielement_l = 0
                            worker%function_info%seed%iproc      = NO_PROC
                        end if



                        call equation_set(eqn_ID)%models(imodel)%model%compute(worker)

                    end do !idepend
                end if !select model

            end do !imodel


        end if ! worker%face_type()





    end subroutine update_model_bc
    !*****************************************************************************************















    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine update_lift_faces_internal(self,worker,equation_set,bc_state_group,differentiate)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),     intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate

        character(:),   allocatable :: field
        real(rk),       allocatable :: ale_g_m(:), ale_g_p(:)
        integer(ik)                 :: idomain_l, ielement_l, iface, idepend, &
                                       ndepend, BC_ID, BC_face, ifield, idiff, eqn_ID

        type(AD_D), allocatable, dimension(:), save   ::    &
            var_m, var_p, var_diff, var_diff_weighted,      &
            var_diff_x,     var_diff_y,     var_diff_z,     &
            rhs_x,          rhs_y,          rhs_z,          &
            lift_modes_x,   lift_modes_y,   lift_modes_z,   &
            lift_gq_face_x, lift_gq_face_y, lift_gq_face_z, &
            lift_gq_vol_x,  lift_gq_vol_y,  lift_gq_vol_z



        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 





        !
        ! For each face, compute the lifting operators associated with each equation for the 
        ! internal and external states and also their linearization.
        !
        do iface = 1,NFACES


            !
            ! Update worker face index
            !
            call worker%set_face(iface)



            associate ( weights          => worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%weights_face(iface),                            &
                        val_face_trans   => transpose(worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%interpolator_face('Value',iface)),    &
                        val_face         => worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%interpolator_face('Value',iface),               &
                        val_vol          => worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%interpolator_element('Value'),                  &
                        invmass          => worker%mesh%domain(idomain_l)%elems(ielement_l)%invmass,                                                &
                        br2_face         => worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%br2_face,                                         &
                        br2_vol          => worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%br2_vol)





            eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
            do ifield = 1,worker%prop(eqn_ID)%nprimary_fields()

                !
                ! Get field
                !
                field = worker%prop(eqn_ID)%get_primary_field_name(ifield)



                !
                ! Compute Interior lift, differentiated wrt Interior
                !

                ! Set differentiation indicator
                if (differentiate) then
                    idiff = DIAG
                else
                    idiff = 0
                end if

                ndepend = 1
                do idepend = 1,ndepend

                    ! Get Seed
                    worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,idiff)
                    worker%function_info%idepend = idepend


                    ! Get interior/exterior state on deformed element
                    var_m = worker%cache%get_data(field,'face interior', 'value', 0, worker%function_info%seed, iface)
                    var_p = worker%cache%get_data(field,'face exterior', 'value', 0, worker%function_info%seed, iface)


                    ! Get ALE transformation
                    ale_g_m = worker%get_det_jacobian_grid_face('value', 'face interior')
                    ale_g_p = worker%get_det_jacobian_grid_face('value', 'face exterior')

                    ! Transform values to undeformed element
                    var_m = var_m*ale_g_m
                    var_p = var_p*ale_g_p


                    ! Difference
                    var_diff = HALF*(var_p - var_m) 

                    ! Multiply by weights
                    var_diff_weighted = var_diff * weights

                    ! Multiply by normal. Note: normal is scaled by face jacobian.
                    var_diff_x = var_diff_weighted * worker%normal(1)
                    var_diff_y = var_diff_weighted * worker%normal(2)
                    var_diff_z = var_diff_weighted * worker%normal(3)


                    !
                    ! Standard Approach breaks the process up into several steps:
                    !   1: Project onto basis
                    !   2: Local solve for lift modes in element basis
                    !   3: Interpolate lift modes to face/volume quadrature nodes
                    !
                    ! Improved approach creates a single matrix that performs the
                    ! three steps in one MV multiply:
                    !
                    !   br2_face = [val_face][invmass][val_face_trans]
                    !   br2_vol  = [val_vol ][invmass][val_face_trans]
                    !
                    lift_gq_face_x = matmul(br2_face,var_diff_x)
                    lift_gq_face_y = matmul(br2_face,var_diff_y)
                    lift_gq_face_z = matmul(br2_face,var_diff_z)

                    
                    ! Store lift
                    call worker%cache%set_data(field,'face interior', lift_gq_face_x, 'lift face', 1, worker%function_info%seed, iface)
                    call worker%cache%set_data(field,'face interior', lift_gq_face_y, 'lift face', 2, worker%function_info%seed, iface)
                    call worker%cache%set_data(field,'face interior', lift_gq_face_z, 'lift face', 3, worker%function_info%seed, iface)


                    ! 1: Project onto element basis
                    ! 2: Interpolate lift modes to volume quadrature nodes
                    lift_gq_vol_x = matmul(br2_vol,var_diff_x)
                    lift_gq_vol_y = matmul(br2_vol,var_diff_y)
                    lift_gq_vol_z = matmul(br2_vol,var_diff_z)

                    
                    ! Store lift
                    call worker%cache%set_data(field,'face interior', lift_gq_vol_x, 'lift element', 1, worker%function_info%seed, iface)
                    call worker%cache%set_data(field,'face interior', lift_gq_vol_y, 'lift element', 2, worker%function_info%seed, iface)
                    call worker%cache%set_data(field,'face interior', lift_gq_vol_z, 'lift element', 3, worker%function_info%seed, iface)



                end do !idepend







                !
                ! Compute Interior lift, differentiated wrt Exterior
                !

                ! Set differentiation indicator
                if (differentiate) then
                    idiff = iface
                else
                    idiff = 0
                end if
                ndepend = get_ndepend_exterior(worker,equation_set,bc_state_group,differentiate)
                do idepend = 1,ndepend

                    ! Get Seed
                    worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,idiff)
                    worker%function_info%idepend = idepend


                    ! Get interior/exterior state
                    var_m = worker%cache%get_data(field,'face interior', 'value', 0, worker%function_info%seed, iface)
                    var_p = worker%cache%get_data(field,'face exterior', 'value', 0, worker%function_info%seed, iface)


                    ! Get ALE transformation
                    ale_g_m = worker%get_det_jacobian_grid_face('value', 'face interior')
                    ale_g_p = worker%get_det_jacobian_grid_face('value', 'face exterior')

                    ! Transform values to undeformed element
                    var_m = var_m*ale_g_m
                    var_p = var_p*ale_g_p


                    ! Difference
                    var_diff = HALF*(var_p - var_m) 


                    ! Multiply by weights
                    var_diff_weighted = var_diff * weights

                    ! Multiply by normal. Note: normal is scaled by face jacobian.
                    var_diff_x = var_diff_weighted * worker%normal(1)
                    var_diff_y = var_diff_weighted * worker%normal(2)
                    var_diff_z = var_diff_weighted * worker%normal(3)

                    ! Project onto element basis, evaluate at face quadrature nodes
                    lift_gq_face_x = matmul(br2_face,var_diff_x)
                    lift_gq_face_y = matmul(br2_face,var_diff_y)
                    lift_gq_face_z = matmul(br2_face,var_diff_z)
                    
                    ! Store lift
                    call worker%cache%set_data(field,'face interior', lift_gq_face_x, 'lift face', 1, worker%function_info%seed, iface)
                    call worker%cache%set_data(field,'face interior', lift_gq_face_y, 'lift face', 2, worker%function_info%seed, iface)
                    call worker%cache%set_data(field,'face interior', lift_gq_face_z, 'lift face', 3, worker%function_info%seed, iface)


                    ! Project onto element basis, evaluate at element quadrature nodes
                    lift_gq_vol_x = matmul(br2_vol,var_diff_x)
                    lift_gq_vol_y = matmul(br2_vol,var_diff_y)
                    lift_gq_vol_z = matmul(br2_vol,var_diff_z)
                    
                    ! Store lift
                    call worker%cache%set_data(field,'face interior', lift_gq_vol_x, 'lift element', 1, worker%function_info%seed, iface)
                    call worker%cache%set_data(field,'face interior', lift_gq_vol_y, 'lift element', 2, worker%function_info%seed, iface)
                    call worker%cache%set_data(field,'face interior', lift_gq_vol_z, 'lift element', 3, worker%function_info%seed, iface)



                end do !idepend

            end do !ifield


            end associate

        end do !iface


    end subroutine update_lift_faces_internal
    !*****************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_lift_faces_external(self,worker,equation_set,bc_state_group,differentiate)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),     intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate

        integer(ik) :: idomain_l, ielement_l, iface, idepend, ieqn, &
                       ndepend, BC_ID, BC_face, idiff
        logical     :: boundary_face, interior_face, chimera_face


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 


        !
        ! For each face, compute the lifting operators associated with each equation for the 
        ! internal and external states and also their linearization.
        !
        do iface = 1,NFACES

            !
            ! Update worker face index
            !
            call worker%set_face(iface)


            !
            ! Check if boundary or interior
            !
            boundary_face = (worker%face_type() == BOUNDARY)
            interior_face = (worker%face_type() == INTERIOR)
            chimera_face  = (worker%face_type() == CHIMERA )



            !
            ! Compute lift for each equation
            !
            do ieqn = 1,worker%mesh%domain(idomain_l)%neqns


                !
                ! Compute External lift, differentiated wrt Interior
                !

                ! Set differentiation indicator
                if (differentiate) then
                    idiff = DIAG
                else
                    idiff = 0
                end if
                ndepend = 1
                do idepend = 1,ndepend

                    ! Get Seed
                    worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,idiff)
                    worker%function_info%idepend = idepend


                    if (interior_face) then
                        call handle_external_lift__interior_face(worker,equation_set,bc_state_group,ieqn)
                    else if (boundary_face) then
                        call handle_external_lift__boundary_face(worker,equation_set,bc_state_group,ieqn)
                    else if (chimera_face) then
                        call handle_external_lift__chimera_face( worker,equation_set,bc_state_group,ieqn)
                    else
                        call chidg_signal(FATAL,"update_lift_faces_external: unsupported face type")
                    end if


                end do !idepend




                !
                ! Compute External lift, differentiated wrt Exterior
                !

                ! Set differentiation indicator
                if (differentiate) then
                    idiff = iface
                else
                    idiff = 0
                end if
                ndepend = get_ndepend_exterior(worker,equation_set,bc_state_group,differentiate)
                do idepend = 1,ndepend

                    ! Get Seed
                    worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,idiff)
                    worker%function_info%idepend = idepend

                    if (interior_face) then
                        call handle_external_lift__interior_face(worker,equation_set,bc_state_group,ieqn)
                    else if (boundary_face) then
                        call handle_external_lift__boundary_face(worker,equation_set,bc_state_group,ieqn)
                    else if (chimera_face) then
                        call handle_external_lift__chimera_face( worker,equation_set,bc_state_group,ieqn)
                    else
                        call chidg_signal(FATAL,"update_lift_faces_external: unsupported face type")
                    end if


                end do !idepend

            end do !ieqn



        end do !iface


    end subroutine update_lift_faces_external
    !*****************************************************************************************
















    !>  Handle computing lift for an external element, when the face is an interior face.
    !!
    !!  In this case, the external element exists and we can just use its data. This is 
    !!  not the case for a boundary condition face, and it is complicated further by a 
    !!  Chimera boundary face.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine handle_external_lift__interior_face(worker,equation_set,bc_state_group,ieqn)
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),     intent(inout)   :: bc_state_group(:)
        integer(ik),                intent(in)      :: ieqn

        integer(ik) :: idomain_l, ielement_l, iface, idomain_l_n, ielement_l_n, iface_n, iproc_n, eqn_ID
        logical     :: boundary_face, interior_face, local_neighbor, remote_neighbor

        type(AD_D), allocatable, dimension(:)   ::          &
            var_m, var_p, var_diff, var_diff_weighted,      &
            var_diff_x,     var_diff_y,     var_diff_z,     &
            rhs_x,          rhs_y,          rhs_z,          &
            lift_modes_x,   lift_modes_y,   lift_modes_z,   &
            lift_gq_face_x, lift_gq_face_y, lift_gq_face_z, &
            lift_gq_vol_x,  lift_gq_vol_y,  lift_gq_vol_z

        character(:),   allocatable                 :: field
        real(rk),       allocatable, dimension(:)   :: normx, normy, normz, weights, ale_g_m, ale_g_p
        real(rk),       allocatable, dimension(:,:) :: val_face_trans, val_face, val_vol, &
                                                       invmass, br2_face


        !
        ! Interior element
        ! 
        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface


        !
        ! Neighbor element
        !
        idomain_l_n  = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_l
        ielement_l_n = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_l
        iface_n      = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_face
        iproc_n      = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_proc

        local_neighbor  = (iproc_n == IRANK)
        remote_neighbor = (iproc_n /= IRANK)


        !
        ! Get field
        !
        eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
        field = worker%prop(eqn_ID)%get_primary_field_name(ieqn)


        if ( local_neighbor ) then
            weights          = worker%mesh%domain(idomain_l_n)%elems(ielement_l_n)%basis_s%weights_face(iface_n)
            val_face_trans   = transpose(worker%mesh%domain(idomain_l_n)%elems(ielement_l_n)%basis_s%interpolator_face('Value',iface_n))
            val_face         = worker%mesh%domain(idomain_l_n)%elems(ielement_l_n)%basis_s%interpolator_face('Value',iface_n)
            val_vol          = worker%mesh%domain(idomain_l_n)%elems(ielement_l_n)%basis_s%interpolator_element('Value')
            invmass          = worker%mesh%domain(idomain_l_n)%elems(ielement_l_n)%invmass
            br2_face         = worker%mesh%domain(idomain_l_n)%faces(ielement_l_n,iface_n)%br2_face


        else if ( remote_neighbor ) then
            ! User local element gq instance. Assumes same order of accuracy.
            weights          = worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%weights_face(iface_n)
            val_face_trans   = transpose(worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%interpolator_face('Value',iface_n))
            val_face         = worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%interpolator_face('Value',iface_n)
            val_vol          = worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%interpolator_element('Value')
            invmass          = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%neighbor_invmass
            br2_face         = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%neighbor_br2_face


        end if



            ! Use reverse of interior element's normal vector
            normx = -worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%norm(:,1)
            normy = -worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%norm(:,2)
            normz = -worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%norm(:,3)

            ! Get interior/exterior state
            var_m = worker%cache%get_data(field,'face interior', 'value', 0, worker%function_info%seed, iface)
            var_p = worker%cache%get_data(field,'face exterior', 'value', 0, worker%function_info%seed, iface)


            ! Get ALE transformation
            ale_g_m = worker%get_det_jacobian_grid_face('value', 'face interior')
            ale_g_p = worker%get_det_jacobian_grid_face('value', 'face exterior')

            ! Transform values to undeformed element
            var_m = var_m*ale_g_m
            var_p = var_p*ale_g_p



            ! Difference. Relative to exterior element, so reversed
            ! Relative to the exterior element, var_m is the exterior state
            ! and var_p is the interior state.
            var_diff = HALF*(var_m - var_p) 

            ! Multiply by weights
            var_diff_weighted = var_diff * weights

            ! Multiply by normal. Note: normal is scaled by face jacobian.
            var_diff_x = var_diff_weighted * normx
            var_diff_y = var_diff_weighted * normy
            var_diff_z = var_diff_weighted * normz

            ! 1: Lift boundary difference. Project into element basis.
            ! 2: Evaluate lift modes at face quadrature nodes
            lift_gq_face_x = matmul(br2_face,var_diff_x)
            lift_gq_face_y = matmul(br2_face,var_diff_y)
            lift_gq_face_z = matmul(br2_face,var_diff_z)
            
            ! Store lift
            call worker%cache%set_data(field,'face exterior', lift_gq_face_x, 'lift face', 1, worker%function_info%seed, iface)
            call worker%cache%set_data(field,'face exterior', lift_gq_face_y, 'lift face', 2, worker%function_info%seed, iface)
            call worker%cache%set_data(field,'face exterior', lift_gq_face_z, 'lift face', 3, worker%function_info%seed, iface)


    end subroutine handle_external_lift__interior_face
    !*****************************************************************************************














    !>  Handle computing lift for an external element, when the face is a boundary face.
    !!
    !!  In this case, the external element does NOT exist, so we use the interior element. 
    !!  !This is kind of like assuming that a boundary element exists of equal size to 
    !!  !the interior element.
    !!
    !!  Actually, on the boundary, we basically just need the interior lift because we aren't
    !!  computing an average flux. Rather we are just computing a boundary flux, so here we
    !!  essentially compute the interior lift and use it for the boundary.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine handle_external_lift__boundary_face(worker,equation_set,bc_state_group,ieqn)
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),     intent(inout)   :: bc_state_group(:)
        integer(ik),                intent(in)      :: ieqn

        integer(ik) :: idomain_l, ielement_l, iface, idomain_l_n, ielement_l_n, iface_n, eqn_ID
        logical     :: boundary_face, interior_face

        type(AD_D), allocatable, dimension(:)   ::          &
            var_m, var_p, var_diff, var_diff_weighted,      &
            var_diff_x,     var_diff_y,     var_diff_z,     &
            rhs_x,          rhs_y,          rhs_z,          &
            lift_modes_x,   lift_modes_y,   lift_modes_z,   &
            lift_gq_x,      lift_gq_y,      lift_gq_z

        character(:),   allocatable                 :: field
        real(rk),       allocatable, dimension(:)   :: normx, normy, normz, ale_g_m, ale_g_p


        !
        ! Interior element
        ! 
        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface


        if (iface == XI_MIN) then
            iface_n = XI_MAX
        else if (iface == ETA_MIN) then
            iface_n = ETA_MAX
        else if (iface == ZETA_MIN) then
            iface_n = ZETA_MAX
        else if (iface == XI_MAX) then
            iface_n = XI_MIN
        else if (iface == ETA_MAX) then
            iface_n = ETA_MIN
        else if (iface == ZETA_MAX) then
            iface_n = ZETA_MIN
        end if


        !
        ! Get field
        !
        eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
        field = worker%prop(eqn_ID)%get_primary_field_name(ieqn)



        !
        ! Neighbor element
        !
        idomain_l_n  = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_l
        ielement_l_n = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_l


        associate ( weights          => worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%weights_face(iface_n),                         &
                    val_face_trans   => transpose(worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%interpolator_face('Value',iface_n)), &
                    val_face         => worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%interpolator_face('Value',iface_n),            &
                    invmass          => worker%mesh%domain(idomain_l)%elems(ielement_l)%invmass,                                        &
                    br2_face         => worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%br2_face)
                    !br2_face         => worker%mesh%domain(idomain_l)%faces(ielement_l,iface_n)%br2_face)

            ! Get normal vector. Use reverse of the normal vector from the interior element since no exterior element exists.
            !normx = -worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%norm(:,1)
            !normy = -worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%norm(:,2)
            !normz = -worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%norm(:,3)
            normx = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%norm(:,1)
            normy = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%norm(:,2)
            normz = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%norm(:,3)

            ! Get interior/exterior state
            var_m = worker%cache%get_data(field,'face interior', 'value', 0, worker%function_info%seed, iface)
            var_p = worker%cache%get_data(field,'face exterior', 'value', 0, worker%function_info%seed, iface)

            ! Get ALE transformation
            ale_g_m = worker%get_det_jacobian_grid_face('value', 'face interior')
            ale_g_p = worker%get_det_jacobian_grid_face('value', 'face exterior')

            ! Transform values to undeformed element
            var_m = var_m*ale_g_m
            var_p = var_p*ale_g_p

            ! Difference. Relative to exterior element, so reversed
            !var_diff = HALF*(var_m - var_p) 
            var_diff = HALF*(var_p - var_m) 


            ! Multiply by weights
            var_diff_weighted = var_diff * weights


            ! Multiply by normal. Note: normal is scaled by face jacobian.
            var_diff_x = var_diff_weighted * normx
            var_diff_y = var_diff_weighted * normy
            var_diff_z = var_diff_weighted * normz


            ! 1: Lift boundary difference. Project into element basis.
            ! 2: Evaluate lift modes at face quadrature nodes
            lift_gq_x = matmul(br2_face,var_diff_x)
            lift_gq_y = matmul(br2_face,var_diff_y)
            lift_gq_z = matmul(br2_face,var_diff_z)
            

            ! Store lift
            call worker%cache%set_data(field,'face exterior', lift_gq_x, 'lift face', 1, worker%function_info%seed, iface)
            call worker%cache%set_data(field,'face exterior', lift_gq_y, 'lift face', 2, worker%function_info%seed, iface)
            call worker%cache%set_data(field,'face exterior', lift_gq_z, 'lift face', 3, worker%function_info%seed, iface)


        end associate






    end subroutine handle_external_lift__boundary_face
    !*****************************************************************************************












    !>  Handle computing lift for an external element, when the face is a Chimera face.
    !!
    !!  In this case, potentially multiple external elements exist, so we don't have just
    !!  a single exterior mass matrix.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine handle_external_lift__chimera_face(worker,equation_set,bc_state_group,ieqn)
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),     intent(inout)   :: bc_state_group(:)
        integer(ik),                intent(in)      :: ieqn

        integer(ik) :: idomain_l, ielement_l, iface, idomain_l_n, ielement_l_n, iface_n, eqn_ID
        logical     :: boundary_face, interior_face

        type(AD_D), allocatable, dimension(:)   ::          &
            var_m, var_p, var_diff, var_diff_weighted,      &
            var_diff_x,     var_diff_y,     var_diff_z,     &
            rhs_x,          rhs_y,          rhs_z,          &
            lift_modes_x,   lift_modes_y,   lift_modes_z,   &
            lift_gq_face_x, lift_gq_face_y, lift_gq_face_z, &
            lift_gq_vol_x,  lift_gq_vol_y,  lift_gq_vol_z

        character(:),   allocatable                 :: field
        real(rk),       allocatable, dimension(:)   :: normx, normy, normz, ale_g_m, ale_g_p


        !
        ! Interior element
        ! 
        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface



        !
        ! Get field
        !
        eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
        field = worker%prop(eqn_ID)%get_primary_field_name(ieqn)

        !
        ! Use components from receiver element since no single element exists to act 
        ! as the exterior element. This implicitly treats the diffusion terms as if 
        ! there were a reflected element like the receiver element that was acting as 
        ! the donor.
        !
        associate ( weights        => worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%weights_face(iface),                            &
                    val_face_trans => transpose(worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%interpolator_face('Value',iface)),    &
                    val_face       => worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%interpolator_face('Value',iface),               &
                    val_vol        => worker%mesh%domain(idomain_l)%elems(ielement_l)%basis_s%interpolator_element('Value'),                  &
                    invmass        => worker%mesh%domain(idomain_l)%elems(ielement_l)%invmass,                                                &
                    br2_face       => worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%br2_face )


            ! Use reversed normal vectors of receiver element
            normx = -worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%norm(:,1)
            normy = -worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%norm(:,2)
            normz = -worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%norm(:,3)

            ! Get interior/exterior state
            var_m = worker%cache%get_data(field,'face interior', 'value', 0, worker%function_info%seed, iface)
            var_p = worker%cache%get_data(field,'face exterior', 'value', 0, worker%function_info%seed, iface)



            ! Get ALE transformation
            ale_g_m = worker%get_det_jacobian_grid_face('value', 'face interior')
            ale_g_p = worker%get_det_jacobian_grid_face('value', 'face exterior')

            ! Transform values to undeformed element
            var_m = var_m*ale_g_m
            var_p = var_p*ale_g_p



            ! Difference. Relative to exterior element, so reversed
            var_diff = HALF*(var_m - var_p) 

            ! Multiply by weights
            var_diff_weighted = var_diff * weights

            ! Multiply by normal. Note: normal is scaled by face jacobian.
            var_diff_x = var_diff_weighted * normx
            var_diff_y = var_diff_weighted * normy
            var_diff_z = var_diff_weighted * normz

            ! 1: Lift boundary difference. Project into element basis.
            ! 2: Evaluate lift modes at face quadrature nodes
            lift_gq_face_x = matmul(br2_face,var_diff_x)
            lift_gq_face_y = matmul(br2_face,var_diff_y)
            lift_gq_face_z = matmul(br2_face,var_diff_z)
            
            ! Store lift
            call worker%cache%set_data(field,'face exterior', lift_gq_face_x, 'lift face', 1, worker%function_info%seed, iface)
            call worker%cache%set_data(field,'face exterior', lift_gq_face_y, 'lift face', 2, worker%function_info%seed, iface)
            call worker%cache%set_data(field,'face exterior', lift_gq_face_z, 'lift face', 3, worker%function_info%seed, iface)


        end associate






    end subroutine handle_external_lift__chimera_face
    !*****************************************************************************************







    !>  For a given state of the chidg_worker(idomain,ielement,iface), return the number
    !!  of exterior dependent elements.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/7/2016
    !!
    !----------------------------------------------------------------------------------------
    function get_ndepend_exterior(worker,equation_set,bc_state_group,differentiate) result(ndepend)
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bc_state_group_t),     intent(inout)   :: bc_state_group(:)
        logical,                    intent(in)      :: differentiate

        integer(ik) :: ndepend, idomain_l, ielement_l, iface, &
                       ChiID, group_ID, patch_ID, face_ID


        if (differentiate) then

            idomain_l  = worker%element_info%idomain_l 
            ielement_l = worker%element_info%ielement_l 
            iface      = worker%iface

            ! 
            ! Compute the number of exterior element dependencies for face exterior state
            !
            if ( worker%face_type() == INTERIOR ) then
                ndepend = 1
                
            else if ( worker%face_type() == CHIMERA ) then
                ChiID   = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%ChiID
                ndepend = worker%mesh%domain(idomain_l)%chimera%recv(ChiID)%ndonors()

            else if ( worker%face_type() == BOUNDARY ) then
                group_ID = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%group_ID
                patch_ID = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%patch_ID
                face_ID  = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%face_ID
                ndepend  = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ncoupled_elements(face_ID)

            end if

        else

            ndepend = 1

        end if

    end function get_ndepend_exterior
    !****************************************************************************************





















end module type_cache_handler
