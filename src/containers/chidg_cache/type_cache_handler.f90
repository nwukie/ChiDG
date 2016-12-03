module type_cache_handler
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: NFACES, INTERIOR, CHIMERA, BOUNDARY, DIAG, ME, NEIGHBOR, HALF, ONE, &
                                  XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use mod_DNAD_tools,     only: face_compute_seed, element_compute_seed
    use mod_interpolate,    only: interpolate_face_autodiff, interpolate_element_autodiff
    use DNAD_D

    use type_chidg_cache,   only: chidg_cache_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_equation_set,  only: equation_set_t
    use type_bcset,         only: bcset_t

    use mod_chidg_mpi,      only: IRANK
    use mpi_f08,            only: MPI_WTime
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

        procedure   :: update

        procedure   :: update_value
        procedure   :: update_derivative
        procedure   :: update_lift
        procedure   :: update_models

        procedure   :: update_lift_faces_internal
        procedure   :: update_lift_faces_external

    end type cache_handler_t
    !****************************************************************************************





contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/7/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update(self,worker,equation_set,bc_set)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bcset_t),              intent(inout)   :: bc_set(:)

        integer(ik) :: idomain_l, ielement_l


        !
        ! Resize cache
        !
        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        call worker%cache%resize(worker%mesh,worker%prop,idomain_l,ielement_l)



        call self%update_value(worker,equation_set,bc_set)
        call self%update_derivative(worker,equation_set,bc_set)
        call self%update_models(worker,equation_set,bc_set)


        !
        ! Update lift terms if diffusive operators are present
        !
        idomain_l = worker%element_info%idomain_l
        if (allocated(equation_set(idomain_l)%volume_diffusive_operator) .or. &
            allocated(equation_set(idomain_l)%boundary_diffusive_operator)) then

            call self%update_lift(worker,equation_set,bc_set)

        end if




    end subroutine update
    !****************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/7/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_value(self,worker,equation_set,bc_set)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bcset_t),              intent(inout)   :: bc_set(:)

        character(:),   allocatable :: field
        integer(ik)                 :: iface, iside, ieqn, idomain_l, ielement_l, idepend, &
                                       ndepend, ChiID, BC_ID, BC_face, ielement_c, istate

        type(AD_D), allocatable, dimension(:) :: value_gq



        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 


        !
        ! Resize cache
        !
        !call worker%cache%resize(worker%mesh,worker%prop,idomain_l,ielement_l)


        !
        ! Loop through faces and cache internal, external interpolated states
        !
        do iface = 1,NFACES



            !
            ! Update worker face index
            !
            call worker%set_face(iface)




            !
            ! Face interior state. 'values' only depends on interior element.
            !
            idepend = 1

            do ieqn = 1,worker%mesh(idomain_l)%neqns


                worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,DIAG)
                worker%function_info%idepend = idepend
                worker%function_info%idiff   = DIAG

                ! Interpolate modes to nodes
                value_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ieqn,'value',ME)

                ! Store gq data in cache
                field = worker%prop(idomain_l)%get_primary_field_name(ieqn)
                call worker%cache%set_data(field,'face interior',value_gq,'value',0,worker%function_info%seed,iface)

            end do !ieqn




            ! 
            ! Compute the number of exterior element dependencies for face exterior state
            !
            if ( worker%face_type() == INTERIOR ) then
                ndepend = 1
                
            else if ( worker%face_type() == CHIMERA ) then
                ChiID   = worker%mesh(idomain_l)%faces(ielement_l,iface)%ChiID
                ndepend = worker%mesh(idomain_l)%chimera%recv%data(ChiID)%ndonors()

            else if ( worker%face_type() == BOUNDARY ) then
                BC_ID   = worker%mesh(idomain_l)%faces(ielement_l,iface)%BC_ID
                BC_face = worker%mesh(idomain_l)%faces(ielement_l,iface)%BC_face
                ndepend = bc_set(idomain_l)%bcs(BC_ID)%get_ncoupled_elems(BC_face)

            end if




            !
            ! Face exterior state
            !
            if ( (worker%face_type() == INTERIOR) .or. (worker%face_type() == CHIMERA) ) then
                
                do ieqn = 1,worker%mesh(idomain_l)%neqns

                    field = worker%prop(idomain_l)%get_primary_field_name(ieqn)

                    do idepend = 1,ndepend

                        worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,iface)
                        worker%function_info%idepend = idepend

                        value_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ieqn,'value',NEIGHBOR)

                        call worker%cache%set_data(field,'face exterior',value_gq,'value',0,worker%function_info%seed,iface)

                    end do !idepend
                end do !ieqn




            else if ( (worker%face_type() == BOUNDARY) ) then
                !
                ! Do nothing here. Boundary condition states(value and derivative) are both updated in 'update_derivative' 
                ! since they are both handled by the bc_state functions.
                !
            end if


        end do !iface





        !
        ! Element volume 'value' cache. Only depends on interior element
        !
        idepend = 1
        do ieqn = 1,worker%mesh(idomain_l)%neqns

                worker%function_info%seed    = element_compute_seed(worker%mesh,idomain_l,ielement_l,idepend,DIAG)
                worker%function_info%idepend = idepend

                value_gq = interpolate_element_autodiff(worker%mesh,worker%solverdata%q,worker%element_info,worker%function_info,ieqn,'value')

                field = worker%prop(idomain_l)%get_primary_field_name(ieqn)
                call worker%cache%set_data(field,'element',value_gq,'value',0,worker%function_info%seed)

        end do !ieqn



    end subroutine update_value
    !************************************************************************************************












    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/15/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine update_derivative(self,worker,equation_set,bc_set)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bcset_t),              intent(inout)   :: bc_set(:)


        character(:),   allocatable :: field
        integer(ik)                 :: iface, iside, ieqn, idomain_l, ielement_l, idepend, &
                                       ndepend, ChiID, BC_ID, BC_face, ielement_c, istate

        type(AD_D), allocatable, dimension(:) :: ddx_gq, ddy_gq, ddz_gq



        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 



        !
        ! Loop through faces and cache internal, external interpolated derivatives
        !
        do iface = 1,NFACES



            !
            ! Update worker face index
            !
            call worker%set_face(iface)




            !
            ! Face interior state. 'values' only depends on interior element.
            !
            idepend = 1

            do ieqn = 1,worker%mesh(idomain_l)%neqns

                worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,DIAG)
                worker%function_info%idepend = idepend
                worker%function_info%idiff   = DIAG

                ! Interpolate modes to nodes
                ddx_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ieqn,'ddx',ME)
                ddy_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ieqn,'ddy',ME)
                ddz_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ieqn,'ddz',ME)

                ! Store gq data in cache
                field = worker%prop(idomain_l)%get_primary_field_name(ieqn)
                call worker%cache%set_data(field,'face interior',ddx_gq,'derivative',1,worker%function_info%seed,iface)
                call worker%cache%set_data(field,'face interior',ddy_gq,'derivative',2,worker%function_info%seed,iface)
                call worker%cache%set_data(field,'face interior',ddz_gq,'derivative',3,worker%function_info%seed,iface)

            end do !ieqn




            ! 
            ! Compute the number of exterior element dependencies for face exterior state
            !
            if ( worker%face_type() == INTERIOR ) then
                ndepend = 1
                
            else if ( worker%face_type() == CHIMERA ) then
                ChiID   = worker%mesh(idomain_l)%faces(ielement_l,iface)%ChiID
                ndepend = worker%mesh(idomain_l)%chimera%recv%data(ChiID)%ndonors()

            else if ( worker%face_type() == BOUNDARY ) then
                BC_ID   = worker%mesh(idomain_l)%faces(ielement_l,iface)%BC_ID
                BC_face = worker%mesh(idomain_l)%faces(ielement_l,iface)%BC_face
                ndepend = bc_set(idomain_l)%bcs(BC_ID)%get_ncoupled_elems(BC_face)

            end if




            !
            ! Face exterior state
            !
            if ( (worker%face_type() == INTERIOR) .or. (worker%face_type() == CHIMERA) ) then
                
                do ieqn = 1,worker%mesh(idomain_l)%neqns

                    field = worker%prop(idomain_l)%get_primary_field_name(ieqn)

                    do idepend = 1,ndepend

                        worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,iface)
                        worker%function_info%idepend = idepend

                        ddx_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ieqn,'ddx',NEIGHBOR)
                        ddy_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ieqn,'ddy',NEIGHBOR)
                        ddz_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ieqn,'ddz',NEIGHBOR)

                        call worker%cache%set_data(field,'face exterior',ddx_gq,'derivative',1,worker%function_info%seed,iface)
                        call worker%cache%set_data(field,'face exterior',ddy_gq,'derivative',2,worker%function_info%seed,iface)
                        call worker%cache%set_data(field,'face exterior',ddz_gq,'derivative',3,worker%function_info%seed,iface)

                    end do !idepend
                end do !ieqn




            else if ( (worker%face_type() == BOUNDARY) ) then


                do istate = 1,size(bc_set(idomain_l)%bcs(BC_ID)%bc_state)
                    do idepend = 1,ndepend


                        !
                        ! Get coupled bc element to linearize against.
                        !
                        ielement_c = bc_set(idomain_l)%bcs(BC_ID)%bc_patch%coupled_elements(BC_face)%at(idepend)
                        worker%function_info%seed%idomain_g  = worker%mesh(idomain_l)%elems(ielement_c)%idomain_g
                        worker%function_info%seed%idomain_l  = worker%mesh(idomain_l)%elems(ielement_c)%idomain_l
                        worker%function_info%seed%ielement_g = worker%mesh(idomain_l)%elems(ielement_c)%ielement_g
                        worker%function_info%seed%ielement_l = worker%mesh(idomain_l)%elems(ielement_c)%ielement_l
                        worker%function_info%seed%iproc      = IRANK

                        call bc_set(idomain_l)%bcs(BC_ID)%bc_state(istate)%state%compute_bc_state(worker,equation_set(idomain_l)%prop)

                    end do !idepend
                end do




            end if


        end do !iface





        !
        ! Element volume 'value' cache. Only depends on interior element
        !
        idepend = 1
        do ieqn = 1,worker%mesh(idomain_l)%neqns

                worker%function_info%seed    = element_compute_seed(worker%mesh,idomain_l,ielement_l,idepend,DIAG)
                worker%function_info%idepend = idepend

                ddx_gq = interpolate_element_autodiff(worker%mesh,worker%solverdata%q,worker%element_info,worker%function_info,ieqn,'ddx')
                ddy_gq = interpolate_element_autodiff(worker%mesh,worker%solverdata%q,worker%element_info,worker%function_info,ieqn,'ddy')
                ddz_gq = interpolate_element_autodiff(worker%mesh,worker%solverdata%q,worker%element_info,worker%function_info,ieqn,'ddz')

                field = worker%prop(idomain_l)%get_primary_field_name(ieqn)
                call worker%cache%set_data(field,"element",ddx_gq,"derivative",1,worker%function_info%seed)
                call worker%cache%set_data(field,"element",ddy_gq,"derivative",2,worker%function_info%seed)
                call worker%cache%set_data(field,"element",ddz_gq,"derivative",3,worker%function_info%seed)

        end do !ieqn



    end subroutine update_derivative
    !************************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/13/2016
    !!
    !!
    !------------------------------------------------------------------------------------------------
    subroutine update_lift(self,worker,equation_set,bc_set)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bcset_t),              intent(inout)   :: bc_set(:)


        call self%update_lift_faces_internal(worker,equation_set,bc_set)

        call self%update_lift_faces_external(worker,equation_set,bc_set)

    end subroutine update_lift
    !************************************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/7/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_models(self,worker,equation_set,bc_set)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bcset_t),              intent(inout)   :: bc_set(:)

        integer(ik)                 :: iface, imodel, idomain_l, ielement_l, idepend, &
                                       ndepend, ChiID, BC_ID, BC_face, ielement_c, istate


        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 



        !
        ! Loop through faces and cache internal, external interpolated states
        !
        do iface = 1,NFACES


            ! Update worker face index
            call worker%set_face(iface)


            ! Update models for face interior. Only depends on interior element.
            idepend = 1
            worker%interpolation_source = 'face interior'
            do imodel = 1,equation_set(idomain_l)%nmodels()

                    worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,DIAG)
                    worker%function_info%idepend = idepend
                    worker%function_info%idiff   = DIAG

                    call equation_set(idomain_l)%models(imodel)%model%compute(worker)
            end do



            ! 
            ! Compute the number of exterior element dependencies for face exterior state
            !
            if ( worker%face_type() == INTERIOR ) then
                ndepend = 1
                
            else if ( worker%face_type() == CHIMERA ) then
                ChiID   = worker%mesh(idomain_l)%faces(ielement_l,iface)%ChiID
                ndepend = worker%mesh(idomain_l)%chimera%recv%data(ChiID)%ndonors()

            else if ( worker%face_type() == BOUNDARY ) then
                BC_ID   = worker%mesh(idomain_l)%faces(ielement_l,iface)%BC_ID
                BC_face = worker%mesh(idomain_l)%faces(ielement_l,iface)%BC_face
                ndepend = bc_set(idomain_l)%bcs(BC_ID)%get_ncoupled_elems(BC_face)

            end if




            !
            ! Face exterior state: interior and chimera
            !
            worker%interpolation_source = 'face exterior'
            if ( (worker%face_type() == INTERIOR) .or. (worker%face_type() == CHIMERA) ) then
                
                do imodel = 1,equation_set(idomain_l)%nmodels()
                    do idepend = 1,ndepend

                        worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,iface)
                        worker%function_info%idepend = idepend

                        call equation_set(idomain_l)%models(imodel)%model%compute(worker)

                    end do !idepend
                end do !imodel




            !
            ! Face exterior state: boundaries
            !
            worker%interpolation_source = 'face exterior'
            else if ( (worker%face_type() == BOUNDARY) ) then

                do imodel = 1,equation_set(idomain_l)%nmodels()
                    do idepend = 1,ndepend

                        ! Get coupled bc element to linearize against.
                        ielement_c = bc_set(idomain_l)%bcs(BC_ID)%bc_patch%coupled_elements(BC_face)%at(idepend)
                        worker%function_info%seed%idomain_g  = worker%mesh(idomain_l)%elems(ielement_c)%idomain_g
                        worker%function_info%seed%idomain_l  = worker%mesh(idomain_l)%elems(ielement_c)%idomain_l
                        worker%function_info%seed%ielement_g = worker%mesh(idomain_l)%elems(ielement_c)%ielement_g
                        worker%function_info%seed%ielement_l = worker%mesh(idomain_l)%elems(ielement_c)%ielement_l
                        worker%function_info%seed%iproc      = IRANK

                        call equation_set(idomain_l)%models(imodel)%model%compute(worker)

                    end do !idepend
                end do !imodel





            end if


        end do !iface






        !
        ! Element volume cache. Models only depend on interior element
        !
        idepend = 1
        worker%interpolation_source = 'element'
        do imodel = 1,equation_set(idomain_l)%nmodels()

                worker%function_info%seed    = element_compute_seed(worker%mesh,idomain_l,ielement_l,idepend,DIAG)
                worker%function_info%idepend = idepend

                call equation_set(idomain_l)%models(imodel)%model%compute(worker)

        end do !imodel


    end subroutine update_models
    !************************************************************************************************































    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------
    subroutine update_lift_faces_internal(self,worker,equation_set,bc_set)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bcset_t),              intent(inout)   :: bc_set(:)

        character(:),   allocatable :: field
        integer(ik)                 :: idomain_l, ielement_l, iface, idepend, &
                                       ndepend, BC_ID, BC_face, ChiID, ieqn

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



            associate ( weights          => worker%mesh(idomain_l)%elems(ielement_l)%gq%face%weights(:,iface),      &
                        val_face_trans   => worker%mesh(idomain_l)%elems(ielement_l)%gq%face%val_trans(:,:,iface),  &
                        val_face         => worker%mesh(idomain_l)%elems(ielement_l)%gq%face%val(:,:,iface),        &
                        val_vol          => worker%mesh(idomain_l)%elems(ielement_l)%gq%vol%val,                    &
                        invmass          => worker%mesh(idomain_l)%elems(ielement_l)%invmass,                       &
                        br2_face         => worker%mesh(idomain_l)%faces(ielement_l,iface)%br2_face,                &
                        br2_vol          => worker%mesh(idomain_l)%faces(ielement_l,iface)%br2_vol)





            do ieqn = 1,worker%mesh(idomain_l)%neqns

                !
                ! Get field
                !
                field = worker%prop(idomain_l)%get_primary_field_name(ieqn)


                !
                ! Compute Interior lift, differentiated wrt Interior
                !
                ndepend = 1
                do idepend = 1,ndepend

                    ! Get Seed
                    worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,DIAG)
                    worker%function_info%idepend = idepend


                    ! Get interior/exterior state
                    var_m = worker%cache%get_data(field,'face interior', 'value', 0, worker%function_info%seed, iface)
                    var_p = worker%cache%get_data(field,'face exterior', 'value', 0, worker%function_info%seed, iface)

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
                    ! Project onto basis
!                    rhs_x = matmul(val_face_trans,var_diff_x)
!                    rhs_y = matmul(val_face_trans,var_diff_y)
!                    rhs_z = matmul(val_face_trans,var_diff_z)
! 
!                    ! Local solve for lift modes in element basis
!                    lift_modes_x = matmul(invmass,rhs_x)
!                    lift_modes_y = matmul(invmass,rhs_y)
!                    lift_modes_z = matmul(invmass,rhs_z)
! 
!                    ! Evaluate lift modes at face quadrature nodes
!                    lift_gq_face_x = matmul(val_face,lift_modes_x)
!                    lift_gq_face_y = matmul(val_face,lift_modes_y)
!                    lift_gq_face_z = matmul(val_face,lift_modes_z)

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


!                    ! Evaluate lift modes at volume quadrature nodes
!                    lift_gq_vol_x = matmul(val_vol,lift_modes_x)
!                    lift_gq_vol_y = matmul(val_vol,lift_modes_y)
!                    lift_gq_vol_z = matmul(val_vol,lift_modes_z)

                    lift_gq_vol_x = matmul(br2_vol,var_diff_x)
                    lift_gq_vol_y = matmul(br2_vol,var_diff_y)
                    lift_gq_vol_z = matmul(br2_vol,var_diff_z)

                    
                    ! Store lift
                    call worker%cache%set_data(field,'face interior', lift_gq_vol_x, 'lift element', 1, worker%function_info%seed, iface)
                    call worker%cache%set_data(field,'face interior', lift_gq_vol_y, 'lift element', 2, worker%function_info%seed, iface)
                    call worker%cache%set_data(field,'face interior', lift_gq_vol_z, 'lift element', 3, worker%function_info%seed, iface)



                end do !idepend





                ! 
                ! Compute the number of dependencies from external state. Sets 'ndepend'
                !
                if ( worker%face_type() == INTERIOR ) then
                    ndepend = 1
                    
                else if ( worker%face_type() == CHIMERA ) then
                    ChiID   = worker%mesh(idomain_l)%faces(ielement_l,iface)%ChiID
                    ndepend = worker%mesh(idomain_l)%chimera%recv%data(ChiID)%ndonors()

                else if ( worker%face_type() == BOUNDARY ) then
                    BC_ID   = worker%mesh(idomain_l)%faces(ielement_l,iface)%BC_ID
                    BC_face = worker%mesh(idomain_l)%faces(ielement_l,iface)%BC_face
                    ndepend = bc_set(idomain_l)%bcs(BC_ID)%get_ncoupled_elems(BC_face)

                end if





                !
                ! Compute Interior lift, differentiated wrt Exterior
                !
                do idepend = 1,ndepend

                    ! Get Seed
                    worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,iface)
                    worker%function_info%idepend = idepend


                    ! Get interior/exterior state
                    var_m = worker%cache%get_data(field,'face interior', 'value', 0, worker%function_info%seed, iface)
                    var_p = worker%cache%get_data(field,'face exterior', 'value', 0, worker%function_info%seed, iface)

                    ! Difference
                    var_diff = HALF*(var_p - var_m) 

                    ! Multiply by weights
                    var_diff_weighted = var_diff * weights

                    ! Multiply by normal. Note: normal is scaled by face jacobian.
                    var_diff_x = var_diff_weighted * worker%normal(1)
                    var_diff_y = var_diff_weighted * worker%normal(2)
                    var_diff_z = var_diff_weighted * worker%normal(3)

!                    ! Project onto basis
!                    rhs_x = matmul(val_face_trans,var_diff_x)
!                    rhs_y = matmul(val_face_trans,var_diff_y)
!                    rhs_z = matmul(val_face_trans,var_diff_z)
!
!                    ! Local solve for lift modes in element basis
!                    lift_modes_x = matmul(invmass,rhs_x)
!                    lift_modes_y = matmul(invmass,rhs_y)
!                    lift_modes_z = matmul(invmass,rhs_z)


                    ! Evaluate lift modes at face quadrature nodes
!                    lift_gq_face_x = matmul(val_face,lift_modes_x)
!                    lift_gq_face_y = matmul(val_face,lift_modes_y)
!                    lift_gq_face_z = matmul(val_face,lift_modes_z)
                    lift_gq_face_x = matmul(br2_face,var_diff_x)
                    lift_gq_face_y = matmul(br2_face,var_diff_y)
                    lift_gq_face_z = matmul(br2_face,var_diff_z)
                    
                    ! Store lift
                    call worker%cache%set_data(field,'face interior', lift_gq_face_x, 'lift face', 1, worker%function_info%seed, iface)
                    call worker%cache%set_data(field,'face interior', lift_gq_face_y, 'lift face', 2, worker%function_info%seed, iface)
                    call worker%cache%set_data(field,'face interior', lift_gq_face_z, 'lift face', 3, worker%function_info%seed, iface)


                    ! Evaluate lift modes at volume quadrature nodes
!                    lift_gq_vol_x = matmul(val_vol,lift_modes_x)
!                    lift_gq_vol_y = matmul(val_vol,lift_modes_y)
!                    lift_gq_vol_z = matmul(val_vol,lift_modes_z)
                    lift_gq_vol_x = matmul(br2_vol,var_diff_x)
                    lift_gq_vol_y = matmul(br2_vol,var_diff_y)
                    lift_gq_vol_z = matmul(br2_vol,var_diff_z)
                    
                    ! Store lift
                    call worker%cache%set_data(field,'face interior', lift_gq_vol_x, 'lift element', 1, worker%function_info%seed, iface)
                    call worker%cache%set_data(field,'face interior', lift_gq_vol_y, 'lift element', 2, worker%function_info%seed, iface)
                    call worker%cache%set_data(field,'face interior', lift_gq_vol_z, 'lift element', 3, worker%function_info%seed, iface)




                end do !idepend

            end do !ieqn


            end associate

        end do !iface


    end subroutine update_lift_faces_internal
    !*************************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------
    subroutine update_lift_faces_external(self,worker,equation_set,bc_set)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bcset_t),              intent(inout)   :: bc_set(:)

        integer(ik) :: idomain_l, ielement_l, iface, idepend, ieqn, ndepend, BC_ID, BC_face, ChiID
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
            do ieqn = 1,worker%mesh(idomain_l)%neqns


                !
                ! Compute External lift, differentiated wrt Interior
                !
                ndepend = 1
                do idepend = 1,ndepend

                    ! Get Seed
                    worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,DIAG)
                    worker%function_info%idepend = idepend


                    if (interior_face) then
                        call handle_external_lift__interior_face(worker,equation_set,bc_set,ieqn)
                    else if (boundary_face) then
                        call handle_external_lift__boundary_face(worker,equation_set,bc_set,ieqn)
                    else if (chimera_face) then
                        call handle_external_lift__chimera_face( worker,equation_set,bc_set,ieqn)
                    else
                        call chidg_signal(FATAL,"update_lift_faces_external: unsupported face type")
                    end if


                end do !idepend





                ! 
                ! Compute the number of dependencies from external state. Sets 'ndepend'
                !
                if ( worker%face_type() == INTERIOR ) then
                    ndepend = 1
                    
                else if ( worker%face_type() == CHIMERA ) then
                    ChiID   = worker%mesh(idomain_l)%faces(ielement_l,iface)%ChiID
                    ndepend = worker%mesh(idomain_l)%chimera%recv%data(ChiID)%ndonors()

                else if ( worker%face_type() == BOUNDARY ) then
                    BC_ID   = worker%mesh(idomain_l)%faces(ielement_l,iface)%BC_ID
                    BC_face = worker%mesh(idomain_l)%faces(ielement_l,iface)%BC_face
                    ndepend = bc_set(idomain_l)%bcs(BC_ID)%get_ncoupled_elems(BC_face)

                end if





                !
                ! Compute External lift, differentiated wrt Exterior
                !
                do idepend = 1,ndepend

                    ! Get Seed
                    worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,iface)
                    worker%function_info%idepend = idepend

                    if (interior_face) then
                        call handle_external_lift__interior_face(worker,equation_set,bc_set,ieqn)
                    else if (boundary_face) then
                        call handle_external_lift__boundary_face(worker,equation_set,bc_set,ieqn)
                    else if (chimera_face) then
                        call handle_external_lift__chimera_face( worker,equation_set,bc_set,ieqn)
                    else
                        call chidg_signal(FATAL,"update_lift_faces_external: unsupported face type")
                    end if


                end do !idepend

            end do !ieqn



        end do !iface


    end subroutine update_lift_faces_external
    !*************************************************************************************************




















    !>  Handle computing lift for an external element, when the face is an interior face.
    !!
    !!  In this case, the external element exists and we can just use its data. This is not the case
    !!  for a boundary condition face, and it is complicated further by a Chimera boundary face.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine handle_external_lift__interior_face(worker,equation_set,bc_set,ieqn)
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bcset_t),              intent(inout)   :: bc_set(:)
        integer(ik),                intent(in)      :: ieqn

        integer(ik) :: idomain_l, ielement_l, iface, idomain_l_n, ielement_l_n, iface_n, iproc_n
        logical     :: boundary_face, interior_face, local_neighbor, remote_neighbor

        type(AD_D), allocatable, dimension(:)   ::          &
            var_m, var_p, var_diff, var_diff_weighted,      &
            var_diff_x,     var_diff_y,     var_diff_z,     &
            rhs_x,          rhs_y,          rhs_z,          &
            lift_modes_x,   lift_modes_y,   lift_modes_z,   &
            lift_gq_face_x, lift_gq_face_y, lift_gq_face_z, &
            lift_gq_vol_x,  lift_gq_vol_y,  lift_gq_vol_z

        character(:),   allocatable                 :: field
        real(rk),       allocatable, dimension(:)   :: normx, normy, normz, weights
        real(rk),       allocatable, dimension(:,:) :: val_face_trans, val_face, val_vol, invmass, br2_face


        !
        ! Interior element
        ! 
        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface


        !
        ! Neighbor element
        !
        idomain_l_n  = worker%mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_l
        ielement_l_n = worker%mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_element_l
        iface_n      = worker%mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_face
        iproc_n      = worker%mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_proc

        local_neighbor  = (iproc_n == IRANK)
        remote_neighbor = (iproc_n /= IRANK)


        !
        ! Get field
        !
        field = worker%prop(idomain_l)%get_primary_field_name(ieqn)


        if ( local_neighbor ) then
            weights          = worker%mesh(idomain_l_n)%elems(ielement_l_n)%gq%face%weights(:,iface_n)
            val_face_trans   = worker%mesh(idomain_l_n)%elems(ielement_l_n)%gq%face%val_trans(:,:,iface_n)
            val_face         = worker%mesh(idomain_l_n)%elems(ielement_l_n)%gq%face%val(:,:,iface_n)
            val_vol          = worker%mesh(idomain_l_n)%elems(ielement_l_n)%gq%vol%val
            invmass          = worker%mesh(idomain_l_n)%elems(ielement_l_n)%invmass
            br2_face         = worker%mesh(idomain_l_n)%faces(ielement_l_n,iface_n)%br2_face


        else if ( remote_neighbor ) then
            ! User local element gq instance. Assumes same order of accuracy.
            weights          = worker%mesh(idomain_l)%elems(ielement_l)%gq%face%weights(:,iface_n)
            val_face_trans   = worker%mesh(idomain_l)%elems(ielement_l)%gq%face%val_trans(:,:,iface_n)
            val_face         = worker%mesh(idomain_l)%elems(ielement_l)%gq%face%val(:,:,iface_n)
            val_vol          = worker%mesh(idomain_l)%elems(ielement_l)%gq%vol%val
            invmass          = worker%mesh(idomain_l)%faces(ielement_l,iface)%neighbor_invmass
            br2_face         = worker%mesh(idomain_l)%faces(ielement_l,iface)%neighbor_br2_face


        end if



            ! Use reverse of interior element's normal vector
            normx = -worker%mesh(idomain_l)%faces(ielement_l,iface)%norm(:,1)
            normy = -worker%mesh(idomain_l)%faces(ielement_l,iface)%norm(:,2)
            normz = -worker%mesh(idomain_l)%faces(ielement_l,iface)%norm(:,3)

            ! Get interior/exterior state
            var_m = worker%cache%get_data(field,'face interior', 'value', 0, worker%function_info%seed, iface)
            var_p = worker%cache%get_data(field,'face exterior', 'value', 0, worker%function_info%seed, iface)

            ! Difference. Relative to exterior element, so reversed
            var_diff = HALF*(var_m - var_p) 

            ! Multiply by weights
            var_diff_weighted = var_diff * weights

            ! Multiply by normal. Note: normal is scaled by face jacobian.
            var_diff_x = var_diff_weighted * normx
            var_diff_y = var_diff_weighted * normy
            var_diff_z = var_diff_weighted * normz

            ! Project onto basis
!            rhs_x = matmul(val_face_trans,var_diff_x)
!            rhs_y = matmul(val_face_trans,var_diff_y)
!            rhs_z = matmul(val_face_trans,var_diff_z)
!
!            ! Local solve for lift modes in element basis
!            lift_modes_x = matmul(invmass,rhs_x)
!            lift_modes_y = matmul(invmass,rhs_y)
!            lift_modes_z = matmul(invmass,rhs_z)
!
!            ! Evaluate lift modes at quadrature nodes
!            lift_gq_face_x = matmul(val_face,lift_modes_x)
!            lift_gq_face_y = matmul(val_face,lift_modes_y)
!            lift_gq_face_z = matmul(val_face,lift_modes_z)
            lift_gq_face_x = matmul(br2_face,var_diff_x)
            lift_gq_face_y = matmul(br2_face,var_diff_y)
            lift_gq_face_z = matmul(br2_face,var_diff_z)
            
            ! Store lift
            call worker%cache%set_data(field,'face exterior', lift_gq_face_x, 'lift face', 1, worker%function_info%seed, iface)
            call worker%cache%set_data(field,'face exterior', lift_gq_face_y, 'lift face', 2, worker%function_info%seed, iface)
            call worker%cache%set_data(field,'face exterior', lift_gq_face_z, 'lift face', 3, worker%function_info%seed, iface)

!            ! Evaluate lift modes at quadrature nodes
!            lift_gq_vol_x = matmul(val_vol,lift_modes_x)
!            lift_gq_vol_y = matmul(val_vol,lift_modes_y)
!            lift_gq_vol_z = matmul(val_vol,lift_modes_z)
!            
!            ! Store lift
!            call worker%cache%set_data('face exterior', lift_gq_vol_x, 'lift element', 1, worker%function_info%seed, ieqn, iface)
!            call worker%cache%set_data('face exterior', lift_gq_vol_y, 'lift element', 2, worker%function_info%seed, ieqn, iface)
!            call worker%cache%set_data('face exterior', lift_gq_vol_z, 'lift element', 3, worker%function_info%seed, ieqn, iface)







    end subroutine handle_external_lift__interior_face
    !*************************************************************************************************














    !>  Handle computing lift for an external element, when the face is a boundary face.
    !!
    !!  In this case, the external element does NOT exist, so we use the interior element. This is kind of like
    !!  Assuming that a boundary element exists of equal size to the interior element.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine handle_external_lift__boundary_face(worker,equation_set,bc_set,ieqn)
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bcset_t),              intent(inout)   :: bc_set(:)
        integer(ik),                intent(in)      :: ieqn

        integer(ik) :: idomain_l, ielement_l, iface, idomain_l_n, ielement_l_n, iface_n
        logical     :: boundary_face, interior_face

        type(AD_D), allocatable, dimension(:)   ::          &
            var_m, var_p, var_diff, var_diff_weighted,      &
            var_diff_x,     var_diff_y,     var_diff_z,     &
            rhs_x,          rhs_y,          rhs_z,          &
            lift_modes_x,   lift_modes_y,   lift_modes_z,   &
            lift_gq_x,      lift_gq_y,      lift_gq_z

        character(:),   allocatable                 :: field
        real(rk),       allocatable, dimension(:)   :: normx, normy, normz


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
        field = worker%prop(idomain_l)%get_primary_field_name(ieqn)



        !
        ! Neighbor element
        !
        idomain_l_n  = worker%mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_l
        ielement_l_n = worker%mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_element_l


        associate ( weights          => worker%mesh(idomain_l)%elems(ielement_l)%gq%face%weights(:,iface_n),        &
                    val_face_trans   => worker%mesh(idomain_l)%elems(ielement_l)%gq%face%val_trans(:,:,iface_n),    &
                    val_face         => worker%mesh(idomain_l)%elems(ielement_l)%gq%face%val(:,:,iface_n),          &
                    invmass          => worker%mesh(idomain_l)%elems(ielement_l)%invmass,                           &
                    br2_face         => worker%mesh(idomain_l)%faces(ielement_l,iface)%br2_face)

            ! Get normal vector. Use reverse of the normal vector from the interior element since no exterior element exists.
            normx = -worker%mesh(idomain_l)%faces(ielement_l,iface)%norm(:,1)
            normy = -worker%mesh(idomain_l)%faces(ielement_l,iface)%norm(:,2)
            normz = -worker%mesh(idomain_l)%faces(ielement_l,iface)%norm(:,3)

            ! Get interior/exterior state
            var_m = worker%cache%get_data(field,'face interior', 'value', 0, worker%function_info%seed, iface)
            var_p = worker%cache%get_data(field,'face exterior', 'value', 0, worker%function_info%seed, iface)


            ! Difference. Relative to exterior element, so reversed
            var_diff = HALF*(var_m - var_p) 


            ! Multiply by weights
            var_diff_weighted = var_diff * weights


            ! Multiply by normal. Note: normal is scaled by face jacobian.
            var_diff_x = var_diff_weighted * normx
            var_diff_y = var_diff_weighted * normy
            var_diff_z = var_diff_weighted * normz


            ! Project onto basis
!            rhs_x = matmul(val_face_trans,var_diff_x)
!            rhs_y = matmul(val_face_trans,var_diff_y)
!            rhs_z = matmul(val_face_trans,var_diff_z)
!
!
!            ! Local solve for lift modes in element basis
!            lift_modes_x = matmul(invmass,rhs_x)
!            lift_modes_y = matmul(invmass,rhs_y)
!            lift_modes_z = matmul(invmass,rhs_z)
!
!
!            ! Evaluate lift modes at quadrature nodes
!            lift_gq_x = matmul(val_face,lift_modes_x)
!            lift_gq_y = matmul(val_face,lift_modes_y)
!            lift_gq_z = matmul(val_face,lift_modes_z)
            lift_gq_x = matmul(br2_face,var_diff_x)
            lift_gq_y = matmul(br2_face,var_diff_y)
            lift_gq_z = matmul(br2_face,var_diff_z)
            

            ! Store lift
            call worker%cache%set_data(field,'face exterior', lift_gq_x, 'lift face', 1, worker%function_info%seed, iface)
            call worker%cache%set_data(field,'face exterior', lift_gq_y, 'lift face', 2, worker%function_info%seed, iface)
            call worker%cache%set_data(field,'face exterior', lift_gq_z, 'lift face', 3, worker%function_info%seed, iface)


        end associate






    end subroutine handle_external_lift__boundary_face
    !*************************************************************************************************












    !>  Handle computing lift for an external element, when the face is an interior face.
    !!
    !!  In this case, the external element exists and we can just use its data. This is not the case
    !!  for a boundary condition face, and it is complicated further by a Chimera boundary face.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine handle_external_lift__chimera_face(worker,equation_set,bc_set,ieqn)
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bcset_t),              intent(inout)   :: bc_set(:)
        integer(ik),                intent(in)      :: ieqn

        integer(ik) :: idomain_l, ielement_l, iface, idomain_l_n, ielement_l_n, iface_n
        logical     :: boundary_face, interior_face

        type(AD_D), allocatable, dimension(:)   ::          &
            var_m, var_p, var_diff, var_diff_weighted,      &
            var_diff_x,     var_diff_y,     var_diff_z,     &
            rhs_x,          rhs_y,          rhs_z,          &
            lift_modes_x,   lift_modes_y,   lift_modes_z,   &
            lift_gq_face_x, lift_gq_face_y, lift_gq_face_z, &
            lift_gq_vol_x,  lift_gq_vol_y,  lift_gq_vol_z

        character(:),   allocatable                 :: field
        real(rk),       allocatable, dimension(:)   :: normx, normy, normz


        !
        ! Interior element
        ! 
        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 
        iface      = worker%iface


        !
        ! Get field
        !
        field = worker%prop(idomain_l)%get_primary_field_name(ieqn)

        !
        ! Use components from receiver element since no single element exists to act as the exterior element.
        ! This implicitly treats the diffusion terms as if there were a reflected element like the receiver
        ! element that was acting as the donor.
        !
        associate ( weights          => worker%mesh(idomain_l)%elems(ielement_l)%gq%face%weights(:,iface),        &
                    val_face_trans   => worker%mesh(idomain_l)%elems(ielement_l)%gq%face%val_trans(:,:,iface),    &
                    val_face         => worker%mesh(idomain_l)%elems(ielement_l)%gq%face%val(:,:,iface),          &
                    val_vol          => worker%mesh(idomain_l)%elems(ielement_l)%gq%vol%val,                      &
                    invmass          => worker%mesh(idomain_l)%elems(ielement_l)%invmass,                         &
                    br2_face         => worker%mesh(idomain_l)%faces(ielement_l,iface)%br2_face)

            ! Use reversed normal vectors of receiver element
            normx = -worker%mesh(idomain_l)%faces(ielement_l,iface)%norm(:,1)
            normy = -worker%mesh(idomain_l)%faces(ielement_l,iface)%norm(:,2)
            normz = -worker%mesh(idomain_l)%faces(ielement_l,iface)%norm(:,3)

            ! Get interior/exterior state
            var_m = worker%cache%get_data(field,'face interior', 'value', 0, worker%function_info%seed, iface)
            var_p = worker%cache%get_data(field,'face exterior', 'value', 0, worker%function_info%seed, iface)

            ! Difference. Relative to exterior element, so reversed
            var_diff = HALF*(var_m - var_p) 

            ! Multiply by weights
            var_diff_weighted = var_diff * weights

            ! Multiply by normal. Note: normal is scaled by face jacobian.
            var_diff_x = var_diff_weighted * normx
            var_diff_y = var_diff_weighted * normy
            var_diff_z = var_diff_weighted * normz

            ! Project onto basis
!            rhs_x = matmul(val_face_trans,var_diff_x)
!            rhs_y = matmul(val_face_trans,var_diff_y)
!            rhs_z = matmul(val_face_trans,var_diff_z)
!
!            ! Local solve for lift modes in element basis
!            lift_modes_x = matmul(invmass,rhs_x)
!            lift_modes_y = matmul(invmass,rhs_y)
!            lift_modes_z = matmul(invmass,rhs_z)
!
!            ! Evaluate lift modes at quadrature nodes
!            lift_gq_face_x = matmul(val_face,lift_modes_x)
!            lift_gq_face_y = matmul(val_face,lift_modes_y)
!            lift_gq_face_z = matmul(val_face,lift_modes_z)
            lift_gq_face_x = matmul(br2_face,var_diff_x)
            lift_gq_face_y = matmul(br2_face,var_diff_y)
            lift_gq_face_z = matmul(br2_face,var_diff_z)
            
            ! Store lift
            call worker%cache%set_data(field,'face exterior', lift_gq_face_x, 'lift face', 1, worker%function_info%seed, iface)
            call worker%cache%set_data(field,'face exterior', lift_gq_face_y, 'lift face', 2, worker%function_info%seed, iface)
            call worker%cache%set_data(field,'face exterior', lift_gq_face_z, 'lift face', 3, worker%function_info%seed, iface)

!            ! Evaluate lift modes at quadrature nodes
!            lift_gq_vol_x = matmul(val_vol,lift_modes_x)
!            lift_gq_vol_y = matmul(val_vol,lift_modes_y)
!            lift_gq_vol_z = matmul(val_vol,lift_modes_z)
!            
!            ! Store lift
!            call worker%cache%set_data('face exterior', lift_gq_vol_x, 'lift element', 1, worker%function_info%seed, ieqn, iface)
!            call worker%cache%set_data('face exterior', lift_gq_vol_y, 'lift element', 2, worker%function_info%seed, ieqn, iface)
!            call worker%cache%set_data('face exterior', lift_gq_vol_z, 'lift element', 3, worker%function_info%seed, ieqn, iface)

        end associate






    end subroutine handle_external_lift__chimera_face
    !*************************************************************************************************






































end module type_cache_handler
