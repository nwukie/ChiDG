module type_cache_handler
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: NFACES, INTERIOR, CHIMERA, BOUNDARY, DIAG, ME, NEIGHBOR
    use mod_DNAD_tools,     only: face_compute_seed, element_compute_seed
    use mod_interpolate,    only: interpolate_face_autodiff, interpolate_element_autodiff
    use DNAD_D

    use type_chidg_cache,   only: chidg_cache_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_equation_set,  only: equation_set_t
    use type_bcset,         only: bcset_t
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
    !---------------------------------------------------------------------------------------
    type, public :: cache_handler_t



    contains

        procedure   :: update

    end type cache_handler_t
    !***************************************************************************************





contains






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/7/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine update(self,worker,equation_set,bc_set)
        class(cache_handler_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(equation_set_t),       intent(inout)   :: equation_set(:)
        type(bcset_t),              intent(inout)   :: bc_set(:)

        integer(ik)         :: iface, iside, ieqn, idomain_l, ielement_l, idepend, ndepend, ChiID

        type(AD_D), allocatable, dimension(:) :: value_gq



        idomain_l  = worker%element_info%idomain_l 
        ielement_l = worker%element_info%ielement_l 


        !
        ! Resize cache
        !
        call worker%cache%resize(worker%mesh,idomain_l,ielement_l)


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
                call worker%cache%set_data('face interior',value_gq,'value',0,idepend,worker%function_info%seed,ieqn,iface)

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
!                BCID = worker%mesh(idomain_l)%faces(ielemen_l,iface)%BCID
!                ndepend = bc_set%bc(BCID)%ndepend()
            end if




            !
            ! Face exterior state
            !
            if ( (worker%face_type() == INTERIOR) .or. (worker%face_type() == CHIMERA) ) then
                do ieqn = 1,worker%mesh(idomain_l)%neqns
                    do idepend = 1,ndepend

                        worker%function_info%seed    = face_compute_seed(worker%mesh,idomain_l,ielement_l,iface,idepend,iface)
                        worker%function_info%idepend = idepend

                        value_gq = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,worker%face_info(),worker%function_info,ieqn,'value',NEIGHBOR)

                        call worker%cache%set_data('face exterior',value_gq,'value',0,idepend,worker%function_info%seed,ieqn,iface)

                    end do !idepend
                end do !ieqn




            else if ( (worker%face_type() == BOUNDARY) ) then

!                do ieqn = 1,worker%mesh(idomain_l)%neqns
!                    idepend = 1,ndepend
!
!
!                    end do !idepend
!                end do !ieqn




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

                call worker%cache%set_data('element',value_gq,'value',0,idepend,worker%function_info%seed,ieqn)

        end do !ieqn













    end subroutine update
    !***************************************************************************************

















end module type_cache_handler
