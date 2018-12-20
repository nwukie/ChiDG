module mod_modify_solution
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES, DIAG, CHIMERA, INTERIOR, NO_ID
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM, GLOBAL_MASTER
    use mod_io,                 only: verbosity
    use mpi_f08,                only: MPI_Barrier
    use mod_interpolate



    use type_chidg_data,        only: chidg_data_t
    use type_chidg_vector,        only: chidg_vector_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_chidg_cache,       only: chidg_cache_t
    use type_cache_handler,     only: cache_handler_t
    use type_element_info,      only: element_info_t
    use type_timer,             only: timer_t

    use DNAD_D
    implicit none


contains





    !>  Spatial loop through domains, elements, and faces. Functions get called for each 
    !!  element/face.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/15/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!  @note   Improved layout, added computation of diffusion terms.
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine modify_solution(data,q0,timing)
        type(chidg_data_t), intent(inout),  target      :: data
        type(chidg_vector_t), intent(inout)                :: q0 
        real(rk),           intent(inout),  optional    :: timing

        integer(ik)                 :: idom, ielem, iface, idiff, ierr, &
                                       diff_min, diff_max, eqn_ID
        type(timer_t)               :: timer, comm_timer, loop_timer

        logical                  :: differentiate, unreal_elem_11, unreal_elem_22, unreal_elem_33, unreal_face_11, unreal_face_22, unreal_face_33
        logical                  :: unreal_elem_12, unreal_elem_13, unreal_elem_23, unreal_face_12, unreal_face_13, unreal_face_23, unreal_elem_det, unreal_face_det
        logical                  :: unreal_elem, unreal_face
        logical                  :: unreal_elem_old, unreal_face_old, unreal_elem_time_avg, unreal_face_time_avg
        logical                  :: reject_new, reject_time_avg, reject_nbr_avg, reject_old, reject_zero
        type(element_info_t)        :: elem_info
        type(chidg_worker_t)        :: worker
        type(chidg_cache_t)         :: cache
        type(cache_handler_t)       :: cache_handler
        real(rk), allocatable :: r11(:), r22(:), r33(:), r12(:), r13(:), r23(:), r11f(:), r22f(:), r33f(:), r12f(:), r13f(:), r23f(:)
        real(rk), allocatable :: r11_old(:), r22_old(:), r33_old(:), r12_old(:), r13_old(:), r23_old(:), r11f_old(:), r22f_old(:), r33f_old(:), r12f_old(:), r13f_old(:), r23f_old(:)
        real(rk), allocatable :: qtemp(:), qtemp0(:),det_elem(:), det_face(:)
        real(rk), allocatable :: avg_vals(:,:)
        integer(ik) :: ivar, nterms, niter, maxiter
        character(:), allocatable :: eqnset_name

        real(rk) :: shrink, nbr_val, nbr_vol, my_val, ctrl_val
        integer(ik) :: ival, nvals
        real(rk) :: eps
        
    
    
        eps = 1.0e-12_rk
        shrink = 0.95_rk

        differentiate = .false.
        
        !
        ! Initialize Chidg Worker references
        !
        call worker%init(data%mesh, data%eqnset(:)%prop, data%sdata, data%time_manager, cache)


        !
        ! Start timer on spatial discretization update
        !
        call timer%start()

    

        !
        ! Clear function_status data. This tracks if a function has already been called. So, 
        ! in this way we can compute a function on a face and apply it to both elements. 
        ! The function is just registered as computed for both. So we need to reset all of 
        ! that data here. This is only tracked for the interior scheme. Boundary condition 
        ! evaluations and Chimera faces are not tracked.
        !
        call data%sdata%function_status%clear()




        !
        ! Communicate solution vector
        !
        call comm_timer%start()
        call data%sdata%q%comm_send()
        call data%sdata%q%comm_recv()
        call data%sdata%q%comm_wait()
        call comm_timer%stop()




        !
        ! Loop through given element compute the residual functions and also the 
        ! linearization of those functions.
        !
        call write_line('-  Enforcing realizability...', io_proc=GLOBAL_MASTER, silence=(verbosity < 3))


        !
        ! Clear function_status data. This tracks if a function has already been called. So, in this way
        ! we can compute a function on a face and apply it to both elements. The function is just registered
        ! as computed for both. So we need to reset all of that data here. This is only tracked for the interior scheme.
        ! Boundary condition evaluations and Chimera faces are not tracked.
        !
        call data%sdata%function_status%clear()



        !
        ! Set time info on chidg_worker
        !
        worker%itime = data%time_manager%itime
        worker%t     = data%time_manager%t





        call loop_timer%start()
        
        ! STEP 1: Check nonnegativity of diagonal RS components
        ! If any nodes have negative values,
        ! we drop to the constant mode.
        ! If the constant mode value is negative, it is reset to zero.

                    
        do idom = 1,data%mesh%ndomains()
            associate ( domain => data%mesh%domain(idom))

            ! Loop through elements in the current domain
            do ielem = 1,domain%nelem
                eqn_ID = worker%mesh%domain(idom)%elems(ielem)%eqn_ID


                elem_info%idomain_g  = domain%elems(ielem)%idomain_g
                elem_info%idomain_l  = domain%elems(ielem)%idomain_l
                elem_info%ielement_g = domain%elems(ielem)%ielement_g
                elem_info%ielement_l = domain%elems(ielem)%ielement_l
                call worker%set_element(elem_info)


                ! Update the element cache
                call cache_handler%update(worker,data%eqnset, data%bc_state_group, components    = 'none',           &
                                                                                   face          = NO_ID,           &
                                                                                   differentiate = differentiate,   &
                                                                                   lift          = .false.)
                !call eqnset%modify_solution(worker, differentiate)
                eqnset_name = data%get_equation_set_name(eqn_ID)
                if (trim(eqnset_name)=='RANS_RSTM') then
!                    ! Omega POSITIVITY
!                    eps = 1.0e-12_rk
!                    unreal_elem = .true.
!                    unreal_face = .true.
!
!                    reject_new      = .true.
!                    reject_time_avg = .true.
!                    reject_nbr_avg  = .true.
!                    reject_old      = .true.
!                    reject_zero     = .true.
!                    
!                    do while ( reject_new .and. reject_time_avg .and. reject_nbr_avg .and. reject_old .and. reject_zero)
!                        ivar = worker%cache%element%get_field_index('Density * Omega') 
!                        r11 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')
!                        r11_old = interpolate_element_standard(worker%mesh,q0,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')
!
!                        unreal_elem     = (any(r11 < eps))
!                        unreal_elem_old = (any(r11_old < eps))
!                        unreal_elem_time_avg = (any(0.5_rk*(r11+r11_old)<eps))
!                        
!
!                        unreal_face = .false.
!                        unreal_face_old = .false.
!                        unreal_face_time_avg = .false.
!                        do iface = 1,NFACES
!
!                            call worker%set_face(iface)
!
!                            if (allocated(r11f)) deallocate(r11f)
!                            if (allocated(r11f_old)) deallocate(r11f_old)
!
!                            ivar = worker%cache%element%get_field_index('Density * Omega') 
!                            r11f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)
!                            r11f_old = interpolate_face_standard(worker%mesh,q0,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)
!
!
!                            if (any(r11f < eps))                     unreal_face     = .true.
!                            if (any(r11f_old < eps))                 unreal_face_old = .true.
!                            if (any(0.5_rk*(r11f+r11f_old) < eps))   unreal_face_time_avg = .true.
!
!                        end do  ! faces loop
!
!                        if (.not. (unreal_elem .or. unreal_face)) then
!                            ! Accept the new solution without modification
!                            reject_new = .false.
!
!                        else if (.not. (unreal_elem_time_avg .or. unreal_face_time_avg)) then
!                            ! Accept the time average of the old and new solution
!                            reject_time_avg = .false.
!                            qtemp0 = q0%dom(idom)%vecs(ielem)%getvar(ivar, 1)
!
!                            qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)
!
!                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, 0.5_rk*(qtemp+qtemp0))
!
!                        !else if (.not. (unreal_elem_old .or. unreal_face_old)) then
!                        !    ! Accept the old solution
!                        !    reject_old = .false.
!                        !    qtemp0 = q0%dom(idom)%vecs(ielem)%getvar(ivar, 1)
!
!
!                        !    call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp0)
!
!
!                        else 
!                            ! Otherwise, get the average values from the neighbor
!                            ! Take the average value of positive values
!                            ! If no neighbors with positive average values, just set it to zero
!                            reject_zero = .false.
!
!                            if (allocated(avg_vals)) deallocate(avg_vals)
!                            avg_vals = get_neighbor_average_values(worker, worker%mesh, data%sdata%q, ivar, 1)
!                            nvals = size(avg_vals,1)
!                            nbr_val = eps
!                            nbr_vol = 0.0_rk
!                            do ival = 1, nvals
!                                if (eps<avg_vals(ival,1)) then
!                                    nbr_val = nbr_val + avg_vals(ival,1)*avg_vals(ival,2)
!                                    nbr_vol = nbr_vol + avg_vals(ival,2)
!                                end if
!                            end do
!
!                            if (nbr_vol > eps) then
!                                my_val = nbr_val/nbr_vol
!                            else
!                                my_val = eps
!                            end if
!
!                            qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)
!                            qtemp = 0.0_rk*qtemp
!                            qtemp(1) = my_val
!
!                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp)
!
!                        end if
!
!                    end do

                    eps = 0.0_rk

                    ! R11 POSITIVITY
                    unreal_elem = .true.
                    unreal_face = .true.

                    reject_new      = .true.
                    reject_time_avg = .true.
                    reject_nbr_avg  = .true.
                    reject_old      = .true.
                    reject_zero     = .true.
                    
                    do while ( reject_new .and. reject_time_avg .and. reject_nbr_avg .and. reject_old .and. reject_zero)
                        ivar = worker%cache%element%get_field_index('Density * Reynolds-11') 
                        r11 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')
                        r11_old = interpolate_element_standard(worker%mesh,q0,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                        unreal_elem     = (any(r11 < eps))
                        unreal_elem_old = (any(r11_old < eps))
                        unreal_elem_time_avg = (any(0.5_rk*(r11+r11_old)<eps))
                        

                        unreal_face = .false.
                        unreal_face_old = .false.
                        unreal_face_time_avg = .false.
                        do iface = 1,NFACES

                            call worker%set_face(iface)

                            if (allocated(r11f)) deallocate(r11f)
                            if (allocated(r11f_old)) deallocate(r11f_old)

                            ivar = worker%cache%element%get_field_index('Density * Reynolds-11') 
                            r11f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)
                            r11f_old = interpolate_face_standard(worker%mesh,q0,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)


                            if (any(r11f < eps))                     unreal_face     = .true.
                            if (any(r11f_old < eps))                 unreal_face_old = .true.
                            if (any(0.5_rk*(r11f+r11f_old) < eps))   unreal_face_time_avg = .true.

                        end do  ! faces loop

                        if (.not. (unreal_elem .or. unreal_face)) then
                            ! Accept the new solution without modification
                            reject_new = .false.

                        else if (.not. (unreal_elem_time_avg .or. unreal_face_time_avg)) then
                            ! Accept the time average of the old and new solution
                            reject_time_avg = .false.
                            qtemp0 = q0%dom(idom)%vecs(ielem)%getvar(ivar, 1)

                            qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)

                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, 0.5_rk*(qtemp+qtemp0))

                        !else if (.not. (unreal_elem_old .or. unreal_face_old)) then
                        !    ! Accept the old solution
                        !    reject_old = .false.
                        !    qtemp0 = q0%dom(idom)%vecs(ielem)%getvar(ivar, 1)


                        !    call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp0)


                        else 
                            ! Otherwise, get the average values from the neighbor
                            ! Take the average value of positive values
                            ! If no neighbors with positive average values, just set it to zero
                            reject_zero = .false.

                            if (allocated(avg_vals)) deallocate(avg_vals)
                            avg_vals = get_neighbor_average_values(worker, worker%mesh, data%sdata%q, ivar, 1)
                            nvals = size(avg_vals,1)
                            nbr_val = eps
                            nbr_vol = 0.0_rk
                            do ival = 1, nvals
                                if (eps<avg_vals(ival,1)) then
                                    nbr_val = nbr_val + avg_vals(ival,1)*avg_vals(ival,2)
                                    nbr_vol = nbr_vol + avg_vals(ival,2)
                                end if
                            end do

                            if (nbr_vol > eps) then
                                my_val = nbr_val/nbr_vol
                            else
                                my_val = eps
                            end if

                            qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)
                            qtemp = 0.0_rk*qtemp
                            qtemp(1) = my_val

                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp)

                        end if

                    end do
 
                    ! R22 POSITIVITY
                    unreal_elem = .true.
                    unreal_face = .true.

                    reject_new      = .true.
                    reject_time_avg = .true.
                    reject_nbr_avg  = .true.
                    reject_old      = .true.
                    reject_zero     = .true.
                    
                    do while ( reject_new .and. reject_time_avg .and. reject_nbr_avg .and. reject_old .and. reject_zero)
                        ivar = worker%cache%element%get_field_index('Density * Reynolds-22') 
                        r22 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')
                        r22_old = interpolate_element_standard(worker%mesh,q0,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                        unreal_elem     = (any(r22 < eps))
                        unreal_elem_old = (any(r22_old < eps))
                        unreal_elem_time_avg = (any(0.5_rk*(r22+r22_old)<eps))
                        

                        unreal_face = .false.
                        unreal_face_old = .false.
                        unreal_face_time_avg = .false.
                        do iface = 1,NFACES

                            call worker%set_face(iface)

                            if (allocated(r22f)) deallocate(r22f)
                            if (allocated(r22f_old)) deallocate(r22f_old)

                            ivar = worker%cache%element%get_field_index('Density * Reynolds-22') 
                            r22f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)
                            r22f_old = interpolate_face_standard(worker%mesh,q0,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)


                            if (any(r22f < eps))                     unreal_face     = .true.
                            if (any(r22f_old < eps))                 unreal_face_old = .true.
                            if (any(0.5_rk*(r22f+r22f_old) < eps))   unreal_face_time_avg = .true.

                        end do  ! faces loop

                        if (.not. (unreal_elem .or. unreal_face)) then
                            ! Accept the new solution without modification
                            reject_new = .false.

                        else if (.not. (unreal_elem_time_avg .or. unreal_face_time_avg)) then
                            ! Accept the time average of the old and new solution
                            reject_time_avg = .false.
                            qtemp0 = q0%dom(idom)%vecs(ielem)%getvar(ivar, 1)

                            qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)

                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, 0.5_rk*(qtemp+qtemp0))

                        !else if (.not. (unreal_elem_old .or. unreal_face_old)) then
                        !    ! Accept the old solution
                        !    reject_old = .false.
                        !    qtemp0 = q0%dom(idom)%vecs(ielem)%getvar(ivar, 1)


                        !    call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp0)


                        else 
                            ! Otherwise, get the average values from the neighbor
                            ! Take the average value of positive values
                            ! If no neighbors with positive average values, just set it to zero
                            reject_zero = .false.

                            if (allocated(avg_vals)) deallocate(avg_vals)
                            avg_vals = get_neighbor_average_values(worker, worker%mesh, data%sdata%q, ivar, 1)
                            nvals = size(avg_vals,1)
                            nbr_val = eps
                            nbr_vol = 0.0_rk
                            do ival = 1, nvals
                                if (eps<avg_vals(ival,1)) then
                                    nbr_val = nbr_val + avg_vals(ival,1)*avg_vals(ival,2)
                                    nbr_vol = nbr_vol + avg_vals(ival,2)
                                end if
                            end do

                            if (nbr_vol > eps) then
                                my_val = nbr_val/nbr_vol
                            else
                                my_val = eps
                            end if

                            qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)
                            qtemp = 0.0_rk*qtemp
                            qtemp(1) = my_val

                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp)


                        end if

                    end do
 
                    ! R33 POSITIVITY
                    unreal_elem = .true.
                    unreal_face = .true.

                    reject_new      = .true.
                    reject_time_avg = .true.
                    reject_nbr_avg  = .true.
                    reject_old      = .true.
                    reject_zero     = .true.
                    
                    do while ( reject_new .and. reject_time_avg .and. reject_nbr_avg .and. reject_old .and. reject_zero)
                        ivar = worker%cache%element%get_field_index('Density * Reynolds-33') 
                        r33 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')
                        r33_old = interpolate_element_standard(worker%mesh,q0,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                        unreal_elem     = (any(r33 < eps))
                        unreal_elem_old = (any(r33_old < eps))
                        unreal_elem_time_avg = (any(0.5_rk*(r33+r33_old)<eps))
                        

                        unreal_face = .false.
                        unreal_face_old = .false.
                        unreal_face_time_avg = .false.
                        do iface = 1,NFACES

                            call worker%set_face(iface)

                            if (allocated(r33f)) deallocate(r33f)
                            if (allocated(r33f_old)) deallocate(r33f_old)

                            ivar = worker%cache%element%get_field_index('Density * Reynolds-33') 
                            r33f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)
                            r33f_old = interpolate_face_standard(worker%mesh,q0,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)


                            if (any(r33f < eps))                     unreal_face     = .true.
                            if (any(r33f_old < eps))                 unreal_face_old = .true.
                            if (any(0.5_rk*(r33f+r33f_old) < eps))   unreal_face_time_avg = .true.

                        end do  ! faces loop

                        if (.not. (unreal_elem .or. unreal_face)) then
                            ! Accept the new solution without modification
                            reject_new = .false.

                        else if (.not. (unreal_elem_time_avg .or. unreal_face_time_avg)) then
                            ! Accept the time average of the old and new solution
                            reject_time_avg = .false.
                            qtemp0 = q0%dom(idom)%vecs(ielem)%getvar(ivar, 1)

                            qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)

                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, 0.5_rk*(qtemp+qtemp0))

                        !else if (.not. (unreal_elem_old .or. unreal_face_old)) then
                        !    ! Accept the old solution
                        !    reject_old = .false.
                        !    qtemp0 = q0%dom(idom)%vecs(ielem)%getvar(ivar, 1)


                        !    call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp0)


                        else 
                            ! Otherwise, get the average values from the neighbor
                            ! Take the average value of positive values
                            ! If no neighbors with positive average values, just set it to zero
                            reject_zero = .false.

                            if (allocated(avg_vals)) deallocate(avg_vals)
                            avg_vals = get_neighbor_average_values(worker, worker%mesh, data%sdata%q, ivar, 1)
                            nvals = size(avg_vals,1)
                            nbr_val = eps
                            nbr_vol = 0.0_rk
                            do ival = 1, nvals
                                if (eps<avg_vals(ival,1)) then
                                    nbr_val = nbr_val + avg_vals(ival,1)*avg_vals(ival,2)
                                    nbr_vol = nbr_vol + avg_vals(ival,2)
                                end if
                            end do

                            if (nbr_vol > eps) then
                                my_val = nbr_val/nbr_vol
                            else
                                my_val = eps
                            end if

                            qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)
                            qtemp = 0.0_rk*qtemp
                            qtemp(1) = my_val

                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp)


                        end if

                    end do
 
                  
                end if


            end do  ! ielem
            end associate
        end do  ! idom

        ! STEP 2: Check the Cauchy-Schwarz identity
        ! If any off-diagonal RS component violates CS at any node,
        ! we drop to the constant mode and set it to zero.
        ! NOTE: We need to use the modified values of diag(RS) from step 1.
        eps = 0.0_rk
        do idom = 1,data%mesh%ndomains()
            associate ( domain => data%mesh%domain(idom))

            ! Loop through elements in the current domain
            do ielem = 1,domain%nelem
                eqn_ID = worker%mesh%domain(idom)%elems(ielem)%eqn_ID


                elem_info%idomain_g  = domain%elems(ielem)%idomain_g
                elem_info%idomain_l  = domain%elems(ielem)%idomain_l
                elem_info%ielement_g = domain%elems(ielem)%ielement_g
                elem_info%ielement_l = domain%elems(ielem)%ielement_l
                call worker%set_element(elem_info)


                ! Update the element cache
                call cache_handler%update(worker,data%eqnset, data%bc_state_group, components    = 'none',           &
                                                                                   face          = NO_ID,           &
                                                                                   differentiate = differentiate,   &
                                                                                   lift          = .false.)
               

                !call eqnset%modify_solution(worker, differentiate)
                eqnset_name = data%get_equation_set_name(eqn_ID)
                if (trim(eqnset_name)=='RANS_RSTM') then

                    ! R12 CS
                    unreal_elem = .true.
                    unreal_face = .true.

                    reject_new      = .true.
                    reject_time_avg = .true.
                    reject_nbr_avg  = .true.
                    reject_old      = .true.
                    reject_zero     = .true.
                    
                    do while ( reject_new .and. reject_time_avg .and. reject_nbr_avg .and. reject_old .and. reject_zero)
                        ivar = worker%cache%element%get_field_index('Density * Reynolds-11') 
                        r11 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-22') 
                        r22 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-12') 
                        r12 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')
                        r12_old = interpolate_element_standard(worker%mesh,q0,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                        unreal_elem                 = (any(abs(r12) > sqrt(r11*r22)))
                        unreal_elem_old             = (any(abs(r12_old) > sqrt(r11*r22)))
                        unreal_elem_time_avg        = (any(abs(0.5_rk*(r12+r12_old)) > sqrt(r11*r22)))
                        
                        ctrl_val = minval(sqrt(r11*r22))

                        unreal_face = .false.
                        unreal_face_old = .false.
                        unreal_face_time_avg = .false.
                        do iface = 1,NFACES

                            call worker%set_face(iface)

                            if (allocated(r11f)) deallocate(r11f)
                            if (allocated(r11f_old)) deallocate(r11f_old)

                            ivar = worker%cache%element%get_field_index('Density * Reynolds-12') 
                            r12f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)
                            r12f_old = interpolate_face_standard(worker%mesh,q0,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)

                            ivar = worker%cache%element%get_field_index('Density * Reynolds-11') 
                            r11f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)

                            ivar = worker%cache%element%get_field_index('Density * Reynolds-22') 
                            r22f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)


                            if (any(abs(r12f) > sqrt(r11f*r22f)))                   unreal_face = .true.
                            if (any(abs(r12f_old) > sqrt(r11f*r22f)))               unreal_face_old = .true.
                            if (any(abs(0.5_rk*(r12f+r12f_old)) > sqrt(r11f*r22f))) unreal_face_time_avg = .true.


                            ctrl_val = min(ctrl_val, minval(sqrt(r11f*r22f)))
                        end do  ! faces loop

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-12') 
                        if (.not. (unreal_elem .or. unreal_face)) then
                            ! Accept the new solution without modification
                            reject_new = .false.

                        else if (.not. (unreal_elem_time_avg .or. unreal_face_time_avg)) then
                            ! Accept the time average of the old and new solution
                            reject_time_avg = .false.
                            qtemp0 = q0%dom(idom)%vecs(ielem)%getvar(ivar, 1)

                            qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)

                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, 0.5_rk*(qtemp+qtemp0))

                        !else if (.not. (unreal_elem_old .or. unreal_face_old)) then
                        !    ! Accept the old solution
                        !    reject_old = .false.
                        !    qtemp0 = q0%dom(idom)%vecs(ielem)%getvar(ivar, 1)


                        !    call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp0)


                        else 
                            ! Reject all of the above - nothing we can do but set the solution to zero
                            reject_zero = .false.

                            !qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)
                            !qtemp = 0.0_rk*qtemp

                            !call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp)

                            if (allocated(avg_vals)) deallocate(avg_vals)
                            avg_vals = get_neighbor_average_values(worker, worker%mesh, data%sdata%q, ivar, 1)
                            nvals = size(avg_vals,1)
                            nbr_val = 0.0_rk 
                            nbr_vol = 0.0_rk
                            do ival = 1, nvals
                                if (ctrl_val>abs(avg_vals(ival,1))) then
                                    nbr_val = nbr_val + avg_vals(ival,1)*avg_vals(ival,2)
                                    nbr_vol = nbr_vol + avg_vals(ival,2)
                                end if
                            end do

                            if (nbr_vol > eps) then
                                my_val = nbr_val/nbr_vol
                            else
                                my_val = 0.0_rk
                            end if

                            qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)
                            qtemp = 0.0_rk*qtemp
                            qtemp(1) = my_val

                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp)


                        end if

                    end do
 
                    ! R13 CS
                    unreal_elem = .true.
                    unreal_face = .true.

                    reject_new      = .true.
                    reject_time_avg = .true.
                    reject_nbr_avg  = .true.
                    reject_old      = .true.
                    reject_zero     = .true.
                    
                    do while ( reject_new .and. reject_time_avg .and. reject_nbr_avg .and. reject_old .and. reject_zero)
                        ivar = worker%cache%element%get_field_index('Density * Reynolds-11') 
                        r11 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-33') 
                        r33 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-13') 
                        r13 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')
                        r13_old = interpolate_element_standard(worker%mesh,q0,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                        unreal_elem                 = (any(abs(r13) > sqrt(r11*r33)))
                        unreal_elem_old             = (any(abs(r13_old) > sqrt(r11*r33)))
                        unreal_elem_time_avg        = (any(abs(0.5_rk*(r13+r13_old)) > sqrt(r11*r33)))
                        

                        ctrl_val = minval(sqrt(r11*r33))
                        unreal_face = .false.
                        unreal_face_old = .false.
                        unreal_face_time_avg = .false.
                        do iface = 1,NFACES

                            call worker%set_face(iface)

                            if (allocated(r11f)) deallocate(r11f)
                            if (allocated(r11f_old)) deallocate(r11f_old)

                            ivar = worker%cache%element%get_field_index('Density * Reynolds-13') 
                            r13f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)
                            r13f_old = interpolate_face_standard(worker%mesh,q0,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)

                            ivar = worker%cache%element%get_field_index('Density * Reynolds-11') 
                            r11f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)

                            ivar = worker%cache%element%get_field_index('Density * Reynolds-33') 
                            r33f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)


                            if (any(abs(r13f) > sqrt(r11f*r33f)))                   unreal_face = .true.
                            if (any(abs(r13f_old) > sqrt(r11f*r33f)))               unreal_face_old = .true.
                            if (any(abs(0.5_rk*(r13f+r13f_old)) > sqrt(r11f*r33f))) unreal_face_time_avg = .true.

                            ctrl_val = min(ctrl_val, minval(sqrt(r11f*r33f)))

                        end do  ! faces loop

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-13') 
                        if (.not. (unreal_elem .or. unreal_face)) then
                            ! Accept the new solution without modification
                            reject_new = .false.

                        else if (.not. (unreal_elem_time_avg .or. unreal_face_time_avg)) then
                            ! Accept the time average of the old and new solution
                            reject_time_avg = .false.
                            qtemp0 = q0%dom(idom)%vecs(ielem)%getvar(ivar, 1)

                            qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)

                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, 0.5_rk*(qtemp+qtemp0))

                        !else if (.not. (unreal_elem_old .or. unreal_face_old)) then
                        !    ! Accept the old solution
                        !    reject_old = .false.
                        !    qtemp0 = q0%dom(idom)%vecs(ielem)%getvar(ivar, 1)


                        !    call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp0)


                        else 
                            ! Reject all of the above - nothing we can do but set the solution to zero
                            reject_zero = .false.

                            !qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)
                            !qtemp = 0.0_rk*qtemp

                            !call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp)

                            if (allocated(avg_vals)) deallocate(avg_vals)
                            avg_vals = get_neighbor_average_values(worker, worker%mesh, data%sdata%q, ivar, 1)
                            nvals = size(avg_vals,1)
                            nbr_val = 0.0_rk
                            nbr_vol = 0.0_rk
                            do ival = 1, nvals
                                if (ctrl_val>abs(avg_vals(ival,1))) then
                                    nbr_val = nbr_val + avg_vals(ival,1)*avg_vals(ival,2)
                                    nbr_vol = nbr_vol + avg_vals(ival,2)
                                end if
                            end do

                            if (nbr_vol > 0.0_rk) then
                                my_val = nbr_val/nbr_vol
                            else
                                my_val = 0.0_rk
                            end if

                            qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)
                            qtemp = 0.0_rk*qtemp
                            qtemp(1) = my_val

                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp)


                        end if

                    end do
 
                    ! R23 CS
                    unreal_elem = .true.
                    unreal_face = .true.

                    reject_new      = .true.
                    reject_time_avg = .true.
                    reject_nbr_avg  = .true.
                    reject_old      = .true.
                    reject_zero     = .true.
                    
                    do while ( reject_new .and. reject_time_avg .and. reject_nbr_avg .and. reject_old .and. reject_zero)
                        ivar = worker%cache%element%get_field_index('Density * Reynolds-22') 
                        r22 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-33') 
                        r33 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-23') 
                        r23 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')
                        r23_old = interpolate_element_standard(worker%mesh,q0,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                        unreal_elem                 = (any(abs(r23) > sqrt(r22*r33)))
                        unreal_elem_old             = (any(abs(r23_old) > sqrt(r22*r33)))
                        unreal_elem_time_avg        = (any(abs(0.5_rk*(r23+r23_old)) > sqrt(r22*r33)))
                        

                        ctrl_val = minval(sqrt(r22*r33))
                        unreal_face = .false.
                        unreal_face_old = .false.
                        unreal_face_time_avg = .false.
                        do iface = 1,NFACES

                            call worker%set_face(iface)

                            if (allocated(r11f)) deallocate(r11f)
                            if (allocated(r11f_old)) deallocate(r11f_old)

                            ivar = worker%cache%element%get_field_index('Density * Reynolds-23') 
                            r23f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)
                            r23f_old = interpolate_face_standard(worker%mesh,q0,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)

                            ivar = worker%cache%element%get_field_index('Density * Reynolds-22') 
                            r22f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)

                            ivar = worker%cache%element%get_field_index('Density * Reynolds-33') 
                            r33f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)


                            if (any(abs(r23f) > sqrt(r22f*r33f)))                   unreal_face = .true.
                            if (any(abs(r23f_old) > sqrt(r22f*r33f)))               unreal_face_old = .true.
                            if (any(abs(0.5_rk*(r23f+r23f_old)) > sqrt(r22f*r33f))) unreal_face_time_avg = .true.


                            ctrl_val = min(ctrl_val, minval(sqrt(r22f*r33f)))
                        end do  ! faces loop

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-23') 
                        if (.not. (unreal_elem .or. unreal_face)) then
                            ! Accept the new solution without modification
                            reject_new = .false.

                        else if (.not. (unreal_elem_time_avg .or. unreal_face_time_avg)) then
                            ! Accept the time average of the old and new solution
                            reject_time_avg = .false.
                            qtemp0 = q0%dom(idom)%vecs(ielem)%getvar(ivar, 1)

                            qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)

                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, 0.5_rk*(qtemp+qtemp0))

                        !else if (.not. (unreal_elem_old .or. unreal_face_old)) then
                        !    ! Accept the old solution
                        !    reject_old = .false.
                        !    qtemp0 = q0%dom(idom)%vecs(ielem)%getvar(ivar, 1)


                        !    call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp0)


                        else 
                            ! Reject all of the above - nothing we can do but set the solution to zero
                            reject_zero = .false.

                            !qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)
                            !qtemp = 0.0_rk*qtemp

                            !call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp)

                            if (allocated(avg_vals)) deallocate(avg_vals)
                            avg_vals = get_neighbor_average_values(worker, worker%mesh, data%sdata%q, ivar, 1)
                            nvals = size(avg_vals,1)
                            nbr_val = 0.0_rk
                            nbr_vol = 0.0_rk
                            do ival = 1, nvals
                                if (ctrl_val>abs(avg_vals(ival,1))) then
                                    nbr_val = nbr_val + avg_vals(ival,1)*avg_vals(ival,2)
                                    nbr_vol = nbr_vol + avg_vals(ival,2)
                                end if
                            end do

                            if (nbr_vol > 0.0_rk) then
                                my_val = nbr_val/nbr_vol
                            else
                                my_val = 0.0_rk
                            end if

                            qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)
                            qtemp = 0.0_rk*qtemp
                            qtemp(1) = my_val

                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp)


                        end if

                    end do
 
                  
                end if


            end do  ! ielem
            end associate
        end do  ! idom


                    
        ! STEP 3: Check the nonnegativity of the RS determinant
        ! If found to be negative, rescale the off-diagonal RS components.
        ! Iterate until the condition is satisfied.

        do idom = 1,data%mesh%ndomains()
            associate ( domain => data%mesh%domain(idom) )

            ! Loop through elements in the current domain
            do ielem = 1,domain%nelem
                eqn_ID = worker%mesh%domain(idom)%elems(ielem)%eqn_ID


                elem_info%idomain_g  = domain%elems(ielem)%idomain_g
                elem_info%idomain_l  = domain%elems(ielem)%idomain_l
                elem_info%ielement_g = domain%elems(ielem)%ielement_g
                elem_info%ielement_l = domain%elems(ielem)%ielement_l
                call worker%set_element(elem_info)


                ! Update the element cache
                call cache_handler%update(worker,data%eqnset, data%bc_state_group, components    = 'none',           &
                                                                                   face          = NO_ID,           &
                                                                                   differentiate = differentiate,   &
                                                                                   lift          = .false.)
               

                !call eqnset%modify_solution(worker, differentiate)
                eqnset_name = data%get_equation_set_name(eqn_ID)
                if (trim(eqnset_name)=='RANS_RSTM') then
                    unreal_elem_det = .true.
                    unreal_face_det = .true.
                    niter = 0
                    maxiter = 200
                    do while (((unreal_elem_det) .or. (unreal_face_det)) .and. (niter<maxiter) )
                    if (allocated(r11)) deallocate(r11)
                    if (allocated(r22)) deallocate(r22)
                    if (allocated(r33)) deallocate(r33)
                    if (allocated(r12)) deallocate(r12)
                    if (allocated(r13)) deallocate(r13)
                    if (allocated(r23)) deallocate(r23)

                    ivar = worker%cache%element%get_field_index('Density * Reynolds-11') 
                    r11 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                    ivar = worker%cache%element%get_field_index('Density * Reynolds-22') 
                    r22 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                    ivar = worker%cache%element%get_field_index('Density * Reynolds-33') 
                    r33 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                    ivar = worker%cache%element%get_field_index('Density * Reynolds-12') 
                    r12 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                    ivar = worker%cache%element%get_field_index('Density * Reynolds-13') 
                    r13 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')

                    ivar = worker%cache%element%get_field_index('Density * Reynolds-23') 
                    r23 = interpolate_element_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, ivar, worker%itime, 'value')


                    det_elem =  r11*(r22*r33-r23*r23) &
                                -r12*(r12*r33-r23*r13) &
                                +r13*(r12*r23-r22*r13) 


                    unreal_elem_det = .false.
                    unreal_elem_det = (any(det_elem<eps))

                    unreal_face_det = .false.

                    do iface = 1,NFACES

                        call worker%set_face(iface)
                        if (allocated(r11f)) deallocate(r11f)
                        if (allocated(r22f)) deallocate(r22f)
                        if (allocated(r33f)) deallocate(r33f)
                        if (allocated(r12f)) deallocate(r12f)
                        if (allocated(r13f)) deallocate(r13f)
                        if (allocated(r23f)) deallocate(r23f)


                        ivar = worker%cache%element%get_field_index('Density * Reynolds-11') 
                        r11f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-22') 
                        r22f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-33') 
                        r33f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-12') 
                        r12f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-13') 
                        r13f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-23') 
                        r23f = interpolate_face_standard(worker%mesh,worker%solverdata%q,elem_info%idomain_l, elem_info%ielement_l, iface, ivar, worker%itime)


                        det_face =  r11f*(r22f*r33f-r23f*r23f) &
                                -r12f*(r12f*r33f-r23f*r13f) &
                                +r13f*(r12f*r23f-r22f*r13f) 


                        if (any(det_face<eps)) unreal_face_det = .true.

                    end do  ! faces loop
 


                    if ((unreal_elem_det) .or. (unreal_face_det)) then

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-12') 
                        if (allocated(qtemp)) deallocate(qtemp)
                        nterms = data%sdata%q%dom(idom)%vecs(ielem)%nterms()
                        qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)


                        if (niter==(maxiter-1)) then
                            qtemp = 0.0_rk
                        else
                            qtemp = shrink*qtemp 
                        end if


                        call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp)

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-13') 
                        if (allocated(qtemp)) deallocate(qtemp)
                        nterms = data%sdata%q%dom(idom)%vecs(ielem)%nterms()
                        qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)


                        if (niter==(maxiter-1)) then
                            qtemp = 0.0_rk
                        else
                            qtemp = shrink*qtemp 
                        end if


                        call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp)

                        ivar = worker%cache%element%get_field_index('Density * Reynolds-23') 
                        if (allocated(qtemp)) deallocate(qtemp)
                        nterms = data%sdata%q%dom(idom)%vecs(ielem)%nterms()
                        qtemp = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ivar, 1)

                        if (niter==(maxiter-1)) then
                            qtemp = 0.0_rk
                        else
                            qtemp = shrink*qtemp 
                        end if



                        call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar, 1, qtemp)
                    end if


                niter = niter + 1
                end do

                end if


            end do  ! ielem
            end associate
        end do  ! idom



        call loop_timer%stop()




        !
        ! Synchronize
        !
        call MPI_Barrier(ChiDG_COMM,ierr)


        call timer%stop()



        !
        ! Timing IO
        !
        call write_line('Realizability Enforcement Time: ', timer%elapsed(), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<3))
        call write_line('   - RE comm time: ', comm_timer%elapsed(), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<3))
        call write_line('   - RE loop time: ', loop_timer%elapsed(), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<3))
!            call timer%report('Spatial Discretization Time')
!            call comm_timer%report('    - Spatial comm time:')
!            call loop_timer%report('    - Spatial loop time:')

        if (present(timing)) then
            timing = timer%elapsed()
        end if


    end subroutine modify_solution
    !****************************************************************************************










    function get_neighbor_average_values(worker,mesh,vector, ifield, itime) result(avg_vals)
        type(chidg_worker_t),           intent(inout)              :: worker
        type(mesh_t),           intent(in)              :: mesh
        type(chidg_vector_t),   intent(in)              :: vector
        integer(ik),            intent(in)              :: ifield
        integer(ik),            intent(in)              :: itime

        type(face_info_t)   :: iface_info
        type(recv_t)        :: recv_info

        real(rk), allocatable :: avg_vals(:,:)

        real(rk),           allocatable :: qtmp(:)
        real(rk),           allocatable :: interpolator(:,:)
        character(:),       allocatable :: interpolation_style

        integer(ik) :: nderiv, set_deriv, iterm, igq, nterms_s, ierr, nnodes
        logical     :: differentiate_me, conforming_interpolation, chimera_interpolation, parallel_interpolation

        real(rk) :: vol


        ! Chimera data
        integer(ik)                 :: ndonors, idonor
        logical,    allocatable     :: mask(:)          ! node mask for distributing Chimera quadrature points
        type(AD_D), allocatable     :: var_gq_chimera(:)
        integer(ik) :: iface, nvals, ival, myfacetype
        
        !nnodes   = mesh%domain(face_info%idomain_l)%elems(face_info%ielement_l)%basis_s%nnodes_if()
        !nterms_s = mesh%domain(face_info%idomain_l)%elems(face_info%ielement_l)%basis_s%nterms_i()

       
        !print *, 'getting nbr vals'
        ! Loop over the faces of the current element to see how many neighbor elements it has.
        nvals = 0
        do iface = 1, NFACES
            call worker%set_face(iface)
            ! Get number of donors for the interpolation
            myfacetype =  worker%mesh%domain(worker%element_info%idomain_l)%faces(worker%element_info%ielement_l, worker%iface)%ftype
            ! Only INTERIOR or CHIMERA faces are associated with neighbor elements.
            if ((myfacetype==INTERIOR) .or. (myfacetype==CHIMERA)) then
            ndonors = get_face_interpolation_ndonors(mesh,worker%face_info(),NEIGHBOR)
            nvals = nvals+ndonors
            end if
        end do


        ! Now, we can allocate the storage
        ! avg_val(ival,1) = average value, avg_val(ival,2) = donor volume
        allocate(avg_vals(nvals,2), stat=ierr)
        if (ierr /= 0) call AllocationError


        ival = 0
        do iface = 1, NFACES
            call worker%set_face(iface)
            myfacetype =  worker%mesh%domain(worker%element_info%idomain_l)%faces(worker%element_info%ielement_l, worker%iface)%ftype
            if ((myfacetype==INTERIOR) .or. (myfacetype==CHIMERA)) then
            ndonors = get_face_interpolation_ndonors(mesh,worker%face_info(),NEIGHBOR)
            !
            ! Get interpolation style. Conforming or Chimera
            !
            interpolation_style = get_face_interpolation_style(mesh,worker%face_info(),NEIGHBOR)
            conforming_interpolation = (interpolation_style == 'conforming')
            chimera_interpolation    = (interpolation_style == 'chimera')




            !
            ! For each donor element to the interpolation. 
            ! (ndonors could be > 1 for Chimera interpolations)
            !
            do idonor = 1,ndonors
                ival = ival + 1

                !
                ! Get face info for face being interpolated to(ME, NEIGHBOR), 
                ! interpolation matrix, and recv data for parallel access
                !
                iface_info   = get_face_interpolation_info(        mesh,worker%face_info(),NEIGHBOR,idonor)
                mask         = get_face_interpolation_mask(        mesh,worker%face_info(),NEIGHBOR,idonor)
                recv_info    = get_face_interpolation_comm(        mesh,worker%face_info(),NEIGHBOR,idonor)
                !interpolator = get_face_interpolation_interpolator(mesh,face_info,NEIGHBOR,idonor,'value',iface_info)

                parallel_interpolation = (recv_info%comm /= 0)

            
                !
                ! Retrieve modal coefficients for ifield from vector
                !
                if (parallel_interpolation) then
                    qtmp = vector%recv%comm(recv_info%comm)%dom(recv_info%domain)%vecs(recv_info%element)%getvar(ifield,itime)
                    ! Apparently, element volume is not yet part of the parallel communication
                    ! We'd want something like this:
                    !vol = vector%recv%comm(recv_info%comm)%dom(recv_info%domain)%elems(recv_info%element)%vol_ale

                    ! To get around this for now, just use the volume of the element itself
                    vol = mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%vol_ale

                else
                    qtmp = vector%dom(iface_info%idomain_l)%vecs(iface_info%ielement_l)%getvar(ifield,itime)
                    vol = mesh%domain(iface_info%idomain_l)%elems(iface_info%ielement_l)%vol_ale
                end if

                avg_vals(ival, 1) = qtmp(1)
                avg_vals(ival, 2) = vol 
                

            end do ! idonor
            end if
        end do ! iface



    end function get_neighbor_average_values
    !*****************************************************************************************






end module mod_modify_solution
