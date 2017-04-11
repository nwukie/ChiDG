module mod_spatial
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES, DIAG, CHIMERA, INTERIOR
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM, GLOBAL_MASTER
    use mpi_f08,                only: MPI_Barrier



    use type_chidg_data,        only: chidg_data_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_chidg_cache,       only: chidg_cache_t
    use type_cache_handler,     only: cache_handler_t
    use type_element_info,      only: element_info_t
    use type_timer,             only: timer_t
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
    subroutine update_space(data,timing,differentiate)
        type(chidg_data_t), intent(inout),  target      :: data
        real(rk),           intent(inout),  optional    :: timing
        logical,            intent(in),     optional    :: differentiate

        integer(ik)                 :: idom, ielem, iface, idiff, itime, ierr, &
                                       diff_min, diff_max, eqn_ID
        type(timer_t)               :: timer, comm_timer, loop_timer

        type(element_info_t)        :: elem_info

        type(chidg_worker_t)        :: worker
        type(chidg_cache_t)         :: cache
        type(cache_handler_t)       :: cache_handler

        logical                     :: differentiate_function, io_spatial
        character(:),   allocatable :: time_string

        
        !
        ! Decide whether to differentiate the discretization or not
        !
        if (present(differentiate)) then
            ! User-specified
            differentiate_function = differentiate
        else
            ! Default, differentiate
            differentiate_function = .true.
        end if


        !
        ! Initialize Chidg Worker references
        !
        call worker%init(data%mesh, data%eqnset(:)%prop, data%sdata, cache)


        !
        ! Start timer on spatial discretization update
        !
        call timer%start()

    

        !--------------------------------------------------------------------------------------
        !                                    Interior Scheme
        !--------------------------------------------------------------------------------------

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
        call write_line('-  Updating spatial scheme', io_proc=GLOBAL_MASTER)


        ! Loop through domains
        do itime = 1,data%ntime()

            !
            ! Clear function_status data. This tracks if a function has already been called. So, in this way
            ! we can compute a function on a face and apply it to both elements. The function is just registered
            ! as computed for both. So we need to reset all of that data here. This is only tracked for the interior scheme.
            ! Boundary condition evaluations and Chimera faces are not tracked.
            !
            call data%sdata%function_status%clear()


            call loop_timer%start()
            do idom = 1,data%mesh%ndomains()
                eqn_ID = worker%mesh%domain(idom)%eqn_ID
                associate ( domain => data%mesh%domain(idom), eqnset => data%eqnset(eqn_ID) )

                ! Loop through elements in the current domain
                do ielem = 1,domain%nelem


                    elem_info%idomain_g  = domain%elems(ielem)%idomain_g
                    elem_info%idomain_l  = domain%elems(ielem)%idomain_l
                    elem_info%ielement_g = domain%elems(ielem)%ielement_g
                    elem_info%ielement_l = domain%elems(ielem)%ielement_l
                    call worker%set_element(elem_info)

                    worker%itime = itime
                    time_string = data%time_manager%get_name()

                    worker%t = data%sdata%t

                    ! Update the element cache
                    call cache_handler%update(worker,data%eqnset, data%bc_state_group, differentiate_function)



                    ! Faces loop. For the current element, compute the 
                    ! contributions from boundary integrals.
                    do iface = 1,NFACES

                        call worker%set_face(iface)
 
                        call eqnset%compute_boundary_advective_operators(worker, differentiate_function)
                        call eqnset%compute_boundary_diffusive_operators(worker, differentiate_function)
                        call eqnset%compute_bc_operators(worker,data%bc_state_group, differentiate_function)

                    end do  ! faces loop
                    


                    !
                    ! Compute contributions from volume integrals
                    !
                    call eqnset%compute_volume_advective_operators(worker, differentiate_function)
                    call eqnset%compute_volume_diffusive_operators(worker, differentiate_function)


                end do  ! ielem
                end associate
            end do  ! idom
        end do ! itime
        call loop_timer%stop()




        !
        ! Synchronize
        !
        call MPI_Barrier(ChiDG_COMM,ierr)





        call timer%stop()

        io_spatial = .true.
        if (io_spatial) then
            call timer%report('Spatial Discretization Time')
            call comm_timer%report('    - Spatial comm time:')
            call loop_timer%report('    - Spatial loop time:')
        end if

        if (present(timing)) then
            timing = timer%elapsed()
        end if


    end subroutine update_space
    !****************************************************************************************















end module mod_spatial
