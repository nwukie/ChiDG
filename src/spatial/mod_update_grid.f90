module mod_update_grid
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES, DIAG, CHIMERA, INTERIOR, NO_PMM_ASSIGNED
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM, GLOBAL_MASTER
    use mod_io,                 only: verbosity
    use mpi_f08,                only: MPI_Barrier



    use type_chidg_data,        only: chidg_data_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_chidg_cache,       only: chidg_cache_t
    use type_cache_handler,     only: cache_handler_t
    use type_element_info,      only: element_info_t
    use type_timer,             only: timer_t
    implicit none


contains





    !>  Spatial loop through domains, elements, and faces. Functions get called for each element/face.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/15/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!  @note   Improved layout, added computation of diffusion terms.
    !!
    !!
    !------------------------------------------------------------------------------------------------------------------
    subroutine update_grid(data,timing,info)
        type(chidg_data_t), intent(inout), target :: data
        real(rk),           optional        :: timing
        integer(ik),        optional        :: info

        type(timer_t)               :: timer, comm_timer
        integer(ik)                 :: inode, pmm_ID, idom, ielem, iface, idiff, ifcn, ibc, ierr, nelem
        logical                     :: interior_face
        logical                     :: chimera_face 
        logical                     :: compute_face 



        type(chidg_worker_t)        :: worker
        type(element_info_t)        :: elem_info

        type(chidg_cache_t)         :: cache
        type(cache_handler_t)       :: cache_handler


        ! Initialize Chidg Worker references
        call worker%init(data%mesh, data%eqnset(:)%prop, data%sdata, cache)


        !
        ! Start timer on spatial discretization update
        !
        call timer%start()

    

        !------------------------------------------------------------------------------------------
        !                                      Interior Scheme
        !------------------------------------------------------------------------------------------

        !
        ! Clear function_status data. This tracks if a function has already been called. So, in this way
        ! we can compute a function on a face and apply it to both elements. The function is just registered
        ! as computed for both. So we need to reset all of that data here. This is only tracked for the interior scheme.
        ! Boundary condition evaluations and Chimera faces are not tracked.
        !
        !call data%sdata%function_status%clear()




        !
        ! Communicate solution vector
        !
        !call comm_timer%start()
        !call data%sdata%q%comm_send()
        !call data%sdata%q%comm_recv()
        !call data%sdata%q%comm_wait()
        !call comm_timer%stop()




        !
        ! Loop through given element compute the residual functions and also the linearization of those functions
        !
    
        call write_line('Updating grid according to mesh motion...', io_proc=GLOBAL_MASTER)


        worker%itime = data%time_manager%itime
        worker%t     = data%time_manager%t
        ! Loop through domains
        do idom = 1,data%mesh%ndomains()
            associate ( mesh => data%mesh)
            pmm_ID = mesh%domain(idom)%pmm_ID

            if (pmm_ID /= NO_PMM_ASSIGNED) then

                        do inode = 1, size(mesh%domain(idom)%nodes,1)
                    mesh%domain(idom)%dnodes(inode,:) = &
                        data%pmm(pmm_ID)%pmmf%compute_pos(data%time_manager%t, mesh%domain(idom)%nodes(inode,:)) - &
                        mesh%domain(idom)%nodes(inode,:)
                    mesh%domain(idom)%vnodes(inode,:) = data%pmm(pmm_ID)%pmmf%compute_vel(data%time_manager%t, mesh%domain(idom)%nodes(inode,:))

                end do

                call mesh%domain(idom)%init_ale(mesh%domain(idom)%dnodes, mesh%domain(idom)%vnodes)
                call mesh%domain(idom)%update_ale()
            end if
            end associate
        end do  ! idom





        !
        ! Synchronize
        !
        call MPI_Barrier(ChiDG_COMM,ierr)





!        call timer%stop()
!        call timer%report('Grid Update Time')
!        !call comm_timer%report('    - Grid update comm time:')
!        if (present(timing)) then
!            timing = timer%elapsed()
!        end if


    end subroutine update_grid
    !******************************************************************************************************************















end module mod_update_grid
