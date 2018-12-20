module mod_update_grid
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NO_MM_ASSIGNED
    use mod_chidg_mpi,          only: GLOBAL_MASTER
    use type_chidg_data,        only: chidg_data_t
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
        type(chidg_data_t), intent(inout)   :: data
        real(rk),           optional        :: timing
        integer(ik),        optional        :: info

        integer(ik) :: mm_ID, idom, ierr, imm


        call write_line('Updating grid according to mesh motion...', io_proc=GLOBAL_MASTER)


        ! Update mesh motion objects
        do imm = 1, size(data%mesh_motion) 
            call data%mesh_motion(imm)%mm%update(data)
            call data%mesh_motion(imm)%mm%apply(data)
        end do

        ! Loop through domains and apply mesh motions
        do idom = 1,data%mesh%ndomains()
            associate ( mesh => data%mesh)
            mm_ID = mesh%domain(idom)%mm_ID

            if (mm_ID /= NO_MM_ASSIGNED) then
                call mesh%domain(idom)%set_displacements_velocities(mesh%domain(idom)%dnodes, mesh%domain(idom)%vnodes)
                call mesh%domain(idom)%update_interpolations_ale()
            end if
            end associate
        end do  ! idom





! OLD
!        ! Loop through domains
!        do idom = 1,data%mesh%ndomains()
!            associate ( mesh => data%mesh)
!            pmm_ID = mesh%domain(idom)%pmm_ID
!
!            if (pmm_ID /= NO_MM_ASSIGNED) then
!
!                do inode = 1, size(mesh%domain(idom)%nodes,1)
!                    mesh%domain(idom)%dnodes(inode,:) = data%pmm(pmm_ID)%pmmf%compute_pos(data%time_manager%t, mesh%domain(idom)%nodes(inode,:)) - mesh%domain(idom)%nodes(inode,:)
!                    mesh%domain(idom)%vnodes(inode,:) = data%pmm(pmm_ID)%pmmf%compute_vel(data%time_manager%t, mesh%domain(idom)%nodes(inode,:))
!                end do
!
!                call mesh%domain(idom)%set_displacements_velocities(mesh%domain(idom)%dnodes, mesh%domain(idom)%vnodes)
!                call mesh%domain(idom)%update_interpolations_ale()
!            end if
!            end associate
!        end do  ! idom


        call data%mesh%comm_send()
        call data%mesh%comm_recv()
        call data%mesh%comm_wait()


    end subroutine update_grid
    !******************************************************************************************************************















end module mod_update_grid
