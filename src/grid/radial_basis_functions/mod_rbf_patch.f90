!>  Compute
!!
!!  @author Nathan A. Wukie (AFRL)
!!  @date   7/12/2017
!!  @note   Modified directly from 'chidg airfoil' action
!!
!!
!---------------------------------------------------------------------------------------------
module mod_rbf_patch
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO, TWO, NO_ID
    use mod_chidg_mpi,          only: ChiDG_COMM
    use type_chidg_data,        only: chidg_data_t
    use type_element_info,      only: element_info_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_chidg_cache,       only: chidg_cache_t
    use type_cache_handler,     only: cache_handler_t
    use mpi_f08,                only: MPI_AllReduce, MPI_REAL8, MPI_SUM
    use DNAD_D
    implicit none







contains



    !>  Compute force integrated over a specified patch group.
    !!
    !!
    !!  F = int[ (tau-p) dot n ] dPatch
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/12/2017
    !!  @note   Modified directly from 'chidg airfoil' action
    !!
    !!
    !!  @param[in]      data            chidg_data instance
    !!  @param[in]      patch_group     Name of patch group over which the force will be integrated.
    !!  @result[out]    force           Integrated force vector: force = [f1, f2, f3]
    !!
    !-----------------------------------------------------------------------------------
    subroutine assemble_rbf_patch(data,patch_groups)
        type(chidg_data_t), intent(inout)               :: data
        type(svector_t),    intent(in)                  :: patch_groups

        character(*), allocatable                       :: patch_group
    
        integer(ik)                 :: group_ID, patch_ID, face_ID, &
                                       idomain_g,  idomain_l,        &
                                       ielement_g, ielement_l, iface, ierr



        call write_line('Assembling RBF source patch...', io_proc=GLOBAL_MASTER)

        npatches = patch_groups%size()
        
        do ipatch = 1, npatches

            patch_group = trim(patch_groups%at(ipatch))
            !
            ! Get patch_group boundary group ID
            !
            group_ID = data%mesh%get_bc_patch_group_id(trim(patch_group))


            !
            ! Loop over domains/elements/faces for "patch_group" 
            !

            if (group_ID /= NO_ID) then
                do patch_ID = 1,data%mesh%bc_patch_group(group_ID)%npatches()
                    do face_ID = 1,data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                        idomain_g  = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_g()
                        idomain_l  = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                        ielement_g = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_g(face_ID)
                        ielement_l = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                        iface      = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)


                        ! Collect the number of nodes
                        nnodes_patch = nnodes_patch + &
                            size(data%mesh%domain(idomain_l)%element(ielement_l)%faces(iface)%interp_coords,1)

                    end do !iface
                end do !ipatch
            end if ! group_ID /= NO_ID
        end do

        allocate(rbf_nodes(nnodes_patch, 3), rbf_idomain_g(nnodes_patch), rbf_idomain_l(nnodes_patch), &
            rbf_ielement_g(nnodes_patch), rbf_ielement_l(nnodes_patch), rbf_iface(nnodes_patch), &
            rbf_node_index(nnodes_patch), rbf_patch_ID(nnodes_patch))

        node_count = 0
        do ipatch = 1, npatches

            patch_group = trim(patch_groups%at(ipatch))
            !
            ! Get patch_group boundary group ID
            !
            group_ID = data%mesh%get_bc_patch_group_id(trim(patch_group))


            if (group_ID /= NO_ID) then
                do patch_ID = 1,data%mesh%bc_patch_group(group_ID)%npatches()
                    do face_ID = 1,data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                        idomain_g  = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_g()
                        idomain_l  = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                        ielement_g = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_g(face_ID)
                        ielement_l = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                        iface      = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)

                        do inode = 1, size(data%mesh%domain(idomain_l)%faces(ielement_l, iface)%interp_coords,1)
                            node_count = node_count + 1

                            rbf_nodes(node_count,:)     = data%mesh%domain(idomain_l)%faces(ielement_l, iface)%interp_coords(inode,:)
                            rbf_idomain_g(node_count)   = idomain_g
                            rbf_idomain_l(node_count)   = idomain_l
                            rbf_ielement_g(node_count)  = ielement_g
                            rbf_ielement_l(node_count)  = ielement_l
                            rbf_iface(node_count)       = iface
                            rbf_node_index(node_count)  = inode
                            rbf_patch_ID(node_count)    = ipatch

                        end do

                    end do !iface
                end do !ipatch

            end if ! group_ID /= NO_ID
        end do
        !
        ! Now, store into the rbf_node_patch data structure
        !

        call rbf_node_patch%init(nnodes_patch)
        rbf_node_patch%nodes        = rbf_nodes
        rbf_node_patch%idomain_g    = rbf_idomain_g
        rbf_node_patch%idomain_l    = rbf_idomain_l
        rbf_node_patch%ielement_g   = rbf_ielement_g
        rbf_node_patch%ielement_l   = rbf_ielement_l
        rbf_node_patch%iface        = rbf_iface
        rbf_node_patch%patch_ID     = rbf_patch_ID




    end subroutine assemble_rbf_patch 
    !******************************************************************************************


    



end module mod_rbf_patch
