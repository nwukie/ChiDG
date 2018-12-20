module type_rbf_node_patch
#include <messenger.h>
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: ZERO, ONE, TWO, NO_ID 
    use type_mesh,          only: mesh_t
    use type_svector,       only: svector_t
    use mod_string
    use mpi_f08
    use mod_chidg_mpi
    implicit none

    type :: rbf_node_patch_t

        integer(ik)                 :: nnodes_patch = 0

        real(rk), allocatable       :: nodes(:,:) !(nnodes_patch, 3)
        integer(ik), allocatable    :: rbf_to_patch_index(:) !inode = rbf_to_patch(inode_rbf) is the patch index
        integer(ik), allocatable    :: patch_to_rbf_index(:) !inode_rbf = patch_to_rbf_index(inode) is the RBF system index
        integer(ik), allocatable    :: idomain_g(:)
        integer(ik), allocatable    :: idomain_l(:)
        integer(ik), allocatable    :: ielement_g(:)
        integer(ik), allocatable    :: ielement_l(:)
        integer(ik), allocatable    :: iface(:)
        integer(ik), allocatable    :: face_node_index(:) ! The index of the node in its face%interp_coords(:,:) 
        integer(ik), allocatable    :: patch_ID(:) ! Index of the origin patch of each node

    contains

        procedure :: init
        procedure :: assemble_rbf_patch
        procedure :: assemble_from_array
        procedure :: reorder_patch_to_rbf
        procedure :: reorder_rbf_to_patch

    end type rbf_node_patch_t

contains

    subroutine init(self,nnodes_patch_in)
        class(rbf_node_patch_t), intent(inout) :: self
        integer(ik),            intent(in)      :: nnodes_patch_in

        self%nnodes_patch = nnodes_patch_in
        allocate(self%nodes(nnodes_patch_in,3), self%rbf_to_patch_index(nnodes_patch_in), &
            self%patch_to_rbf_index(nnodes_patch_in), &
            self%idomain_g(nnodes_patch_in), self%idomain_l(nnodes_patch_in),&
            self%ielement_g(nnodes_patch_in), self%ielement_l(nnodes_patch_in),&
            self%iface(nnodes_patch_in), self%face_node_index(nnodes_patch_in), &
            self%patch_ID(nnodes_patch_in))

    end subroutine init

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
    subroutine assemble_rbf_patch(self, mesh, patch_groups)
        class(rbf_node_patch_t), intent(inout)               ::  self
        type(mesh_t), intent(inout)                         :: mesh 
        type(svector_t),    intent(inout)                  :: patch_groups

        character(:), allocatable                       :: patch_group
    
        integer(ik)                 :: group_ID, patch_ID, face_ID, &
                                       idomain_g,  idomain_l,        &
                                       ielement_g, ielement_l, iface, ierr, &
                                       ipatch, npatches, inode, nnodes_patch, node_count

        real(rk),   allocatable :: rbf_nodes(:,:)
        integer(ik), allocatable :: rbf_idomain_g(:), rbf_idomain_l(:), rbf_ielement_g(:), rbf_ielement_l(:), &
            rbf_iface(:), rbf_node_index(:), rbf_patch_ID(:)

        type(string_t) :: tmpstr

        call write_line('Assembling RBF source patch...', io_proc=GLOBAL_MASTER)

        npatches = patch_groups%size()
        
        do ipatch = 1, npatches

            if (allocated(patch_group)) deallocate(patch_group)
            tmpstr = patch_groups%at(ipatch)
            patch_group = tmpstr%get()
            !
            ! Get patch_group boundary group ID
            !
            group_ID = mesh%get_bc_patch_group_id(trim(patch_group))


            !
            ! Loop over domains/elements/faces for "patch_group" 
            !

            if (group_ID /= NO_ID) then
                do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
                    do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                        idomain_g  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_g()
                        idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                        ielement_g = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_g(face_ID)
                        ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                        iface      = mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)


                        ! Collect the number of nodes
                        nnodes_patch = nnodes_patch + &
                            size(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords,1)

                    end do !iface
                end do !ipatch
            end if ! group_ID /= NO_ID
        end do

        allocate(rbf_nodes(nnodes_patch, 3), rbf_idomain_g(nnodes_patch), rbf_idomain_l(nnodes_patch), &
            rbf_ielement_g(nnodes_patch), rbf_ielement_l(nnodes_patch), rbf_iface(nnodes_patch), &
            rbf_node_index(nnodes_patch), rbf_patch_ID(nnodes_patch))

        node_count = 0
        do ipatch = 1, npatches

            if (allocated(patch_group)) deallocate(patch_group)
            tmpstr = patch_groups%at(ipatch)
            patch_group = tmpstr%get()
            !
            ! Get patch_group boundary group ID
            !
            group_ID = mesh%get_bc_patch_group_id(trim(patch_group))


            if (group_ID /= NO_ID) then
                do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
                    do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                        idomain_g  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_g()
                        idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                        ielement_g = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_g(face_ID)
                        ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                        iface      = mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)

                        do inode = 1, size(mesh%domain(idomain_l)%faces(ielement_l, iface)%interp_coords,1)
                            node_count = node_count + 1

                            rbf_nodes(node_count,:)     = mesh%domain(idomain_l)%faces(ielement_l, iface)%interp_coords(inode,:)
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

        call self%init(nnodes_patch)
        self%nodes        = rbf_nodes
        self%idomain_g    = rbf_idomain_g
        self%idomain_l    = rbf_idomain_l
        self%ielement_g   = rbf_ielement_g
        self%ielement_l   = rbf_ielement_l
        self%iface        = rbf_iface
        self%patch_ID     = rbf_patch_ID




    end subroutine assemble_rbf_patch 
    !******************************************************************************************


!    !> Assemble an RBF node patch from a BC patch_group 
!    !!
!    !!
!    !!
!    !!
!    !!  @author Eric Wolf 
!    !!  @date   10/19/2017
!    !!  @note   Modified directly from 'report_aerodynamics' 
!    !!
!    !!
!    !!  @param[in]      data            chidg_data instance
!    !!  @param[in]      patch_group     Name of patch group over which the force will be integrated.
!    !!
!    !-----------------------------------------------------------------------------------
!    subroutine assemble(self,data,patch_group)
!        class(rbf_node_patch_t), intent(inout) :: self
!        type(chidg_data_t), intent(inout)               :: data
!        character(*),       intent(in)                  :: patch_group
!    
!        integer(ik)                 :: group_ID, patch_ID, face_ID, &
!                                       idomain_g,  idomain_l,        &
!                                       ielement_g, ielement_l, iface, ierr
!
!        integer(ik) :: inode, nnodes_patch, node_count
!
!
!        real(rk), allocatable :: rbf_nodes(:,:)
!
!        integer(rk), allocatable  :: rbf_idomain_g(:), rbf_idomain_l(:), rbf_ielement_g(:), rbf_ielement_l(:), &
!            rbf_iface(:), rbf_face_node_index(:)
!
!        !call write_line('Assembling RBF source patch...', io_proc=GLOBAL_MASTER)
!
!        !
!        ! Get patch_group boundary group ID
!        !
!        group_ID = data%mesh%get_bc_patch_group_id(trim(patch_group))
!
!
!        !
!        ! Loop over domains/elements/faces for "patch_group" 
!        !
!
!        if (group_ID /= NO_ID) then
!            do patch_ID = 1,data%mesh%bc_patch_group(group_ID)%npatches()
!                do face_ID = 1,data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()
!
!                    idomain_g  = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_g()
!                    idomain_l  = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
!                    ielement_g = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_g(face_ID)
!                    ielement_l = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
!                    iface      = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)
!
!
!                    ! Collect the number of nodes
!                    nnodes_patch = nnodes_patch + &
!                        size(data%mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords,1)
!
!                end do !iface
!            end do !ipatch
!
!            allocate(rbf_nodes(nnodes_patch, 3), rbf_idomain_g(nnodes_patch), rbf_idomain_l(nnodes_patch), &
!                rbf_ielement_g(nnodes_patch), rbf_ielement_l(nnodes_patch), rbf_iface(nnodes_patch), rbf_face_node_index(nnodes_patch))
!
!            node_count = 0
!            do patch_ID = 1,data%mesh%bc_patch_group(group_ID)%npatches()
!                do face_ID = 1,data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()
!
!                    idomain_g  = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_g()
!                    idomain_l  = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
!                    ielement_g = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_g(face_ID)
!                    ielement_l = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
!                    iface      = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)
!
!                    do inode = 1, size(data%mesh%domain(idomain_l)%faces(ielement_l, iface)%interp_coords,1)
!                        node_count = node_count + 1
!
!                        rbf_nodes(node_count,:)     = data%mesh%domain(idomain_l)%faces(ielement_l, iface)%interp_coords(inode,:)
!                        rbf_idomain_g(node_count)   = idomain_g
!                        rbf_idomain_l(node_count)   = idomain_l
!                        rbf_ielement_g(node_count)  = ielement_g
!                        rbf_ielement_l(node_count)  = ielement_l
!                        rbf_iface(node_count)       = iface
!                        rbf_face_node_index(node_count)  = inode
!
!                    end do
!
!                end do !iface
!            end do !ipatch
!
!            !
!            ! Now, store into the rbf_node_patch data structure
!            !
!
!            call self%init(nnodes_patch)
!            self%nodes      = rbf_nodes
!            self%idomain_g  = rbf_idomain_g
!            self%idomain_l  = rbf_idomain_l
!            self%ielement_g = rbf_ielement_g
!            self%ielement_l = rbf_ielement_l
!            self%iface      = rbf_iface
!            self%face_node_index = rbf_face_node_index
!
!        end if ! group_ID /= NO_ID
!
!
!
!    end subroutine assemble
!    !******************************************************************************************


    !> Assemble an RBF node patch from a BC patch_group 
    !!
    !!
    !!
    !!
    !!  @author Eric Wolf 
    !!  @date   10/19/2017
    !!  @note   Modified directly from 'report_aerodynamics' 
    !!
    !!
    !!  @param[in]      data            chidg_data instance
    !!  @param[in]      patch_group     Name of patch group over which the force will be integrated.
    !!
    !-----------------------------------------------------------------------------------
    subroutine assemble_from_array(self,source_nodes_in)
        class(rbf_node_patch_t), intent(inout) :: self
        real(rk),                   intent(in)  :: source_nodes_in(:,:)
    
        integer(ik)                 :: nnodes_patch



        !call write_line('Assembling RBF source patch...', io_proc=GLOBAL_MASTER)

        nnodes_patch = size(source_nodes_in,1)

        call self%init(nnodes_patch)
        self%nodes      = source_nodes_in




    end subroutine assemble_from_array
    !******************************************************************************************

    !>  Reorders a rank-1 array from patch node ordering to RBF node ordering.
    !!
    !!  @author Eric Wolf 
    !!  @date   10/19/2017
    !!
    !-----------------------------------------------------------------------------------
    function reorder_patch_to_rbf(self,array_in) result(array_out)
        class(rbf_node_patch_t), intent(inout) :: self
        real(rk),                   intent(in)  :: array_in(:)

        real(rk), allocatable :: array_out(:)


        integer(ik)                 :: inode, inode_rbf, nnodes_patch

        nnodes_patch = size(array_in,1)

        allocate(array_out(nnodes_patch))

        do inode = 1, nnodes_patch

            inode_rbf = self%patch_to_rbf_index(inode)
            array_out(inode_rbf) = array_in(inode)

        end do

    end function reorder_patch_to_rbf
    !******************************************************************************************
    
    !>  Reorders a rank-1 array from RBF node ordering to patch node ordering.
    !!
    !!  @author Eric Wolf 
    !!  @date   10/19/2017
    !!
    !-----------------------------------------------------------------------------------
    function reorder_rbf_to_patch(self,array_in) result(array_out)
        class(rbf_node_patch_t), intent(inout) :: self
        real(rk),                   intent(in)  :: array_in(:)

        real(rk), allocatable :: array_out(:)


        integer(ik)                 :: inode, inode_rbf, nnodes_patch

        nnodes_patch = size(array_in,1)

        allocate(array_out(nnodes_patch))

        do inode_rbf = 1, nnodes_patch

            inode = self%rbf_to_patch_index(inode_rbf)
            array_out(inode) = array_in(inode_rbf)

        end do

    end function reorder_rbf_to_patch
    !******************************************************************************************
 
end module type_rbf_node_patch
