module type_bc_patch_group
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: NO_ID, BOUNDARY, ORPHAN, NFACES, NO_FACE
    use mod_grid,                   only: FACE_CORNERS
    use type_point,                 only: point_t
    use type_bc_patch,              only: bc_patch_t
    use type_domain,                only: domain_t
    use type_boundary_connectivity, only: boundary_connectivity_t
    use type_ivector,               only: ivector_t
    implicit none




    !>  An object containing all the bc_patch_t instances for a single
    !!  boundary condition.
    !!
    !!  The bc_patch_group can contain patches from multiple domains. The individual 
    !!  patches are on self%patch(:).
    !!
    !!  The self%group_ID component is used externally to access the correct group
    !!  as mesh%bc_patch_group(group_ID). When a new bc_patch is initialized on the
    !!  group, the bc_patch_group sets the group_ID on the appropriate domain faces
    !!  so they know where to access the patch group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !-----------------------------------------------------------------------------------------
    type, public :: bc_patch_group_t

        character(:),       allocatable :: name

        type(bc_patch_t),   allocatable :: patch(:)

        integer(ik)                     :: group_ID

    contains

        ! Patch routines
        procedure           :: add_bc_patch
        procedure           :: npatches
        procedure           :: nfaces => get_nfaces
        procedure,  private :: new_bc_patch

        ! Parallel communication patterns
        procedure           :: get_recv_procs
        procedure           :: get_send_procs

    end type bc_patch_group_t
    !*****************************************************************************************







contains


    !>  Initialize boundary condition patch.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/19/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_bc_patch(self,domain,patch_name,bc_connectivity,bc_ID)
        class(bc_patch_group_t),        intent(inout)           :: self
        type(domain_t),                 intent(inout)           :: domain
        character(*),                   intent(in)              :: patch_name
        type(boundary_connectivity_t),  intent(in)              :: bc_connectivity
        integer(ik),                    intent(in)              :: bc_ID


        type(point_t)               :: pnt, point_one, point_two, point_three
        character(:),   allocatable :: bcname, user_msg
        real(rk)                    :: time, x, y, z
        integer(ik)                 :: nelem_xi, nelem_eta, nelem_zeta, nelem_bc, ielem_bc,         & 
                                       xi_begin, eta_begin, zeta_begin, xi_end, eta_end, zeta_end,  & 
                                       ixi, ieta, izeta, ierr, ielem, ielem_test, nface_nodes,      &
                                       iface, inode, i, nfaces_bc, iface_bc, face_ID, patch_ID,     &
                                       ncoupled_elements, try_face, nterms_1d, mapping

        logical,        allocatable :: node_matched(:), xi_face, eta_face, zeta_face
        integer(ik),    allocatable :: element_nodes(:), face_nodes(:)
        integer(ik)                 :: face_node
        logical                     :: includes_corner_one, includes_corner_two, &
                                       includes_corner_three, includes_corner_four, face_matches_boundary
        integer(ik)                 :: corner_one, corner_two, corner_three, corner_four
        integer(ik)                 :: corner_indices(4)




        !
        ! Get number of elements/faces associated with boundary condition.
        !
        nelem_bc = bc_connectivity%nfaces()


        !
        ! Loop through each face in bc connectivity:
        !
        !   - Find owner element, determine iface
        !
        iface_bc = 0
        patch_ID = NO_ID
        do ielem_bc = 1,nelem_bc

            !
            ! Get nodes from face
            !
            face_nodes = bc_connectivity%data(ielem_bc)%get_face_nodes()


            nface_nodes = size(face_nodes)
            if ( allocated(node_matched) ) deallocate(node_matched)
            allocate(node_matched(nface_nodes), stat=ierr)
            if (ierr /= 0) call AllocationError

            
            !
            ! Search for local element with common nodes
            !
            do ielem = 1,domain%nelem

                
                !
                ! Get nodes from element being tested
                !
                element_nodes = domain%elems(ielem)%connectivity


                !
                ! Loop through bc nodes and see if current element has them all
                !
                node_matched = .false.
                do inode = 1,nface_nodes
                    if ( any(element_nodes == face_nodes(inode) ) ) then
                        node_matched(inode) = .true.
                    end if
                end do


                !
                ! Determine element mapping index
                !
                nterms_1d = 0
                do while (nterms_1d*nterms_1d*nterms_1d < domain%elems(ielem)%nterms_c)
                    nterms_1d = nterms_1d + 1
                end do
                mapping = nterms_1d - 1



        
                !
                ! If all match: 
                !   - determine face index, iface, corresponding to the element boundary
                !   - set element/face indices in boundary condition list.
                !
                if ( all(node_matched) ) then
                    iface_bc = iface_bc + 1

                    !
                    ! Since a geometry region was detected on the current processor,
                    ! create a new bc_patch if one hasn't already been created.
                    ! We want this creation statement in the inner loop here so that
                    ! only those processor that contain faces on the boundary 
                    ! actually allocate a patch.
                    !
                    ! Also, store global connectivity list for IO purposes.
                    !
                    if (patch_ID == NO_ID) then
                        patch_ID = self%new_bc_patch()
                        call self%patch(patch_ID)%init(patch_name,patch_ID,domain%idomain_g,domain%idomain_l,bc_connectivity)
                    end if


                    !
                    ! Detect element face index associated with the boundary condition
                    !   Goal: determine iface
                    !
                    !   Loop through element faces until a match is found
                    !
                    iface = NO_FACE
                    do try_face = 1,NFACES 

                        !
                        ! Get corner nodes for face, try_face
                        !
                        corner_one   = FACE_CORNERS(try_face,1,mapping)
                        corner_two   = FACE_CORNERS(try_face,2,mapping)
                        corner_three = FACE_CORNERS(try_face,3,mapping)
                        corner_four  = FACE_CORNERS(try_face,4,mapping)


                        !
                        ! For the element, get the global indices of the corner nodes for face, try_face
                        !
                        corner_indices(1) = domain%elems(ielem)%connectivity(corner_one)
                        corner_indices(2) = domain%elems(ielem)%connectivity(corner_two)
                        corner_indices(3) = domain%elems(ielem)%connectivity(corner_three)
                        corner_indices(4) = domain%elems(ielem)%connectivity(corner_four)


                        !
                        ! Test each corner from try_face for match with bc_connectivity
                        !
                        includes_corner_one   = any( face_nodes == corner_indices(1) )
                        includes_corner_two   = any( face_nodes == corner_indices(2) )
                        includes_corner_three = any( face_nodes == corner_indices(3) )
                        includes_corner_four  = any( face_nodes == corner_indices(4) )

                        
                        !
                        ! If try_face corners are all in bc_connectivity:
                        !   - done
                        !   - set iface = try_face
                        !
                        face_matches_boundary = ( includes_corner_one   .and. &
                                                  includes_corner_two   .and. &
                                                  includes_corner_three .and. &
                                                  includes_corner_four )

                        ! Early exit condition
                        if (face_matches_boundary) then
                            iface=try_face
                            exit
                        end if

                    end do ! find iface


                    !
                    ! Check face was detected
                    !
                    user_msg = "bc%init: Could not determine element face associated with the boundary."
                    if (iface == NO_FACE) call chidg_signal(FATAL,user_msg)


                    !
                    ! Add domain, element, face index. Get face_ID, index of where the face exists in the bc_patch
                    !
                    face_ID = self%patch(patch_ID)%add_face(domain%elems(ielem)%ielement_g, &
                                                            domain%elems(ielem)%ielement_l, &
                                                            iface)

                    !
                    ! Inform domain face about:
                    !   - bc_ID    the boundary state group it is associated with
                    !   - group_ID the boundary patch group it is associated with
                    !   - patch_ID the boundary patch is is associated with
                    !   - face_ID  the patch face it is associated with
                    !   - ftype    the face type set to BOUNDARY
                    !
                    domain%faces(ielem,iface)%bc_ID    = bc_ID
                    domain%faces(ielem,iface)%group_ID = self%group_ID
                    domain%faces(ielem,iface)%patch_ID = patch_ID
                    domain%faces(ielem,iface)%face_ID  = face_ID

                    ! If associated with a boundary condition state group, reset face type.
                    if (bc_ID /= NO_ID) then
                        domain%faces(ielem,iface)%ftype = BOUNDARY
                    end if


                    ! End ielem search
                    exit


                end if ! all(nodes_matched)

            end do ! ielem

        end do !ielem_bc



    end subroutine add_bc_patch
    !******************************************************************************************













    !>  Extend self%patch allocation by one and return identifier for 
    !!  accesing the new object.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !!
    !---------------------------------------------------------------------
    function new_bc_patch(self) result(patch_ID)
        class(bc_patch_group_t),    intent(inout)   :: self

        integer(ik)                     :: patch_ID, ierr
        type(bc_patch_t),   allocatable :: temp_patches(:)



        !
        ! Resize array storage
        !
        allocate(temp_patches(self%npatches() + 1), stat=ierr)
        if (ierr /= 0) call AllocationError



        ! Copy previously initialized instances to new array. Be careful about pointers 
        ! components here. For example, a pointer from a face to an element would no 
        ! longer be valid in the new array.
        if (self%npatches() > 0) then
            temp_patches(1:size(self%patch)) = self%patch(1:size(self%patch))
        end if



        !
        ! Move resized temp allocation back to bc_patch_group container. 
        ! Be careful about pointer components here! Their location in 
        ! memory has changed.
        !
        call move_alloc(temp_patches,self%patch)
        


        !
        ! Set patch identifier of newly allocated patch that will be returned
        !
        patch_ID = self%npatches()



    end function new_bc_patch
    !*********************************************************************







    !>  Return the number of bc_patch's in the group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !----------------------------------------------------------------------
    function npatches(self) result(n)
        class(bc_patch_group_t),    intent(in)  :: self

        integer(ik) :: n


        if (allocated(self%patch)) then
            n = size(self%patch)
        else
            n = 0
        end if

    end function npatches
    !**********************************************************************






    !>  Return the number of faces in the group contributed from all patches.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/18/2017
    !!
    !----------------------------------------------------------------------
    function get_nfaces(self) result(n)
        class(bc_patch_group_t),    intent(in)  :: self

        integer(ik) :: n, patch_ID

        n = 0
        do patch_ID = 1,self%npatches()
            n = n + self%patch(patch_ID)%nfaces()
        end do

    end function get_nfaces
    !**********************************************************************





    !>  Return the processors that the boundary patch group is receiving
    !!  information from.
    !!
    !!  Accumulate from individual patches.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/10/2017
    !!
    !-----------------------------------------------------------------------
    function get_recv_procs(self) result(recv_procs_array)
        class(bc_patch_group_t),    intent(in)  :: self

        type(ivector_t)             :: recv_procs
        integer(ik),    allocatable :: recv_procs_patch(:), recv_procs_array(:)
        integer(ik)                 :: ipatch, iproc


        !
        ! Accumulate recv procs from each patch. Add uniquely
        !
        do ipatch = 1,self%npatches()
            recv_procs_patch = self%patch(ipatch)%get_recv_procs()

            do iproc = 1,size(recv_procs_patch)
                call recv_procs%push_back_unique(recv_procs_patch(iproc))
            end do

        end do !ipatch


        !
        ! Return as integer array
        !
        recv_procs_array = recv_procs%data()


    end function get_recv_procs
    !************************************************************************










    !>  Return the processors that the boundary patch group is sending
    !!  information to.
    !!
    !!  Currently, this pattern is the same as the receive pattern, so
    !!  we just call get_recv_procs.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/10/2017
    !!
    !-----------------------------------------------------------------------
    function get_send_procs(self) result(send_procs_array)
        class(bc_patch_group_t),    intent(in)  :: self

        integer(ik),    allocatable :: send_procs_array(:)

        !
        ! Return as integer array
        !
        send_procs_array = self%get_recv_procs()


    end function get_send_procs
    !************************************************************************









end module type_bc_patch_group
