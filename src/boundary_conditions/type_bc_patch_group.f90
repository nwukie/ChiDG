module type_bc_patch_group
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: NO_ID, BOUNDARY, ORPHAN
    use mod_grid,                   only: FACE_CORNERS
    use type_point,                 only: point_t
    use type_bc_patch,              only: bc_patch_t
    use type_domain,                only: domain_t
    use type_boundary_connectivity, only: boundary_connectivity_t
    implicit none




    !>  An object containing all the bc_patch_t instances for a single
    !!  boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !--------------------------------------------------------------------
    type, public :: bc_patch_group_t

        character(:),       allocatable :: name

        type(bc_patch_t),   allocatable :: patch(:)

        integer(ik)                     :: group_ID

    contains

        procedure           :: add_bc_patch
        procedure           :: npatches
        procedure,  private :: new_bc_patch

    end type bc_patch_group_t
    !********************************************************************







contains


    !>  Initialize boundary condition patch.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/19/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_bc_patch(self,domain,bc_connectivity)
        class(bc_patch_group_t),        intent(inout)           :: self
        type(domain_t),                 intent(inout)           :: domain
        type(boundary_connectivity_t),  intent(in)              :: bc_connectivity


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





!        !
!        ! Get number of elements/faces associated with boundary condition.
!        !
!        nelem_bc = bc_connectivity%get_nfaces()
!
!
!        !
!        ! Loop through each face in bc connectivity:
!        !
!        !   - Find owner element, determine iface
!        !
!        iface_bc = 0
!        patch_ID = NO_ID
!        do ielem_bc = 1,nelem_bc
!
!            !
!            ! Get nodes from face
!            !
!            face_nodes = bc_connectivity%data(ielem_bc)%get_face_nodes()
!
!
!            nface_nodes = size(face_nodes)
!            if ( allocated(node_matched) ) deallocate(node_matched)
!            allocate(node_matched(nface_nodes), stat=ierr)
!            if (ierr /= 0) call AllocationError
!
!            
!            !
!            ! Search for local element with common nodes
!            !
!            do ielem = 1,domain%nelem
!
!                
!                !
!                ! Get nodes from element being tested
!                !
!                element_nodes = domain%elems(ielem)%connectivity%get_element_nodes()
!
!
!                !
!                ! Loop through bc nodes and see if current element has them all
!                !
!                node_matched = .false.
!                do inode = 1,nface_nodes
!                    if ( any(element_nodes == face_nodes(inode) ) ) then
!                        node_matched(inode) = .true.
!                    end if
!                end do
!
!
!                !
!                ! Determine element mapping index
!                !
!                nterms_1d = 0
!                do while (nterms_1d*nterms_1d*nterms_1d < domain%elems(ielem)%nterms_c)
!                    nterms_1d = nterms_1d + 1
!                end do
!                mapping = nterms_1d - 1
!
!
!
!        
!                !
!                ! If all match: 
!                !   - determine face index, iface, corresponding to the element boundary
!                !   - set element/face indices in boundary condition list.
!                !
!                if ( all(node_matched) ) then
!                    iface_bc = iface_bc + 1
!
!                    !
!                    ! Since a geometry region was detected on the current processor,
!                    ! create a new bc_patch if one hasn't already been created.
!                    ! We want this creation statement in the inner loop here so that
!                    ! only those processor that contain faces on the boundary 
!                    ! actually allocate a patch.
!                    !
!                    if (patch_ID == NO_ID) patch_ID = self%new_bc_patch()
!
!
!                    !
!                    ! Detect element face index associated with the boundary condition
!                    !   Goal: determine iface
!                    !
!                    !   Loop through element faces until a match is found
!                    !
!                    iface = NO_FACE
!                    do try_face = 1,NFACES 
!
!                        !
!                        ! Get corner nodes for face, try_face
!                        !
!                        corner_one   = FACE_CORNERS(try_face,1,mapping)
!                        corner_two   = FACE_CORNERS(try_face,2,mapping)
!                        corner_three = FACE_CORNERS(try_face,3,mapping)
!                        corner_four  = FACE_CORNERS(try_face,4,mapping)
!
!
!                        !
!                        ! For the element, get the global indices of the corner nodes for face, try_face
!                        !
!                        corner_indices(1) = domain%elems(ielem)%connectivity%get_element_node(corner_one)
!                        corner_indices(2) = domain%elems(ielem)%connectivity%get_element_node(corner_two)
!                        corner_indices(3) = domain%elems(ielem)%connectivity%get_element_node(corner_three)
!                        corner_indices(4) = domain%elems(ielem)%connectivity%get_element_node(corner_four)
!
!
!                        !
!                        ! Test each corner from try_face for match with bc_connectivity
!                        !
!                        includes_corner_one   = any( face_nodes == corner_indices(1) )
!                        includes_corner_two   = any( face_nodes == corner_indices(2) )
!                        includes_corner_three = any( face_nodes == corner_indices(3) )
!                        includes_corner_four  = any( face_nodes == corner_indices(4) )
!
!                        
!                        !
!                        ! If try_face corners are all in bc_connectivity:
!                        !   - done
!                        !   - set iface = try_face
!                        !
!                        face_matches_boundary = ( includes_corner_one   .and. &
!                                                  includes_corner_two   .and. &
!                                                  includes_corner_three .and. &
!                                                  includes_corner_four )
!
!                        ! Early exit condition
!                        if (face_matches_boundary) then
!                            iface=try_face
!                            exit
!                        end if
!
!                    end do ! find iface
!
!
!                    !
!                    ! Check face was detected
!                    !
!                    user_msg = "bc%init: Could not determine element face associated with the boundary."
!                    if (iface == NO_FACE) call chidg_signal(FATAL,user_msg)
!
!
!                    !
!                    ! Add domain, element, face index. Get face_ID, index of where the face exists in the bc_patch
!                    !
!                    face_ID = self%patch(patch_ID)%add_face(domain%elems(ielem)%idomain_g,  &
!                                                            domain%elems(ielem)%idomain_l,  &
!                                                            domain%elems(ielem)%ielement_g, &
!                                                            domain%elems(ielem)%ielement_l, &
!                                                            iface)
!
!
!                    !
!                    ! Inform domain face about bc_ID it is associated with, patch_ID it belongs to, and the location, face_ID, in the bc_patch.
!                    !
!                    domain%faces(ielem,iface)%bc_ID      = self%bc_ID
!                    domain%faces(ielem,iface)%group_ID   = self%group_ID
!                    domain%faces(ielem,iface)%patch_ID   = patch_ID
!                    domain%faces(ielem,iface)%face_ID    = face_ID
!                    !domain%faces(ielem,iface)%patch_face = patch_face
!
!
!
!
!                    !
!                    ! Set face type - 'ftype'
!                    !
!                    if ( self%get_family() == 'Periodic' ) then
!
!                        ! Set to ORPHAN face so it will be recognized as chimera in the detection process.
!                        domain%faces(ielem,iface)%ftype = ORPHAN
!
!                        ! time, pnt do nothing here, but interface for function requires them.
!                        domain%faces(ielem,iface)%periodic_offset  = .true.
!                        domain%faces(ielem,iface)%chimera_offset_1 = self%bc_state(1)%state%bcproperties%compute('Offset-1',time,pnt)
!                        domain%faces(ielem,iface)%chimera_offset_2 = self%bc_state(1)%state%bcproperties%compute('Offset-2',time,pnt)
!                        domain%faces(ielem,iface)%chimera_offset_3 = self%bc_state(1)%state%bcproperties%compute('Offset-3',time,pnt)
!
!                    else if ( allocated(self%bc_state) .and. (.not. self%get_family() == 'Periodic') ) then
!                        domain%faces(ielem,iface)%ftype = BOUNDARY
!
!                    else
!                        domain%faces(ielem,iface)%ftype = ORPHAN
!
!                    end if
!
!
!
!                    ! End search
!                    exit
!
!
!                end if
!
!            end do ! ielem
!
!        end do !ielem_bc



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



        ! Copy previously initialized instances to new array. Be careful about pointers 
        ! components here. For example, a pointer from a face to an element would no 
        ! longer be valid in the new array.
        if (self%npatches() > 0) then
            temp_patches(1:size(self%patch)) = self%patch(1:size(self%patch))
        end if



        !
        ! Move resized temp allocation back to mesh container. 
        ! Be careful about pointer components here! Their location in memory has changed.
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




end module type_bc_patch_group
