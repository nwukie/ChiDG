module type_bc
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, &
                                          BOUNDARY, ORPHAN, ZERO, ONE, TWO, NO_ID, NO_FACE, NFACES
    use mod_grid,                   only: FACE_CORNERS

    use type_bc_patch,              only: bc_patch_t
    use type_bc_group,              only: bc_group_t
    use type_bc_state,              only: bc_state_t
    use type_bc_state_wrapper,      only: bc_state_wrapper_t
    use type_mesh,              only: mesh_t
    use type_domain,                only: domain_t
    use type_point,                 only: point_t
    use type_boundary_connectivity, only: boundary_connectivity_t
    use mod_chidg_mpi,              only: IRANK, NRANK, ChiDG_COMM
    use mpi_f08,                    only: mpi_comm, MPI_Comm_split, MPI_Comm_size, &
                                          MPI_Allgather, MPI_LOGICAL, MPI_UNDEFINED
    implicit none




    !>  Primary boundary condition container.
    !!
    !!      : contains a bc_name;               bc_t takes the name of the bc_group initialized
    !!      : contains a bc_family;             defining a general classification for the bc
    !!      : contains an array of bc_patch's;  defining the bc geometry. potentially across domains
    !!      : contains an array of bc_state's;  defining the solution state computed by the bc
    !!
    !!  Convention is that a given bc_t will exists in an array of bc_t's. Maybe something
    !!  like:
    !!
    !!      type(bc_t), allocatable :: bcs(:)
    !!
    !!  The bc_ID component of a boundary condition, is the location of the boundary condition
    !!  in such an array. In this way, the bc_t knows where it is located and can inform
    !!  other entities about where it is located.
    !!
    !!  If some other component has a bc_ID, they can access the particular boundary condition
    !!  associated with that identifier with:
    !!
    !!      bcs(bc_ID)
    !!
    !!
    !!  bc_family may be:
    !!      'Wall', 'Inlet', 'Outlet', 'Symmetry', 'Periodic', 'Farfield', 'Scalar'
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!  @date   8/30/2016 (AFRL)    Reorganized into bc_patch and bc_state
    !!  @date   2/28/2017           Extended to include multiple patches from different domains
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public :: bc_t

        ! Index of the boundary condition in a set. So a bc knows its location
        integer(ik)                             :: bc_ID

        ! MPI communicator that gets setup for parallel bc's
        type(mpi_comm)                          :: bc_COMM

        ! Boundary condition family
        character(:),               allocatable :: bc_name
        character(:),               allocatable :: bc_family

        ! Boundary condition patches
        type(bc_patch_t),           allocatable :: bc_patch(:)

        ! Boundary condition function group
        class(bc_state_wrapper_t),  allocatable :: bc_state(:)

    contains

        ! bc component addition routines
        procedure   :: init_bc_group            ! Initialize boundary condition group
        procedure   :: init_bc_patch            ! Initialize boundary condition patch

        ! infrastructure initialization routine
        procedure   :: init_bc_comm             ! Setup parallel boundary condition infrastructure
        procedure   :: init_bc_specialized      ! Optional User-specialized initialization routine.
        procedure   :: init_bc_coupling         ! Initialize coupling interaction between bc elements.
        procedure   :: propagate_bc_coupling    ! Propagate coupling information to mesh



        procedure   :: new_bc_state         ! Add a bc_state instance to the boundary condition.
        procedure   :: new_bc_patch         ! Add a bc_patch instance to the boundary condition.

        procedure   :: set_name             ! Set the boundary condition name.
        procedure   :: get_name             ! Return the boundary condition name.
        procedure   :: set_family           ! Set the boundary condition family.
        procedure   :: get_family           ! Return the boundary condition family.

        procedure   :: nfaces => get_nfaces ! Return total number of faces in 
        procedure   :: ncoupled_elements    ! Return the number of elements coupled with a specified boundary element.


    end type bc_t
    !******************************************************************************************





contains



    !>  Initialize array of boundary condition state functions, bc_state_t, from an 
    !!  incoming boundary condition group, bc_group_t. Also sets the Family type of the 
    !!  boundary condition.
    !!
    !!  Also gives the option to override particular boundary condition families. This is
    !!  provided through the optional bc_state_t functions that may be passed in. If these
    !!  are present, they will override any corresponding bc_states of the same family
    !!  in the bc_group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/19/2016
    !!  @date   2/27/2017   modifications for more general boundary conditions.
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_bc_group(self,bc_group,bc_wall,bc_inlet,bc_outlet,bc_symmetry,bc_farfield,bc_periodic,bc_scalar)
        class(bc_t),        intent(inout)           :: self
        type(bc_group_t),   intent(in)              :: bc_group
        class(bc_state_t),  intent(in), optional    :: bc_wall
        class(bc_state_t),  intent(in), optional    :: bc_inlet
        class(bc_state_t),  intent(in), optional    :: bc_outlet
        class(bc_state_t),  intent(in), optional    :: bc_symmetry
        class(bc_state_t),  intent(in), optional    :: bc_farfield
        class(bc_state_t),  intent(in), optional    :: bc_periodic
        class(bc_state_t),  intent(in), optional    :: bc_scalar

        integer(ik) :: istate, ierr, state_ID



        !
        ! Set Name/Family
        !
        call self%set_name(bc_group%name)
        call self%set_family(bc_group%family)

        
        !
        ! Set override boundary condition states if they were passed in:
        !
        if ( present(bc_wall) .and. (trim(bc_group%family) == 'Wall') ) then
            state_ID = self%new_bc_state()
            allocate(self%bc_state(state_ID)%state, source=bc_wall, stat=ierr)

        else if ( present(bc_inlet) .and. (trim(bc_group%family) == 'Inlet') ) then
            state_ID = self%new_bc_state()
            allocate(self%bc_state(state_ID)%state, source=bc_inlet, stat=ierr)

        else if ( present(bc_outlet) .and. (trim(bc_group%family) == 'Outlet') ) then
            state_ID = self%new_bc_state()
            allocate(self%bc_state(state_ID)%state, source=bc_outlet, stat=ierr)

        else if ( present(bc_symmetry) .and. (trim(bc_group%family) == 'Symmetry') ) then
            state_ID = self%new_bc_state()
            allocate(self%bc_state(state_ID)%state, source=bc_symmetry, stat=ierr)

        else if ( present(bc_farfield) .and. (trim(bc_group%family) == 'Farfield') ) then
            state_ID = self%new_bc_state()
            allocate(self%bc_state(state_ID)%state, source=bc_farfield, stat=ierr)

        else if ( present(bc_periodic) .and. (trim(bc_group%family) == 'Periodic') ) then
            state_ID = self%new_bc_state()
            allocate(self%bc_state(state_ID)%state, source=bc_periodic, stat=ierr)

        else if ( present(bc_scalar) .and. (trim(bc_group%family) == 'Scalar') ) then
            state_ID = self%new_bc_state()
            allocate(self%bc_state(state_ID)%state, source=bc_scalar, stat=ierr)


        !
        ! If no default boundary condition was set for the group, add the states from the file:
        !
        else

            ! Add all bc_states in the group to the boundary condition
            do istate = 1,bc_group%bc_states%size()

                ! Add boundary condition state
                state_ID = self%new_bc_state()
                allocate(self%bc_state(state_ID)%state, source=bc_group%bc_states%at(istate), stat=ierr)
                if (ierr /= 0) call AllocationError


            end do !istate

        end if


        !
        ! Catch any allocation error from overriding boundary condition group
        !
        if (ierr /= 0) call AllocationError


    end subroutine init_bc_group
    !******************************************************************************************








    !>  Initialize boundary condition patch.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/19/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_bc_patch(self,domain,bc_connectivity)
        class(bc_t),                    intent(inout)           :: self
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


        !
        ! Get number of elements/faces associated with boundary condition.
        !
        nelem_bc = bc_connectivity%get_nfaces()


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
            !do ielem = 1,mesh%nelem
            do ielem = 1,domain%nelem

                
                !
                ! Get nodes from element being tested
                !
                element_nodes = domain%elems(ielem)%connectivity%get_element_nodes()


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
                    ! create a new bc_patch if one hasn't already been created:
                    !
                    if (patch_ID == NO_ID) patch_ID = self%new_bc_patch()


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
                        corner_indices(1) = domain%elems(ielem)%connectivity%get_element_node(corner_one)
                        corner_indices(2) = domain%elems(ielem)%connectivity%get_element_node(corner_two)
                        corner_indices(3) = domain%elems(ielem)%connectivity%get_element_node(corner_three)
                        corner_indices(4) = domain%elems(ielem)%connectivity%get_element_node(corner_four)


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
                    face_ID = self%bc_patch(patch_ID)%add_face(domain%elems(ielem)%idomain_g,     &
                                                                  domain%elems(ielem)%idomain_l,     &
                                                                  domain%elems(ielem)%ielement_g,    &
                                                                  domain%elems(ielem)%ielement_l,    &
                                                                  iface)


                    !
                    ! Inform domain face about bc_ID it is associated with, patch_ID it belongs to, and the location, face_ID, in the bc_patch.
                    !
                    domain%faces(ielem,iface)%bc_ID    = self%bc_ID
                    domain%faces(ielem,iface)%patch_ID = patch_ID
                    domain%faces(ielem,iface)%face_ID  = face_ID




                    !
                    ! Set face type - 'ftype'
                    !
                    if ( self%get_family() == 'Periodic' ) then

                        ! Set to ORPHAN face so it will be recognized as chimera in the detection process.
                        domain%faces(ielem,iface)%ftype = ORPHAN

                        ! time, pnt do nothing here, but interface for function requires them.
                        domain%faces(ielem,iface)%periodic_offset  = .true.
                        domain%faces(ielem,iface)%chimera_offset_1 = self%bc_state(1)%state%bcproperties%compute('Offset-1',time,pnt)
                        domain%faces(ielem,iface)%chimera_offset_2 = self%bc_state(1)%state%bcproperties%compute('Offset-2',time,pnt)
                        domain%faces(ielem,iface)%chimera_offset_3 = self%bc_state(1)%state%bcproperties%compute('Offset-3',time,pnt)

                    else if ( allocated(self%bc_state) .and. (.not. self%get_family() == 'Periodic') ) then
                        domain%faces(ielem,iface)%ftype = BOUNDARY

                    else
                        domain%faces(ielem,iface)%ftype = ORPHAN

                    end if



                    ! End search
                    exit


                end if

            end do ! ielem

        end do !ielem_bc



    end subroutine init_bc_patch
    !******************************************************************************************












    !>  Initialize parallel boundary condition communicators.
    !!
    !!  When init_bc_comm is called, it is expected that it is also being called on all 
    !!  other processors for the same bc_t. Convention is that all processors have all 
    !!  the same boundary condition groups added. However, the patches added to these 
    !!  groups will be different on each processor. Some processors that contain no part
    !!  of the boundary condition geometry will have no patch data.
    !!
    !!  Example: Processors 0,1,2 both have all the boundary conditions. However, not all
    !!           the boundary condition objects on each processor have geometry associated
    !!           with them.
    !!
    !!      Proc 0                  Proc 1                  Proc 2              bc_COMM
    !!      ------                  ------                  ------              -------
    !!
    !!      bc_one                  bc_one                  bc_one               [0,2]
    !!        bc_patch(1)             bc_patch(empty)         bc_patch(1)
    !!
    !!
    !!      bc_two                  bc_two                  bc_two               [1,2]
    !!        bc_patch(empty)         bc_patch(1)             bc_patch(1)
    !!
    !!
    !!  The goal of init_bc_comm is for each boundary condition (ex. bc_one, bc_two) to know
    !!  what other processors have geometry allocated for the same boundary condition.
    !!  In this way, if some boundary condition implementation required interaction between
    !!  elements, boundary conditions could talk directly to the processors that also contain
    !!  portions of a boundary condition patch.
    !!
    !!
    !!  init_bc_comm creats an MPI communicator, bc%bc_COMM, that includes the processors
    !!  with bc_patch data that has been added, indicating they contain part of the boundary.
    !!  For those processors that do not contain geometry for a particular bc_t, the bc%bc_COMM
    !!  component is set to MPI_COMM_NULL and should not be used because it has no part in 
    !!  communication for that boundary condition.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/1/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_bc_comm(self,mesh)
        class(bc_t),    intent(inout)   :: self
        type(mesh_t),   intent(inout)   :: mesh

        logical                     :: irank_has_geometry, ranks_have_geometry(NRANK)
        integer(ik)                 :: ierr, color
        character(:),   allocatable :: user_msg



        !
        ! Check if current processor contains geometry associated with the bc_t
        !
        irank_has_geometry = allocated(self%bc_patch)



        !
        ! Send this single information to all and receive back ordered information from all
        !
        call MPI_Allgather(irank_has_geometry,1,MPI_LOGICAL,ranks_have_geometry,1,MPI_LOGICAL,ChiDG_COMM,ierr)
        user_msg = "bc%init_bc_comm: Error in collective MPI_Allgather for determining which &
                    MPI ranks contain portions of a boundary condition bc_patch."
        if (ierr /= 0) call chidg_signal(FATAL,user_msg)


        !
        ! Create a new MPI communicator for the current boundary condition 
        ! that includes only those processors with bc_patch data that
        ! has been allocated; indicating they contain a portion of the 
        ! bc geometry.
        !
        if (irank_has_geometry) then
            color = 1
        else
            color = MPI_UNDEFINED
        end if


        call MPI_Comm_split(ChiDG_COMM, color, IRANK, self%bc_COMM, ierr)
        user_msg = "bc%init_bc_comm: Error in collective MPI_Comm_split when trying &
                    to create a communicator for exchanging boundary condition data &
                    between processors."
        if (ierr /= 0) call chidg_signal(FATAL,user_msg)



    end subroutine init_bc_comm
    !******************************************************************************************









    !>  Call init_bc_specialized for all bc_state_t objects that have been attached to the 
    !!  boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !!  @param[inout]   mesh        mesh_t object containing elements and faces
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_bc_specialized(self,mesh)
        class(bc_t),    intent(inout)   :: self
        type(mesh_t),   intent(inout)   :: mesh

        integer(ik) :: iop


!        !
!        ! Have bc_operators initialize the boundary condition coupling
!        !
!        if (allocated(self%bc_state)) then
!            if (allocated(self%bc_patch)) then
!
!                do iop = 1,size(self%bc_state)
!                    call self%bc_state(iop)%state%init_bc_specialized(mesh,self%bc_patch, self%bc_COMM)
!                end do !iop
!
!            end if !bc_patch
!        end if !bc_state


    end subroutine init_bc_specialized
    !******************************************************************************************









    !>  Default boundary coupling initialization routine. 
    !!
    !!  Default initializes coupling for a given element to just itself and no coupling with 
    !!  other elements on the boundary. For a boundary condition that is coupled across the face
    !!  this routine can be overwritten to set the coupling information specific to the boundary 
    !!  condition.
    !!  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/16/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_bc_coupling(self,mesh)
        class(bc_t),        intent(inout)   :: self
        type(mesh_t),   intent(in)      :: mesh

        integer(ik) :: iop

!        !
!        ! Have bc_operators initialize the boundary condition coupling
!        !
!        if (allocated(self%bc_state)) then
!            if (allocated(self%bc_patch)) then
!
!                do iop = 1,size(self%bc_state)
!                    call self%bc_state(iop)%state%init_bc_coupling(mesh,self%bc_patch)
!                end do !iop
!
!            end if !bc_patch
!        end if !bc_state

    end subroutine init_bc_coupling
    !******************************************************************************************





    !>  Propagate boundary condition coupling information to mesh.
    !!
    !!  This informs mesh boundary faces how many element coupling dependencies they have 
    !!  and in-turn informs the infrastructure how many times the boundary condition
    !!  state functions need computed.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/28/2017
    !!
    !-----------------------------------------------------------------------------------------
    subroutine propagate_bc_coupling(self,mesh)
        class(bc_t),    intent(in)      :: self
        type(mesh_t),   intent(inout)   :: mesh

        integer(ik) :: ipatch, iface_bc, idom, ielem, iface


!        !
!        ! set ncoupled elements back to mesh face
!        !
!        if (allocated(self%bc_patch)) then
!
!            do ipatch = 1,size(self%bc_patch)
!                do iface_bc = 1,self%bc_patch(ipatch)%nfaces()
!
!                    idom  = self%bc_patch(ipatch)%idomain_l(iface_bc)
!                    ielem = self%bc_patch(ipatch)%ielement_l(iface_bc)
!                    iface = self%bc_patch(ipatch)%iface(iface_bc)
!
!                    mesh%domain(idom)%faces(ielem,iface)%bc_ndepend = self%bc_patch(ipatch)%ncoupled_elements(iface_bc)
!
!                end do !iface_bc
!            end do !ipatch
!
!        end if !bc_patch



    end subroutine propagate_bc_coupling
    !*****************************************************************************************








    !>  Function for returning the number of elements coupled with a specified boundary element.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/16/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    function ncoupled_elements(self,ipatch,iface) result(ncoupled_elems)
        class(bc_t),    intent(inout)   :: self
        integer(ik),    intent(in)      :: ipatch
        integer(ik),    intent(in)      :: iface

        integer(ik) :: ncoupled_elems


        ncoupled_elems = self%bc_patch(ipatch)%ncoupled_elements(iface)


    end function ncoupled_elements
    !******************************************************************************************












    !>  Add a bc_state function to the boundary condition.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    function new_bc_state(self) result(state_ID)
        class(bc_t),    intent(inout)   :: self

        integer(ik)                             :: iop, ierr, state_ID
        class(bc_state_wrapper_t),  allocatable :: temp(:) 


        !
        ! Allocate temp storage for (size+1), copy current states to temp
        !        
        if (allocated(self%bc_state)) then

            allocate(temp(size(self%bc_state) + 1), stat=ierr)
            if (ierr /= 0) call AllocationError

            ! Copy previously added states to temp
            do iop = 1,size(self%bc_state)
                allocate(temp(iop)%state, source=self%bc_state(iop)%state, stat=ierr)
                if (ierr /= 0) call AllocationError
            end do

        else

            allocate(temp(1), stat=ierr)
            if (ierr /= 0) call AllocationError

        end if


        !
        ! Move temp allocation to bc
        !
        call move_alloc(temp, self%bc_state)


        !
        ! Set state ID to last entry, which is the new one.
        !
        state_ID = size(self%bc_state)


    end function new_bc_state
    !******************************************************************************************








    !>  Add a bc_patch instance to the boundary condition. Return its PATCH_ID index of 
    !!  its location as bc%bc_patch(PATCH_ID).
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/30/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    function new_bc_patch(self) result(patch_ID)
        class(bc_t),        intent(inout)   :: self

        integer(ik)                     :: iop, ierr, patch_ID
        type(bc_patch_t),   allocatable :: temp(:) 


        !
        ! Allocate temp storage for (size+1), copy current states to temp
        !        
        if (allocated(self%bc_patch)) then

            allocate(temp(size(self%bc_patch) + 1), stat=ierr)
            if (ierr /= 0) call AllocationError

            temp(1:size(self%bc_patch)) = self%bc_patch(1:size(self%bc_patch))

        else

            allocate(temp(1), stat=ierr)
            if (ierr /= 0) call AllocationError

        end if

        
        !
        ! Move temp allocation to bc
        !
        call move_alloc(temp, self%bc_patch)


        !
        ! 1: Set patch ID to last entry, which is the new one.
        ! 2: Give new bc_patch its identifier.
        !
        patch_ID = size(self%bc_patch)
        self%bc_patch(patch_ID)%patch_ID = patch_ID


    end function new_bc_patch
    !******************************************************************************************











    !>  Set the bc_family.
    !!
    !!  bc_family may be:
    !!      'Wall', 'Inlet', 'Outlet', 'Symmetry', 'Periodic', 'Farfield', 'Scalar'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/4/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine set_family(self,family)
        class(bc_t),        intent(inout)   :: self
        character(*),       intent(in)      :: family

        character(:),   allocatable :: user_msg

        if ( (trim(family) == 'Wall'    )       .or. &
             (trim(family) == 'Inlet'   )       .or. &
             (trim(family) == 'Outlet'  )       .or. &
             (trim(family) == 'Symmetry')       .or. &
             (trim(family) == 'Periodic')       .or. &
             (trim(family) == 'Farfield')       .or. &
             (trim(family) == 'Scalar'  )       .or. &
             (trim(family) == 'Extrapolation')  .or. &
             (trim(family) == 'Empty'   ) ) then

            self%bc_family = family

        else
            user_msg = "bc%set_family: The string passed in to set the boundary condition family did &
                        not match any of valid boundary condition families. These include: 'Wall', &
                        'Inlet', 'Outlet', 'Symmetry', 'Periodic', 'Farfield', 'Scalar', 'Extrapolation'"
            call chidg_signal_one(FATAL,user_msg,family)
        end if

    end subroutine set_family
    !******************************************************************************************







    !>  Get the bc_family.
    !!
    !!  bc_family may be:
    !!      'Wall', 'Inlet', 'Outlet', 'Symmetry', 'Periodic', 'Farfield', 'Scalar'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/4/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    function get_family(self) result(family)
        class(bc_t),        intent(inout)   :: self

        character(:),   allocatable :: family, user_msg


        if (allocated(self%bc_family)) then
            family = self%bc_family
        else
            user_msg = "bc%get_family: It looks like the boundary condition family was never set.&
                        Make sure bc%set_family gets called in the boundary condition initialization&
                        routine"
            call chidg_signal(FATAL,user_msg)
        end if

    end function get_family
    !******************************************************************************************







    !>  Set the bc_name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   02/28/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine set_name(self,bc_name)
        class(bc_t),        intent(inout)   :: self
        character(*),       intent(in)      :: bc_name

        character(:),   allocatable :: user_msg


        self%bc_name = bc_name


    end subroutine set_name
    !******************************************************************************************





    !>  Return the bc_name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   02/28/2017
    !!
    !------------------------------------------------------------------------------------------
    function get_name(self) result(bc_name)
        class(bc_t),        intent(inout)   :: self

        character(:),   allocatable :: bc_name

        bc_name = self%bc_name

    end function get_name
    !******************************************************************************************




    !>  Return number of faces attached to the boundary condition. Includes contributions from
    !!  all bc_patches attached to the boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !------------------------------------------------------------------------------------------
    function get_nfaces(self) result(nfaces_)
        class(bc_t),    intent(in)  :: self

        integer(ik) :: ipatch, nfaces_

        nfaces_ = 0
        do ipatch = 1,size(self%bc_patch)

            nfaces_ = nfaces_ + self%bc_patch(ipatch)%nfaces()

        end do !ipatch

    end function get_nfaces
    !******************************************************************************************











end module type_bc
