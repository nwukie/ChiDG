module type_bc
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, &
                                          BOUNDARY, ORPHAN, ZERO, ONE, TWO, RKTOL

    use type_bc_patch,              only: bc_patch_t
    use type_bc_group,              only: bc_group_t
    use type_bc_state,              only: bc_state_t
    use type_bc_state_wrapper,      only: bc_state_wrapper_t
    use type_mesh,                  only: mesh_t
    use type_point,                 only: point_t
    use type_boundary_connectivity, only: boundary_connectivity_t
    use mpi_f08,                    only: mpi_comm
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

        procedure   :: init_bc_group            ! Initialize boundary condition group
        procedure   :: init_bc_patch            ! Initialize boundary condition patch
        procedure   :: init_bc_specialized      ! Optional User-specialized initialization routine.
        procedure   :: init_bc_coupling         ! Initialize coupling interaction between bc elements.
        procedure   :: propagate_bc_coupling    ! Propagate coupling information to mesh



        procedure   :: new_bc_state         ! Add a bc_state instance to the boundary condition.
        procedure   :: new_bc_patch         ! Add a bc_patch instance to the boundary condition.

        procedure   :: set_name             ! Set the boundary condition name.
        procedure   :: get_name             ! Return the boundary condition name.
        procedure   :: set_family           ! Set the boundary condition family.
        procedure   :: get_family           ! Return the boundary condition family.

        procedure   :: nfaces               ! Return total number of faces in 
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
    subroutine init_bc_patch(self,mesh,bc_connectivity)
        class(bc_t),                    intent(inout)           :: self
        type(mesh_t),                   intent(inout)           :: mesh
        type(boundary_connectivity_t),  intent(in)              :: bc_connectivity


        type(point_t)               :: pnt, point_one, point_two, point_three
        character(:),   allocatable :: bcname, user_msg
        real(rk)                    :: time, x, y, z
        integer(ik)                 :: nelem_xi, nelem_eta, nelem_zeta, nelem_bc, ielem_bc,         & 
                                       xi_begin, eta_begin, zeta_begin, xi_end, eta_end, zeta_end,  & 
                                       ixi, ieta, izeta, ierr, ielem, ielem_test, nface_nodes,      &
                                       iface, inode, i, nfaces_bc, iface_bc, patch_face, patch_ID,  &
                                       ncoupled_elements

        logical,        allocatable :: node_matched(:), xi_face, eta_face, zeta_face
        integer(ik),    allocatable :: element_nodes(:)
        integer(ik)                 :: face_node


        !
        ! Create a new bc_patch
        !
        patch_ID = self%new_bc_patch()


        !
        ! Get number of elements/faces associated with boundary condition.
        !
        nelem_bc = bc_connectivity%get_nfaces()


        !
        ! Loop through each face in bc connectivity and call initialization for each face in local mesh
        !
        ! Find owner element, determine iface
        !
        nelem_bc = bc_connectivity%get_nfaces()
        iface_bc = 0
        do ielem_bc = 1,nelem_bc

            nface_nodes = size(bc_connectivity%data(ielem_bc)%data)
            if ( allocated(node_matched) ) deallocate(node_matched)
            allocate(node_matched(nface_nodes), stat=ierr)
            if (ierr /= 0) call AllocationError

            
            !
            ! Search for local element with common nodes
            !
            do ielem = 1,mesh%nelem


                ! Loop through bc nodes and see if current element has them all
                node_matched = .false.
                do inode = 1,nface_nodes
                    element_nodes = mesh%elems(ielem)%connectivity%get_element_nodes()
                    face_node     = bc_connectivity%data(ielem_bc)%get_face_node(inode)
                    if ( any(element_nodes == face_node) ) then
                        node_matched(inode) = .true.
                    end if
                end do

        

                ! If all match, set element/face indices in boundary condition list.
                if ( all(node_matched) ) then
                    iface_bc = iface_bc + 1


                    !
                    ! Detect element face index associated with the boundary condition
                    !
                    ! Get xi,eta,zeta for three points, defining a face
                    !
                    x = mesh%nodes(bc_connectivity%data(ielem_bc)%data(1))%c1_
                    y = mesh%nodes(bc_connectivity%data(ielem_bc)%data(1))%c2_
                    z = mesh%nodes(bc_connectivity%data(ielem_bc)%data(1))%c3_
                    point_one = mesh%elems(ielem)%computational_point(x,y,z)
                    

                    x = mesh%nodes(bc_connectivity%data(ielem_bc)%data(2))%c1_
                    y = mesh%nodes(bc_connectivity%data(ielem_bc)%data(2))%c2_
                    z = mesh%nodes(bc_connectivity%data(ielem_bc)%data(2))%c3_
                    point_two = mesh%elems(ielem)%computational_point(x,y,z)



                    x = mesh%nodes(bc_connectivity%data(ielem_bc)%data(nface_nodes))%c1_
                    y = mesh%nodes(bc_connectivity%data(ielem_bc)%data(nface_nodes))%c2_
                    z = mesh%nodes(bc_connectivity%data(ielem_bc)%data(nface_nodes))%c3_
                    point_three = mesh%elems(ielem)%computational_point(x,y,z)

                    ! Check to make sure the computational point coordinates were all valid and found.
                    user_msg = "bc%init_bc: BC connectivity node not found on element face. Invalid point detected."
                    if ( (.not. point_one%valid()) .or. &
                         (.not. point_two%valid()) .or. &
                         (.not. point_three%valid()) ) call chidg_signal(FATAL,user_msg)


                    xi_face   = ( (abs(point_one%c1_ - point_two%c1_  ) < 1.e-5_rk)  .and. &
                                  (abs(point_one%c1_ - point_three%c1_) < 1.e-5_rk) )
                    eta_face  = ( (abs(point_one%c2_ - point_two%c2_  ) < 1.e-5_rk)  .and. &
                                  (abs(point_one%c2_ - point_three%c2_) < 1.e-5_rk) )
                    zeta_face = ( (abs(point_one%c3_ - point_two%c3_  ) < 1.e-5_rk)  .and. &
                                  (abs(point_one%c3_ - point_three%c3_) < 1.e-5_rk) )



                    ! Determine iface by value of xi,eta,zeta
                    if ( xi_face ) then
                        if (point_one%c1_ < ZERO) then
                            iface = XI_MIN
                        else
                            iface = XI_MAX
                        end if
                    else if ( eta_face ) then
                        if (point_one%c2_ < ZERO) then
                            iface = ETA_MIN
                        else
                            iface = ETA_MAX
                        end if
                    else if ( zeta_face ) then
                        if (point_one%c3_ < ZERO) then
                            iface = ZETA_MIN
                        else
                            iface = ZETA_MAX
                        end if
                    else
                        call chidg_signal(FATAL,"bc%init: could not determine element face associated with the boundary")
                    end if


                    !
                    ! Add domain, element, face index. Get patch_face, index of where the face exists in the bc_patch
                    !
                    patch_face = self%bc_patch(patch_ID)%add_face(mesh%elems(ielem)%idomain_g,     &
                                                                  mesh%elems(ielem)%idomain_l,     &
                                                                  mesh%elems(ielem)%ielement_g,    &
                                                                  mesh%elems(ielem)%ielement_l,    &
                                                                  iface)


                    !
                    ! Inform mesh face about bc_ID it is associated with, patch_ID it belongs to, and the location, patch_face, in the bc_patch.
                    !
                    mesh%faces(ielem,iface)%bc_ID      = self%bc_ID
                    mesh%faces(ielem,iface)%patch_ID   = patch_ID
                    mesh%faces(ielem,iface)%patch_face = patch_face




                    !
                    ! Set face type - 'ftype'
                    !
                    if ( self%get_family() == 'Periodic' ) then

                        ! Set to ORPHAN face so it will be recognized as chimera in the detection process.
                        mesh%faces(ielem,iface)%ftype = ORPHAN

                        ! time, pnt do nothing here, but interface for function requires them.
                        mesh%faces(ielem,iface)%periodic_offset  = .true.
                        mesh%faces(ielem,iface)%chimera_offset_1 = self%bc_state(1)%state%bcproperties%compute('Offset-1',time,pnt)
                        mesh%faces(ielem,iface)%chimera_offset_2 = self%bc_state(1)%state%bcproperties%compute('Offset-2',time,pnt)
                        mesh%faces(ielem,iface)%chimera_offset_3 = self%bc_state(1)%state%bcproperties%compute('Offset-3',time,pnt)

                    else if ( allocated(self%bc_state) .and. (.not. self%get_family() == 'Periodic') ) then
                        mesh%faces(ielem,iface)%ftype = BOUNDARY

                    else
                        mesh%faces(ielem,iface)%ftype = ORPHAN

                    end if



                    ! End search
                    exit


                end if

            end do ! ielem

        end do !ielem_bc



    end subroutine init_bc_patch
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
        type(mesh_t),   intent(inout)   :: mesh(:)

        integer(ik) :: iop


        !
        ! Have bc_operators initialize the boundary condition coupling
        !
        if (allocated(self%bc_state)) then
            if (allocated(self%bc_patch)) then

                do iop = 1,size(self%bc_state)
                    call self%bc_state(iop)%state%init_bc_specialized(mesh,self%bc_patch)
                end do !iop

            end if !bc_patch
        end if !bc_state


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
        type(mesh_t),       intent(in)      :: mesh(:)

        integer(ik) :: iop

        !
        ! Have bc_operators initialize the boundary condition coupling
        !
        if (allocated(self%bc_state)) then
            if (allocated(self%bc_patch)) then

                do iop = 1,size(self%bc_state)
                    call self%bc_state(iop)%state%init_bc_coupling(mesh,self%bc_patch)
                end do !iop

            end if !bc_patch
        end if !bc_state

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
        type(mesh_t),   intent(inout)   :: mesh(:)

        integer(ik) :: ipatch, iface_bc, idom, ielem, iface


        !
        ! set ncoupled elements back to mesh face
        !
        if (allocated(self%bc_patch)) then

            do ipatch = 1,size(self%bc_patch)
                do iface_bc = 1,self%bc_patch(ipatch)%nfaces()

                    idom  = self%bc_patch(ipatch)%idomain_l(iface_bc)
                    ielem = self%bc_patch(ipatch)%ielement_l(iface_bc)
                    iface = self%bc_patch(ipatch)%iface(iface_bc)

                    mesh(idom)%faces(ielem,iface)%bc_ndepend = self%bc_patch(ipatch)%ncoupled_elements(iface_bc)

                end do !iface_bc
            end do !ipatch

        end if !bc_patch



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
    function nfaces(self) result(nfaces_)
        class(bc_t),    intent(in)  :: self

        integer(ik) :: ipatch, nfaces_

        nfaces_ = 0
        do ipatch = 1,size(self%bc_patch)

            nfaces_ = nfaces_ + self%bc_patch(ipatch)%nfaces()

        end do !ipatch

    end function nfaces
    !******************************************************************************************











end module type_bc
