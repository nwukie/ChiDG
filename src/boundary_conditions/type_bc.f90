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
    implicit none




    !>  Primary boundary condition container.
    !!
    !!  - contains a bc_family;             defining a general classification for the bc
    !!  - contains a bc_patch;              defining the bc geometry
    !!  - contains an array of bc_state's;  defining the solution state computed by the bc
    !!
    !!  Convention is that a given bc_t will exists in an array of bc_t's. Maybe something
    !!  like:
    !!
    !!      type(bc_t), allocatable :: bcs(:)
    !!
    !!  The BC_ID component of a boundary condition, is the location of the boundary condition
    !!  in such an array. In this way, the bc_t knows where it is located and can inform
    !!  other entities about where it is located.
    !!
    !!  bc_family may be:
    !!      'Wall', 'Inlet', 'Outlet', 'Symmetry', 'Periodic', 'Farfield', 'Scalar'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!  @note   Reorganized into bc patch, and bc operators
    !!
    !--------------------------------------------------------------------------------------------
    type, public :: bc_t

        ! Index of the boundary condition in a set. So a bc knows its location
        integer(ik)                             :: BC_ID

        ! Boundary condition family
        character(:),               allocatable :: bc_family

        ! Boundary condition patch
        type(bc_patch_t)                        :: bc_patch

        ! Boundary condition state
        class(bc_state_wrapper_t),  allocatable :: bc_state(:)

    contains

        procedure           :: init_bc              !< Main boundary condition initialization routine.

        procedure           :: init_bc_group        !< Initialize boundary condition group
        procedure           :: init_bc_patch        !< Initialize boundary condition patch
        procedure, private  :: init_bc_spec         !< Optional User-specialized initialization routine.
        procedure, private  :: init_bc_coupling     !< Initialize coupling interaction between bc elements.



        procedure   :: set_family           !< Set the boundary condition family.
        procedure   :: get_family           !< Return the boundary condition family.

        procedure   :: add_bc_state         !< Add a bc_state function to the boundary condition.

        procedure   :: get_ncoupled_elems   !< Return the number of elements coupled with a specified boundary element.


    end type bc_t
    !*********************************************************************************************





contains

    !> Initialize boundary condition
    !!
    !!  mesh and boundary_connectivity get passed in:
    !!      - Process the boundary_connectivity and search mesh for matching faces
    !!      - If a face matches, set it as a boundary face in mesh
    !!      - If a face matches, add it to the self%bc_patch
    !!      - Initialize the boundary coupling information
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!  @param[inout]   mesh            mesh_t object containing elements and faces
    !!  @param[in]      bconnectivity   Connectivity information for faces defining a boundary condition
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine init_bc(self,mesh,bconnectivity,bc_group,bc_groups,bc_wall,bc_inlet,bc_outlet,bc_symmetry,bc_farfield,bc_periodic)
        class(bc_t),                    intent(inout)           :: self
        type(mesh_t),                   intent(inout)           :: mesh
        type(boundary_connectivity_t),  intent(in)              :: bconnectivity
        character(*),                   intent(in)              :: bc_group
        type(bc_group_t),               intent(in)              :: bc_groups(:)
        class(bc_state_t),              intent(in), optional    :: bc_wall
        class(bc_state_t),              intent(in), optional    :: bc_inlet
        class(bc_state_t),              intent(in), optional    :: bc_outlet
        class(bc_state_t),              intent(in), optional    :: bc_symmetry
        class(bc_state_t),              intent(in), optional    :: bc_farfield
        class(bc_state_t),              intent(in), optional    :: bc_periodic



        !
        ! Boundary condition group initialization
        !
        call self%init_bc_group(mesh,bconnectivity,bc_group,bc_groups,bc_wall,bc_inlet,bc_outlet,bc_symmetry,bc_farfield,bc_periodic)



        !
        ! Boundary condition patch initialization
        !
        call self%init_bc_patch(mesh,bconnectivity,bc_group,bc_groups,bc_wall,bc_inlet,bc_outlet,bc_symmetry,bc_farfield,bc_periodic)


    end subroutine init_bc
    !**********************************************************************************************
    








    !>  Initialize boundary condition group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/19/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine init_bc_group(self,mesh,bconnectivity,bc_group,bc_groups,bc_wall,bc_inlet,bc_outlet,bc_symmetry,bc_farfield,bc_periodic)
        class(bc_t),                    intent(inout)           :: self
        type(mesh_t),                   intent(inout)           :: mesh
        type(boundary_connectivity_t),  intent(in)              :: bconnectivity
        character(*),                   intent(in)              :: bc_group
        type(bc_group_t),               intent(in)              :: bc_groups(:)
        class(bc_state_t),              intent(in), optional    :: bc_wall
        class(bc_state_t),              intent(in), optional    :: bc_inlet
        class(bc_state_t),              intent(in), optional    :: bc_outlet
        class(bc_state_t),              intent(in), optional    :: bc_symmetry
        class(bc_state_t),              intent(in), optional    :: bc_farfield
        class(bc_state_t),              intent(in), optional    :: bc_periodic


        character(:),       allocatable     :: user_msg
        class(bc_state_t),  allocatable     :: bc_state
        integer(ik)                         :: istate, igroup, ierr
        logical                             :: group_found, group_set




        !
        ! Find the correct bc_group in bc_groups(:)
        !
        group_set = .false.
        do igroup = 1,size(bc_groups)

            group_found = (trim(bc_group) == trim(bc_groups(igroup)%name) )
            if (group_found .and. (.not. group_set)) then


                !
                ! Set Family
                !
                call self%set_family(bc_groups(igroup)%family)

                
                !
                ! Set default boundary condition states if they were pass in:
                !
                if ( present(bc_wall) .and. (trim(bc_groups(igroup)%family) == 'Wall') ) then
                    if (allocated(bc_state)) deallocate(bc_state)
                    allocate(bc_state, source=bc_wall, stat=ierr)
                    call self%add_bc_state(bc_state)

                else if ( present(bc_inlet) .and. (trim(bc_groups(igroup)%family) == 'Inlet') ) then
                    if (allocated(bc_state)) deallocate(bc_state)
                    allocate(bc_state, source=bc_inlet, stat=ierr)
                    call self%add_bc_state(bc_state)

                else if ( present(bc_outlet) .and. (trim(bc_groups(igroup)%family) == 'Outlet') ) then
                    if (allocated(bc_state)) deallocate(bc_state)
                    allocate(bc_state, source=bc_outlet, stat=ierr)
                    call self%add_bc_state(bc_state)

                else if ( present(bc_symmetry) .and. (trim(bc_groups(igroup)%family) == 'Symmetry') ) then
                    if (allocated(bc_state)) deallocate(bc_state)
                    allocate(bc_state, source=bc_symmetry, stat=ierr)
                    call self%add_bc_state(bc_state)

                else if ( present(bc_farfield) .and. (trim(bc_groups(igroup)%family) == 'Farfield') ) then
                    if (allocated(bc_state)) deallocate(bc_state)
                    allocate(bc_state, source=bc_farfield, stat=ierr)
                    call self%add_bc_state(bc_state)
                else if ( present(bc_periodic) .and. (trim(bc_groups(igroup)%family) == 'Periodic') ) then
                    if (allocated(bc_state)) deallocate(bc_state)
                    allocate(bc_state, source=bc_periodic, stat=ierr)
                    call self%add_bc_state(bc_state)


                !
                ! If no default boundary condition was set for the group, add the states from the file:
                !
                else

                    ! Add all bc_states in the group to the boundary condition
                    do istate = 1,bc_groups(igroup)%bc_states%size()

                        ! Get boundary condition state
                        if (allocated(bc_state)) deallocate(bc_state)
                        allocate(bc_state, source=bc_groups(igroup)%bc_states%at(istate), stat=ierr)
                        if (ierr /= 0) call AllocationError

                        ! Add boundary condition state
                        call self%add_bc_state(bc_state)

                    end do !istate

                end if

                group_set = .true.
            end if

        end do




        user_msg = "chidg_data%add_bc: It looks like we didn't find a boundary state group that &
                    matches with the string indicated in a boundary patch. Make sure that a &
                    boundary state group with the correct name exists. Also make sure that the name &
                    set on the boundary patch corresponds to one of the boundary state groups that exists."
        if ((.not. group_set) .and. ((trim(bc_group) /= 'empty') .and. (trim(bc_group) /= 'Empty')) ) &
            call chidg_signal_one(FATAL,user_msg,trim(bc_group))


        if ( (trim(bc_group) == 'empty') .or. &
             (trim(bc_group) == 'Empty') ) call self%set_family('Empty')


    end subroutine init_bc_group
    !***********************************************************************************************












    !>  Initialize boundary condition patch.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/19/2016
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine init_bc_patch(self,mesh,bconnectivity,bc_group,bc_groups,bc_wall,bc_inlet,bc_outlet,bc_symmetry,bc_farfield,bc_periodic)
        class(bc_t),                    intent(inout)           :: self
        type(mesh_t),                   intent(inout)           :: mesh
        type(boundary_connectivity_t),  intent(in)              :: bconnectivity
        character(*),                   intent(in)              :: bc_group
        type(bc_group_t),               intent(in)              :: bc_groups(:)
        class(bc_state_t),              intent(in), optional    :: bc_wall
        class(bc_state_t),              intent(in), optional    :: bc_inlet
        class(bc_state_t),              intent(in), optional    :: bc_outlet
        class(bc_state_t),              intent(in), optional    :: bc_symmetry
        class(bc_state_t),              intent(in), optional    :: bc_farfield
        class(bc_state_t),              intent(in), optional    :: bc_periodic


        type(point_t)                   :: pnt, point_one, point_two, point_three
        character(len=:),   allocatable :: bcname, user_msg
        real(rk)                        :: time, x, y, z
        integer(ik)                     :: nelem_xi, nelem_eta, nelem_zeta, nelem_bc, ielem_bc,         & 
                                           xi_begin, eta_begin, zeta_begin, xi_end, eta_end, zeta_end,  & 
                                           ixi, ieta, izeta, ierr, ielem, ielem_test, nface_nodes,      &
                                           iface, inode, i, nfaces_bc, iface_bc, BC_face, ncoupled_elements

        logical,        allocatable :: node_matched(:), xi_face, eta_face, zeta_face
        integer(ik),    allocatable :: element_nodes(:)
        integer(ik)                 :: face_node


        !
        ! Get number of elements/faces associated with boundary condition.
        !
        nelem_bc = bconnectivity%get_nfaces()


        !
        ! Loop through each face in bc connectivity and call initialization for each face in local mesh
        !
        ! Find owner element, determine iface
        !
        nelem_bc = bconnectivity%get_nfaces()
        iface_bc = 0
        do ielem_bc = 1,nelem_bc

            nface_nodes = size(bconnectivity%data(ielem_bc)%data)
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
                    face_node     = bconnectivity%data(ielem_bc)%get_face_node(inode)
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
                    x = mesh%nodes(bconnectivity%data(ielem_bc)%data(1))%c1_
                    y = mesh%nodes(bconnectivity%data(ielem_bc)%data(1))%c2_
                    z = mesh%nodes(bconnectivity%data(ielem_bc)%data(1))%c3_
                    point_one = mesh%elems(ielem)%computational_point(x,y,z)
                    

                    x = mesh%nodes(bconnectivity%data(ielem_bc)%data(2))%c1_
                    y = mesh%nodes(bconnectivity%data(ielem_bc)%data(2))%c2_
                    z = mesh%nodes(bconnectivity%data(ielem_bc)%data(2))%c3_
                    point_two = mesh%elems(ielem)%computational_point(x,y,z)



                    x = mesh%nodes(bconnectivity%data(ielem_bc)%data(nface_nodes))%c1_
                    y = mesh%nodes(bconnectivity%data(ielem_bc)%data(nface_nodes))%c2_
                    z = mesh%nodes(bconnectivity%data(ielem_bc)%data(nface_nodes))%c3_
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
                    ! Add domain, element, face index. Get BC_face, index of where the face exists in the bc_patch
                    !
                    BC_face = self%bc_patch%add_face(mesh%idomain_l,ielem,iface)


                    !
                    ! Inform mesh face about BC_ID it is associated with and the location, BC_face, in the boundary condition.
                    !
                    mesh%faces(ielem,iface)%BC_ID   = self%BC_ID
                    mesh%faces(ielem,iface)%BC_face = BC_face




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





        !
        ! Call user-specialized boundary condition initialization
        !
        call self%init_bc_spec(mesh)


        !
        ! Call user-specialized boundary coupling initialization
        !
        call self%init_bc_coupling(mesh,self%bc_patch)




        !
        ! Push ncoupled elements back to mesh face
        !
        do iface_bc = 1,self%bc_patch%nfaces()

            !idom  = self%bc_patch%idomain_l(iface_bc)
            ielem = self%bc_patch%ielement_l(iface_bc)
            iface = self%bc_patch%iface(iface_bc)

            ncoupled_elements = self%bc_patch%coupled_elements(iface_bc)%size()

            mesh%faces(ielem,iface)%BC_ndepend = ncoupled_elements

        end do !iface_bc





    end subroutine init_bc_patch
    !**********************************************************************************************












    !>  Default specialized initialization procedure. This is called from the base bc%init procedure
    !!  and can be overwritten by derived types to implement specialized initiailization details.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !!  @param[inout]   mesh        mesh_t object containing elements and faces
    !!  @param[in]      iface       block face index to which the boundary condition is being applied
    !!
    !----------------------------------------------------------------------------------------------
    subroutine init_bc_spec(self,mesh)
        class(bc_t),    intent(inout)   :: self
        type(mesh_t),   intent(inout)   :: mesh




    end subroutine init_bc_spec
    !**********************************************************************************************









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
    !----------------------------------------------------------------------------------------------
    subroutine init_bc_coupling(self,mesh,bc_patch)
        class(bc_t),        intent(inout)   :: self
        type(mesh_t),       intent(in)      :: mesh
        type(bc_patch_t),   intent(inout)   :: bc_patch

        integer(ik) :: iop

        !
        ! Have bc_operators initialize the boundary condition coupling
        !
        if (allocated(self%bc_state)) then
            do iop = 1,size(self%bc_state)

                call self%bc_state(iop)%state%init_bc_coupling(mesh,bc_patch)

            end do !iop
        end if

    end subroutine init_bc_coupling
    !**********************************************************************************************









    !>  Function for returning the number of elements coupled with a specified boundary element.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/16/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------
    function get_ncoupled_elems(self,ielem) result(ncoupled_elems)
        class(bc_t),    intent(inout)   :: self
        integer(ik),    intent(in)      :: ielem

        integer(ik) :: ncoupled_elems


        ncoupled_elems = self%bc_patch%coupled_elements(ielem)%size()


    end function get_ncoupled_elems
    !**********************************************************************************************












    !>  Add a bc_state function to the boundary condition.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine add_bc_state(self,bc_state)
        class(bc_t),        intent(inout)   :: self
        class(bc_state_t),  intent(in)      :: bc_state

        integer(ik) :: iop, ierr
        class(bc_state_wrapper_t),   allocatable :: temp(:) 


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
        ! Allocate new state to end
        !
        allocate(temp(size(temp))%state, source=bc_state, stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Move temp allocation to bc
        !
        call move_alloc(temp, self%bc_state)


    end subroutine add_bc_state
    !**************************************************************************************************














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
    !---------------------------------------------------------------------------------------------------
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
    !***************************************************************************************************







    !>  Get the bc_family.
    !!
    !!  bc_family may be:
    !!      'Wall', 'Inlet', 'Outlet', 'Symmetry', 'Periodic', 'Farfield', 'Scalar'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/4/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------------
    function get_family(self) result(family)
        class(bc_t),        intent(inout)   :: self

        character(:),   allocatable :: family

        character(:),   allocatable :: user_msg


        if (allocated(self%bc_family)) then
            family = self%bc_family
        else
            user_msg = "bc%get_family: It looks like the boundary condition family was never set.&
                        Make sure bc%set_family gets called in the boundary condition initialization&
                        routine"
            call chidg_signal(FATAL,user_msg)
        end if

    end function get_family
    !***************************************************************************************************





end module type_bc
