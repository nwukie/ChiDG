module type_bc
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, &
                                          BOUNDARY, ORPHAN, ZERO, ONE, TWO, RKTOL

    use type_bc_patch,              only: bc_patch_t
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

        procedure   :: init_bc              !< Boundary condition initialization.
        procedure   :: init_bc_spec         !< Call specialized initialization routine.
        procedure   :: init_bc_coupling     !< Initialize book-keeping for coupling interaction between elements.

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
    subroutine init_bc(self,mesh,bconnectivity)
        class(bc_t),                    intent(inout)   :: self
        type(mesh_t),                   intent(inout)   :: mesh
        type(boundary_connectivity_t),  intent(in)      :: bconnectivity

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




!                    ! Set face type - 'ftype'
!                    bcname = self%get_name()
!                    if ( trim(bcname) == 'periodic' ) then
!                        !
!                        ! Set to ORPHAN face so it will be recognized as chimera in the detection process.
!                        !
!                        mesh%faces(ielem,iface)%ftype = ORPHAN
!
!                        !
!                        ! Set periodic offset from boundary condition to the face. To be used in detection of gq_donor.
!                        !
!                        if ( self%bcproperties%compute('type', time, pnt) == ONE ) then
!                            mesh%faces(ielem,iface)%periodic_type = 'cartesian'
!                        else if ( self%bcproperties%compute('type', time, pnt) == TWO ) then
!                            mesh%faces(ielem,iface)%periodic_type = 'cylindrical'
!                        end if
!
!                        ! time, pnt do nothing here, but interface for function requires them.
!                        mesh%faces(ielem,iface)%chimera_offset_x     = self%bcproperties%compute('offset_x',     time, pnt)
!                        mesh%faces(ielem,iface)%chimera_offset_y     = self%bcproperties%compute('offset_y',     time, pnt)
!                        mesh%faces(ielem,iface)%chimera_offset_z     = self%bcproperties%compute('offset_z',     time, pnt)
!                        mesh%faces(ielem,iface)%chimera_offset_theta = self%bcproperties%compute('offset_theta', time, pnt)
!
!                    else
!                        !
!                        ! Set face to boundary condition face
!                        !
!                        mesh%faces(ielem,iface)%ftype = BOUNDARY
!                    end if


                    if ( allocated(self%bc_state) ) then
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





    end subroutine init_bc
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

        if ( (trim(family) == 'Wall'    ) .or. &
             (trim(family) == 'Inlet'   ) .or. &
             (trim(family) == 'Outlet'  ) .or. &
             (trim(family) == 'Symmetry') .or. &
             (trim(family) == 'Periodic') .or. &
             (trim(family) == 'Farfield') .or. &
             (trim(family) == 'Scalar'  ) ) then

            self%bc_family = family

        else
            user_msg = "bc%set_family: The string passed in to set the boundary condition family did &
                        not match any of valid boundary condition families. These include: 'Wall', &
                        'Inlet', 'Outlet', 'Symmetry', 'Periodic', 'Farfield', 'Scalar'"
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
