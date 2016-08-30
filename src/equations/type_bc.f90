module type_bc
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_chidg_mpi,              only: IRANK
    use mod_constants,              only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, &
                                          BOUNDARY, CHIMERA, ORPHAN, BC_BLK, ZERO, ONE, TWO, RKTOL

    use type_chidg_worker,          only: chidg_worker_t
    use type_equation_set,          only: equation_set_t
    use type_bc_patch,              only: bc_patch_t
    use type_bc_operator,           only: bc_operator_t
    use type_bc_operator_wrapper,   only: bc_operator_wrapper_t
    use type_mesh,                  only: mesh_t
    use type_point,                 only: point_t
    use type_ivector,               only: ivector_t
    use type_solverdata,            only: solverdata_t
    use type_properties,            only: properties_t
    use type_bcproperty_set,        only: bcproperty_set_t
    use type_boundary_connectivity, only: boundary_connectivity_t
    implicit none




    !> Abstract base-type for boundary conditions
    !!  - contains a list of associated element indices
    !!  - contains a list of face indices
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!  @note   Reorganized into bc patch, and bc operators
    !!
    !--------------------------------------------------------------------------------------------
    type, public, extends(equation_set_t) :: bc_t

        ! Boundary condition patch
        type(bc_patch_t)                            :: bc_patch

        ! Boundary condition operators
        class(bc_operator_wrapper_t),   allocatable :: bc_advective_operator(:)

    contains

        procedure   :: init_bc              !< Boundary condition initialization
        procedure   :: init_bc_spec         !< Call specialized initialization routine
        procedure   :: init_bc_coupling     !< Initialize book-keeping for coupling interaction between elements.

        procedure   :: add_operator

        procedure   :: compute_bc_operators     !< Apply bc function over bc elements
        procedure   :: get_ncoupled_elems       !< Return the number of elements coupled with a specified boundary element.


    end type bc_t
    !*********************************************************************************************





contains

    !> Initialize boundary condition routine
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!  @param[in]  mesh    mesh_t object containing elements and faces
    !!  @param[in]  iface   block face index to which the boundary condition is being applied
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_bc(self,mesh,bconnectivity)
        class(bc_t),                    intent(inout)   :: self
        type(mesh_t),                   intent(inout)   :: mesh
        type(boundary_connectivity_t),  intent(in)      :: bconnectivity

        type(point_t)                   :: pnt, point_one, point_two, point_three
        character(len=:),   allocatable :: bcname        
        real(rk)                        :: time, x, y, z
        integer(ik)                     :: nelem_xi, nelem_eta, nelem_zeta, nelem_bc, ielem_bc, & 
                                           xi_begin, eta_begin, zeta_begin, xi_end, eta_end, zeta_end, & 
                                           ixi, ieta, izeta, ierr, ielem, ielem_test, nface_nodes, iface, inode, i, nfaces_bc, iface_bc

        logical,        allocatable :: node_matched(:), xi_face, eta_face, zeta_face
        integer(ik),    allocatable :: element_nodes(:)
        integer(ik)                 :: face_node


        !
        ! Get number of elements/faces associated with boundary condition.
        !
        nelem_bc = bconnectivity%get_nfaces()




!        !
!        ! Detect number of bcfaces included in the mesh. Could only be a partial match because of parallel partitioning. Sets nfaces_bc
!        !
!        nfaces_bc = 0
!        do ielem_bc = 1,nelem_bc
!
!            ! Allocate array to register if bc node is included in element node list
!            nface_nodes = size(bconnectivity%data(ielem_bc)%data)
!            if ( allocated(node_matched) ) deallocate(node_matched)
!            allocate(node_matched(nface_nodes), stat=ierr)
!            if (ierr /= 0) call AllocationError
!
!
!            ! Search for element with common nodes
!            do ielem = 1,mesh%nelem
!
!                ! Loop through bc nodes and see if current element contains them
!                node_matched = .false.
!                do inode = 1,nface_nodes
!                    element_nodes = mesh%elems(ielem)%connectivity%get_element_nodes()
!                    face_node     = bconnectivity%data(ielem_bc)%get_face_node(inode)
!                    if ( any(element_nodes == face_node) ) then
!                        node_matched(inode) = .true.
!                    end if
!                end do
!
!                ! If all match, increment number of boundary faces to allocate
!                if ( all(node_matched) ) then
!                    nfaces_bc = nfaces_bc + 1 
!                end if
!
!            end do ! ielem
!
!        end do ! ielem_bc
!
!
!
!
!
!
!
!
!        !
!        ! Allocate storage for element and face indices
!        !
!        allocate(self%dom(nfaces_bc), self%elems(nfaces_bc), self%faces(nfaces_bc), self%coupled_elems(nfaces_bc), stat=ierr)
!        if (ierr /= 0) call AllocationError



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
                    ! Set element index
                    !
!                    self%dom(iface_bc)   = mesh%idomain_l
!                    self%elems(iface_bc) = ielem



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


                    xi_face   = ( (abs(point_one%c1_ - point_two%c1_) < 1000.*RKTOL)  .and. (abs(point_one%c1_ - point_three%c1_) < 1000.*RKTOL) )
                    eta_face  = ( (abs(point_one%c2_ - point_two%c2_) < 1000.*RKTOL)  .and. (abs(point_one%c2_ - point_three%c2_) < 1000.*RKTOL) )
                    zeta_face = ( (abs(point_one%c3_ - point_two%c3_) < 1000.*RKTOL)  .and. (abs(point_one%c3_ - point_three%c3_) < 1000.*RKTOL) )


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


                    ! Set face index
!                    self%faces(iface_bc) = iface
                    call self%bc_patch%add_face(mesh%idomain_l,ielem,iface)






                    ! Set face type - 'ftype'
                    bcname = self%get_name()
                    if ( trim(bcname) == 'periodic' ) then
                        !
                        ! Set to ORPHAN face so it will be recognized as chimera in the detection process.
                        !
                        mesh%faces(ielem,iface)%ftype = ORPHAN

                        !
                        ! Set periodic offset from boundary condition to the face. To be used in detection of gq_donor.
                        !
                        if ( self%bcproperties%compute('type', time, pnt) == ONE ) then
                            mesh%faces(ielem,iface)%periodic_type = 'cartesian'
                        else if ( self%bcproperties%compute('type', time, pnt) == TWO ) then
                            mesh%faces(ielem,iface)%periodic_type = 'cylindrical'
                        end if

                        ! time, pnt do nothing here, but interface for function requires them.
                        mesh%faces(ielem,iface)%chimera_offset_x     = self%bcproperties%compute('offset_x',     time, pnt)
                        mesh%faces(ielem,iface)%chimera_offset_y     = self%bcproperties%compute('offset_y',     time, pnt)
                        mesh%faces(ielem,iface)%chimera_offset_z     = self%bcproperties%compute('offset_z',     time, pnt)
                        mesh%faces(ielem,iface)%chimera_offset_theta = self%bcproperties%compute('offset_theta', time, pnt)

                    else
                        !
                        ! Set face to boundary condition face
                        !
                        mesh%faces(ielem,iface)%ftype = BOUNDARY
                    end if



                    ! End search
                    exit



                end if

            end do ! ielem

        end do !ibc_face





        !
        ! Call user-specialized boundary condition initialization
        !
        call self%init_bc_spec(mesh)


        !
        ! Call user-specialized boundary coupling initialization
        !
        call self%init_bc_coupling(mesh,self%bc_patch)


    end subroutine init_bc
    !**********************************************************************************************
    








    !> Default specialized initialization procedure. This is called from the base bc%init procedure
    !! and can be overwritten by derived types to implement specialized initiailization details.
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

        integer(ik) :: ielem_bc, ielem

        !
        ! Have bc_operators initialize the boundary condition coupling
        !
        do iop = 1,size(self%bc_advective_operators)

            call self%bc_advective_operator(iop)%init_boundary_coupling(mesh,bc_patch)

        end do !iop

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


        ncoupled_elems = self%coupled_elems(ielem)%size()


    end function get_ncoupled_elems
    !**********************************************************************************************
















    !>  Apply boundary condition to the mesh and solution
    !!      - Loops through the associated elements(faces) and calls the specialized bc_t%compute
    !!        procedure for computing the rhs and linearization.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !!  @param[in]      mesh    mesh_t defining elements and faces
    !!  @param[inout]   sdata   solverdata_t containing solution, rhs, and linearization(lin) data
    !!  @param[in]      iblk    Block of the linearization for the current element that is being computed (XI_MIN, XI_MAX, eta.)
    !!  @param[inout]   prop    properties_t object containing equationset properties and material_t objects
    !!
    !---------------------------------------------------------------------------------------------
    subroutine apply(self,mesh,sdata,prop)
        class(bc_t),            intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        class(solverdata_t),    intent(inout)   :: sdata
        class(properties_t),    intent(inout)   :: prop

        integer(ik) :: ielem_bc, idomain_l, ielement_l, ielement_c, iface, idonor, iflux, icoupled_elem, ncoupled_elems

        type(chidg_worker_t)    :: worker


        call worker%init(mesh,sdata)


        !
        ! Loop through associated boundary condition elements and call compute routine for the boundary flux calculation
        !
        do ielem_bc = 1,size(self%elems)
            idomain_l   = self%dom(ielem_bc)
            ielement_l  = self%elems(ielem_bc)   ! Get index of the element being operated on
            iface       = self%faces(ielem_bc)   ! Get face index of element 'ielem' that is being operated on


            worker%face_info%idomain_g  = mesh(idomain_l)%elems(ielement_l)%idomain_g
            worker%face_info%idomain_l  = mesh(idomain_l)%elems(ielement_l)%idomain_l
            worker%face_info%ielement_g = mesh(idomain_l)%elems(ielement_l)%ielement_g
            worker%face_info%ielement_l = mesh(idomain_l)%elems(ielement_l)%ielement_l
            worker%face_info%iface      = iface

            worker%function_info%ifcn     = 0       ! Boundary conditions are not tracked.
            worker%function_info%idepend  = 0       ! Chimera interface not applicable on boundary condition.
            worker%function_info%idiff    = BC_BLK  ! Indicates to storage routine in LHS to store in BC section.


            do iop = 1,size(self%bc_advective_operator)
            
                ! For current element, get number of coupled elements.
                ncoupled_elems = self%get_ncoupled_elems(ielem_bc)


                ! Compute current element function enough times to linearize all the coupled elements.
                ! If no coupling accross the face, the ncoupled_elems=1 for just the local interior element.
                do icoupled_elem = 1,ncoupled_elems

                    !
                    ! Get coupled element to linearize against.
                    !
                    ielement_c = self%coupled_elems(ielem_bc)%at(icoupled_elem)
                    worker%function_info%seed%idomain_g  = mesh(idomain_l)%elems(ielement_c)%idomain_g
                    worker%function_info%seed%idomain_l  = mesh(idomain_l)%elems(ielement_c)%idomain_l
                    worker%function_info%seed%ielement_g = mesh(idomain_l)%elems(ielement_c)%ielement_g
                    worker%function_info%seed%ielement_l = mesh(idomain_l)%elems(ielement_c)%ielement_l
                    worker%function_info%seed%iproc      = IRANK

                    !
                    ! For the current boundary element(face), call specialized compute procedure.
                    !
                    call self%compute(worker,prop)

                end do !ielem_c

            end do !iop


        end do !ielem_bc


    end subroutine apply
    !********************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine add_operator(self,bc_operator)
        class(bc_t),            intent(inout)   :: self
        class(bc_operator_t),   intent(in)      :: bc_operator

        integer(ik) :: iop, ierr
        class(bc_operator_wrapper_t),   allocatable :: temp(:) 


        !
        ! Allocate temp storage for (size+1), copy current operators to temp
        !        
        if (allocated(self%bc_advective_operator)) then

            allocate(temp(size(self%bc_advective_operator) + 1), stat=ierr)
            if (ierr /= 0) call AllocationError

            ! Copy previously added operators to temp
            do iop = 1,size(bc_advective_operator)
                allocate(temp(iop)%op, source=self%bc_advective_operator(iop)%op, stat=ierr)
                if (ierr /= 0) call AllocationError
            end do

        else

            allocate(temp(1), stat=ierr)
            if (ierr /= 0) call AllocationError

        end if


        !
        ! Allocate new operator to end
        !
        allocate(temp(size(temp))%op, source=bc_operator, stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Move temp allocation to bc
        !
        call move_alloc(temp, self%bc_advective_operator)


    end subroutine add_operator
    !**************************************************************************************************






end module type_bc
