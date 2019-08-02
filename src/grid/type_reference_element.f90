module type_reference_element
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: NFACES, NEDGES, ZERO, ONE, XI_DIR, ETA_DIR, ZETA_DIR
    use mod_gauss_legendre, only: quadrature_nodes, quadrature_weights
    use mod_nodes_uniform,  only: uniform_nodes, uniform_weights
    use mod_polynomial,     only: polynomial_val, dpolynomial_val, ddpolynomial_val
    use mod_inv,            only: inv
    use mod_gridspace,      only: linspace
    use ieee_arithmetic,    only: ieee_is_nan
    implicit none




    !>  A reference element data type.
    !!
    !!  Contains references nodes/weights, interpolators, and projectors.
    !!  
    !!
    !! nodes_r(.'s) for P1 element     nodes_r(.'s) for P2 element
    !!        .-----------.                  .-----.-----.
    !!       /           /|                 /     /     /|
    !!      /           / |                .-----.-----. |
    !!     /           /  |               /     /     /| .
    !!    .-----------.   |              .-----.-----. |/|
    !!    |           |   |              |     |     | . |
    !!    |           |   .              |     |     |/| .
    !!    |           |  /               .-----.-----. |/
    !!    |           | /                |     |     | .
    !!    |           |/                 |     |     |/
    !!    .-----------.                  .-----.-----.
    !!
    !!
    !!  
    !!  Example of quadrature interpolation nodes(x's) on faces. Note,
    !!  the interpolation nodes are independent of the element reference 
    !!  nodes, self%nodes_r.
    !!
    !!        .-----------.                  .-----.-----.
    !!       /  x    x   /|                 /  x  /  x  /|
    !!      /           / |                .-----.-----. |
    !!     /  x    x   / x|               /  x  /  x  /|x.
    !!    .-----------.   |              .-----.-----. |/|
    !!    |           |x x|              |     |     |x.x|
    !!    |  x    x   |   .              |  x  |  x  |/| .
    !!    |           |x /               .-----.-----.x|/
    !!    |  x    x   | /                |  x  |  x  | .
    !!    |           |/                 |     |     |/
    !!    .-----------.                  .-----.-----.
    !!  
    !!
    !!  The interpolator matrices compute the interpolation of polynomial
    !!  quantities at the interpolation nodes(x's).
    !!
    !!      values_{at inodes} = matmul(val_e, modes)
    !!
    !!  The nodes_to_modes projector matrix, acts on values at the reference
    !!  nodes(nodes_r, .'s). When the nodes_to_modes projector matrix,
    !!  acting on values located at nodes_r, computes modal coefficients
    !!  representing the original quantity as a polynomial expansion on 
    !!  the reference element.
    !!
    !!      modes = matmul(nodes_to_modes, values_{at nodes_r})
    !!
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/26/2017
    !!
    !---------------------------------------------------------------------
    type, public :: reference_element_t

        ! Reference element type/nodes
        integer(ik)                 :: element_type
        real(rk),       allocatable :: nodes_r(:,:)         ! Nodes defining the parametric reference element.
        real(rk),       allocatable :: nodes_to_modes(:,:)  ! linear projector. Takes values at 'nodes_r' and computes modal coefficients
                                                            ! NOTE: nterms of the nodes_to_mdoes projector does not correspond to
                                                            !       the number of terms in the interpolator matrices. The number of
                                                            !       terms in the projector is equal to the number of reference
                                                            !       nodes.

        ! Interpolator settings
        character(:),   allocatable :: polynomial           ! polynomial type being interpolated
        character(:),   allocatable :: node_set             ! node set defined for the interpolators
        integer(ik)                 :: level
        integer(ik)                 :: nterms_rule

        ! Interpolation nodes/weights
        real(rk),       allocatable :: nodes_elem_(:,:)     ! nodes_3d(nnodes, 3). xi = nodes(:,1), eta = nodes(:,2), zeta = nodes(:,3)
        real(rk),       allocatable :: nodes_face_(:,:,:)   ! nodes_2d(nnodes, 3, nfaces). xi = nodes(:,1,iface), eta = nodes(:,2,iface), zeta = nodes(:,3,iface)
        real(rk),       allocatable :: nodes_edge_(:,:,:)   ! nodes_1d(nnodes, 3, nedges). xi = nodes(:,1,iedge), eta = nodes(:,2,iedge), zeta = nodes(:,3,iedge)
        real(rk),       allocatable :: weights_elem_(:)    ! weights_3d(nnodes)
        real(rk),       allocatable :: weights_face_(:)    ! weights_2d(nnodes)
        real(rk),       allocatable :: weights_edge_(:)    ! weights_1d(nnodes)


        ! Interpolators
        real(rk),       allocatable :: val_e(:,:)           ! element interpolator, expansion value  to elem nodes
        real(rk),       allocatable :: ddxi_e(:,:)          ! element interpolator, expansion ddxi   to elem nodes
        real(rk),       allocatable :: ddeta_e(:,:)         ! element interpolator, expansion ddeta  to elem nodes
        real(rk),       allocatable :: ddzeta_e(:,:)        ! element interpolator, expansion ddzeta to elem nodes

        real(rk),       allocatable :: dd_dxidxi_e(:,:)     ! element interpolator
        real(rk),       allocatable :: dd_dxideta_e(:,:)    ! element interpolator
        real(rk),       allocatable :: dd_dxidzeta_e(:,:)   ! element interpolator
        real(rk),       allocatable :: dd_detadeta_e(:,:)   ! element interpolator
        real(rk),       allocatable :: dd_detadzeta_e(:,:)  ! element interpolator
        real(rk),       allocatable :: dd_dzetadzeta_e(:,:) ! element interpolator

        real(rk),       allocatable :: val_f(:,:,:)             ! face interpolator, expansion value  to face nodes
        real(rk),       allocatable :: ddxi_f(:,:,:)            ! face interpolator, expansion ddxi   to face nodes
        real(rk),       allocatable :: ddeta_f(:,:,:)           ! face interpolator, expansion ddeta  to face nodes
        real(rk),       allocatable :: ddzeta_f(:,:,:)          ! face interpolator, expansion ddzeta to face nodes

        real(rk),       allocatable :: dd_dxidxi_f(:,:,:)       ! face interpolator, dxidxi
        real(rk),       allocatable :: dd_dxideta_f(:,:,:)      ! face interpolator, dxideta
        real(rk),       allocatable :: dd_dxidzeta_f(:,:,:)     ! face interpolator, dxidzeta
        real(rk),       allocatable :: dd_detadeta_f(:,:,:)     ! face interpolator, detadeta
        real(rk),       allocatable :: dd_detadzeta_f(:,:,:)    ! face interpolator, detadzeta
        real(rk),       allocatable :: dd_dzetadzeta_f(:,:,:)   ! face interpolator, dzetadzeta

        real(rk),       allocatable :: val_edge(:,:,:)          ! edge interpolator, expansion value  to face nodes
        real(rk),       allocatable :: ddxi_edge(:,:,:)         ! edge interpolator, expansion ddxi   to face nodes
        real(rk),       allocatable :: ddeta_edge(:,:,:)        ! edge interpolator, expansion ddeta  to face nodes
        real(rk),       allocatable :: ddzeta_edge(:,:,:)       ! edge interpolator, expansion ddzeta to face nodes

        real(rk),       allocatable :: dd_dxidxi_edge(:,:,:)     ! edge interpolator, dxidxi
        real(rk),       allocatable :: dd_dxideta_edge(:,:,:)    ! edge interpolator, dxideta
        real(rk),       allocatable :: dd_dxidzeta_edge(:,:,:)   ! edge interpolator, dxidzeta
        real(rk),       allocatable :: dd_detadeta_edge(:,:,:)   ! edge interpolator, detadeta
        real(rk),       allocatable :: dd_detadzeta_edge(:,:,:)  ! edge interpolator, detadzeta
        real(rk),       allocatable :: dd_dzetadzeta_edge(:,:,:) ! edge interpolator, dzetadzeta

        logical :: element_initialized       = .false.
        logical :: nodes_initialized         = .false.
        logical :: weights_initialized       = .false.
        logical :: interpolation_initialized = .false.
    contains

        ! Initialize
        procedure           :: init_element
        procedure           :: init_interpolator
        procedure,  private :: init_interpolator_nodes
        procedure,  private :: init_interpolator_matrices

        ! Produce 
        !generic             :: nodes        => nodes_element, nodes_face, nodes_edge
        procedure   :: nodes_element
        procedure   :: nodes_face
        procedure   :: nodes_edge

        !generic             :: weights      => weights_element, weights_face, weights_edge
        procedure   :: weights_element
        procedure   :: weights_face
        procedure   :: weights_edge

        !generic             :: interpolator => interpolate_element, interpolate_face, interpolate_edge
        procedure   :: interpolator_element
        procedure   :: interpolator_face
        procedure   :: interpolator_edge
        
        ! Query
        procedure   :: nnodes_r     ! number of nodes in the reference node set
        procedure   :: nnodes_elem  ! number of interpolation nodes in the element set
        procedure   :: nnodes_face  ! number of interpolation nodes in the face set
        procedure   :: nnodes_edge  ! number of interpolation nodes in the face set

        procedure   :: nterms_r     ! number of terms in the reference node data
        procedure   :: nterms_i     ! number of terms in the interpolator node data

    end type reference_element_t
    !********************************************************************






contains




!    !>  Initialize the reference element.
!    !!
!    !!  @author Nathan A. Wukie (AFRL)
!    !!  @date   6/29/2017
!    !!
!    !-------------------------------------------------------------------
!    subroutine init_element(self,element_type)
!        class(reference_element_t), intent(inout)   :: self
!        integer(ik),                intent(in)      :: element_type
!
!        real(rk),   allocatable :: temp(:,:)
!        integer(ik)             :: nnodes, iterm, inode, ierr
!
!        !
!        ! Set element type, initialize reference nodes.
!        !
!        self%element_type = element_type
!        self%nodes_r      = uniform_nodes(element_type,dim=3)
!
!
!
!        !
!        ! Compute nodes-to-modes projection matrix
!        !
!        !
!        ! Compute the values of each mapping term at each mesh point
!        !
!        nnodes = self%nnodes_r()
!        if (allocated(temp)) deallocate(temp)
!        allocate(temp(nnodes,nnodes), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!        do iterm = 1,nnodes
!            do inode = 1,nnodes
!                temp(inode,iterm) = polynomial_val(3,nnodes,iterm,[self%nodes_r(inode,1),self%nodes_r(inode,2),self%nodes_r(inode,3)])
!            end do
!        end do
!        ! Invert matrix so that it can multiply a vector of
!        ! element points to compute the mode amplitudes of the x,y mappings
!        self%nodes_to_modes = inv(temp)
!
!
!        self%element_initialized = .true.
!
!
!    end subroutine init_element
!    !*******************************************************************


    !>  Initialize the reference element.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/29/2017
    !!
    !-------------------------------------------------------------------
    subroutine init_element(self,element_type)
        class(reference_element_t), intent(inout)   :: self
        integer(ik),                intent(in)      :: element_type

        real(rk),   allocatable :: temp(:,:)
        integer(ik)             :: nnodes, iterm, inode, ierr

        !
        ! Set element type, initialize reference nodes.
        !
        self%element_type = element_type
        self%nodes_r      = get_reference_nodes(element_type,dim=3)


        !
        ! Compute nodes-to-modes projection matrix
        !
        !
        ! Compute the values of each mapping term at each mesh point
        !
        nnodes = self%nnodes_r()
        if (allocated(temp)) deallocate(temp)
        allocate(temp(nnodes,nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError

        do iterm = 1,nnodes
            do inode = 1,nnodes
                temp(inode,iterm) = polynomial_val(3,nnodes,iterm,[self%nodes_r(inode,1),self%nodes_r(inode,2),self%nodes_r(inode,3)])
            end do
        end do


        ! Invert matrix so that it can multiply a vector of
        ! element points to compute the mode amplitudes of the x,y mappings
        self%nodes_to_modes = inv(temp)


        self%element_initialized = .true.


    end subroutine init_element
    !*******************************************************************








    !>  Initialize interpolation nodes.
    !!
    !!  @param[in]  polynomial  String indicating the polynomial basis
    !!  @param[in]  nterms      Number of terms in the polynomial expansion being interpolated from
    !!  @param[in]  node_set    String indicating the node set being interpolated to
    !!  @param[in]  level       Integer indicator for the level of resolution in the interpolation node set
    !!  @param[in]  nterms_rule Number of terms used as a rule to construct the node set resolution
    !!
    !!  In the case where we want an interpolation from the mesh coordinate expansion
    !!  to the quadrature nodes, we could use the number of terms in the mesh
    !!  expansion for 'nterms' and the number of terms in the solution expansion
    !!  for 'nterms_rule'. In this way, the quadrature node set used to interpolate
    !!  the mesh coordinate expansion will coincide with the interpolation set
    !!  for the solution expansion.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/28/2017
    !!
    !-------------------------------------------------------------------
    subroutine init_interpolator(self,polynomial,nterms,node_set,level,nterms_rule)
        class(reference_element_t), intent(inout)   :: self
        character(*),               intent(in)      :: polynomial
        integer(ik),                intent(in)      :: nterms
        character(*),               intent(in)      :: node_set
        integer(ik),                intent(in)      :: level
        integer(ik),                intent(in)      :: nterms_rule

        character(:),   allocatable :: user_msg

        call self%init_interpolator_nodes(node_set,level,nterms_rule)

        call self%init_interpolator_matrices(polynomial,nterms)

        self%interpolation_initialized = .true.

    end subroutine init_interpolator
    !*******************************************************************







    !>  Initialize the interpolation node set.
    !!
    !!  nodes_elem
    !!  nodes_face
    !!  nodes_edge
    !!  weights_elem_ (maybe)
    !!  weights_face_ (maybe)
    !!  weights_edge_ (maybe)
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/29/2017
    !!
    !-----------------------------------------------------------------------
    subroutine init_interpolator_nodes(self,node_set,level,nterms_rule)
        class(reference_element_t), intent(inout)   :: self
        character(*),               intent(in)      :: node_set
        integer(ik),                intent(in)      :: level
        integer(ik),                intent(in)      :: nterms_rule

        integer(ik)             :: iedge, iface, ierr
        real(rk),   allocatable :: nodes_3d(:,:), nodes_2d(:,:), nodes_1d(:,:)
        real(rk),   allocatable :: weights_3d(:), weights_2d(:), weights_1d(:)



        !
        ! Reset node/weight initialization flags
        !
        self%nodes_initialized   = .false.
        self%weights_initialized = .false.

        self%node_set            = trim(node_set)
        self%level               = level
        self%nterms_rule         = nterms_rule


        !
        ! Generate 3D/2D node sets
        !
        select case(trim(node_set))
            case('Uniform')
                nodes_3d   = uniform_nodes(level,dim=3)
                nodes_2d   = uniform_nodes(level,dim=2)
                nodes_1d   = uniform_nodes(level,dim=1)
                weights_3d = uniform_weights(level,dim=3)
                weights_2d = uniform_weights(level,dim=2)
                weights_1d = uniform_weights(level,dim=1)


            case('Quadrature')
                nodes_3d   = quadrature_nodes(  nterms=nterms_rule, level=level,dim=3)
                nodes_2d   = quadrature_nodes(  nterms=nterms_rule, level=level,dim=2)
                nodes_1d   = quadrature_nodes(  nterms=nterms_rule, level=level,dim=1)
                weights_3d = quadrature_weights(nterms=nterms_rule, level=level,dim=3)
                weights_2d = quadrature_weights(nterms=nterms_rule, level=level,dim=2)
                weights_1d = quadrature_weights(nterms=nterms_rule, level=level,dim=1)

                self%weights_initialized = .true.
            case default
                call chidg_signal_one(FATAL,"reference_element%init_nodes: requested node set not implemented.", trim(node_set))
        end select


        !
        ! Set nodes (if provided)
        !
        if ( allocated(nodes_3d) .and. allocated(nodes_2d) .and. allocated(nodes_1d)) then
            ! Set 3D element nodes
            self%nodes_elem_ = nodes_3d

            ! Set 2D face nodes
            if (allocated(self%nodes_face_)) deallocate(self%nodes_face_)
            allocate(self%nodes_face_(size(nodes_2d,1),3,NFACES), stat=ierr)
            if (ierr /= 0) call AllocationError

            iface = 1
            self%nodes_face_(:,1,iface) = -ONE
            self%nodes_face_(:,2,iface) = nodes_2d(:,1)
            self%nodes_face_(:,3,iface) = nodes_2d(:,2)

            iface = 2
            self%nodes_face_(:,1,iface) = ONE
            self%nodes_face_(:,2,iface) = nodes_2d(:,1)
            self%nodes_face_(:,3,iface) = nodes_2d(:,2)

            iface = 3
            self%nodes_face_(:,1,iface) = nodes_2d(:,1)
            self%nodes_face_(:,2,iface) = -ONE
            self%nodes_face_(:,3,iface) = nodes_2d(:,2)

            iface = 4
            self%nodes_face_(:,1,iface) = nodes_2d(:,1)
            self%nodes_face_(:,2,iface) = ONE
            self%nodes_face_(:,3,iface) = nodes_2d(:,2)

            iface = 5
            self%nodes_face_(:,1,iface) = nodes_2d(:,1)
            self%nodes_face_(:,2,iface) = nodes_2d(:,2)
            self%nodes_face_(:,3,iface) = -ONE

            iface = 6
            self%nodes_face_(:,1,iface) = nodes_2d(:,1)
            self%nodes_face_(:,2,iface) = nodes_2d(:,2)
            self%nodes_face_(:,3,iface) = ONE

            ! Set 1D edge nodes
            if (allocated(self%nodes_edge_)) deallocate(self%nodes_edge_)
            allocate(self%nodes_edge_(size(nodes_1d,1),3,NEDGES), stat=ierr)
            if (ierr /= 0) call AllocationError

            iedge = 1
            self%nodes_edge_(:,1,iedge) = -ONE
            self%nodes_edge_(:,2,iedge) = -ONE
            self%nodes_edge_(:,3,iedge) = nodes_1d(:,1)

            iedge = 2
            self%nodes_edge_(:,1,iedge) = -ONE
            self%nodes_edge_(:,2,iedge) =  ONE
            self%nodes_edge_(:,3,iedge) = nodes_1d(:,1)

            iedge = 3
            self%nodes_edge_(:,1,iedge) =  ONE
            self%nodes_edge_(:,2,iedge) = -ONE
            self%nodes_edge_(:,3,iedge) = nodes_1d(:,1)

            iedge = 4
            self%nodes_edge_(:,1,iedge) = ONE
            self%nodes_edge_(:,2,iedge) = ONE
            self%nodes_edge_(:,3,iedge) = nodes_1d(:,1)

            iedge = 5
            self%nodes_edge_(:,1,iedge) = nodes_1d(:,1)
            self%nodes_edge_(:,2,iedge) = -ONE
            self%nodes_edge_(:,3,iedge) =  ONE

            iedge = 6
            self%nodes_edge_(:,1,iedge) = nodes_1d(:,1)
            self%nodes_edge_(:,2,iedge) = -ONE
            self%nodes_edge_(:,3,iedge) = -ONE

            iedge = 7
            self%nodes_edge_(:,1,iedge) = nodes_1d(:,1)
            self%nodes_edge_(:,2,iedge) = ONE
            self%nodes_edge_(:,3,iedge) = ONE

            iedge = 8
            self%nodes_edge_(:,1,iedge) = nodes_1d(:,1)
            self%nodes_edge_(:,2,iedge) =  ONE
            self%nodes_edge_(:,3,iedge) = -ONE

            iedge = 9
            self%nodes_edge_(:,1,iedge) = -ONE
            self%nodes_edge_(:,2,iedge) = nodes_1d(:,1)
            self%nodes_edge_(:,3,iedge) = -ONE

            iedge = 10
            self%nodes_edge_(:,1,iedge) =  ONE
            self%nodes_edge_(:,2,iedge) = nodes_1d(:,1)
            self%nodes_edge_(:,3,iedge) = -ONE

            iedge = 11
            self%nodes_edge_(:,1,iedge) = -ONE
            self%nodes_edge_(:,2,iedge) = nodes_1d(:,1)
            self%nodes_edge_(:,3,iedge) =  ONE

            iedge = 12
            self%nodes_edge_(:,1,iedge) = ONE
            self%nodes_edge_(:,2,iedge) = nodes_1d(:,1)
            self%nodes_edge_(:,3,iedge) = ONE


            self%nodes_initialized = .true.


        end if !nodes



        !
        ! Set weights (if provided)
        !
        if (allocated(self%weights_elem_)) deallocate(self%weights_elem_)
        if (allocated(self%weights_face_)) deallocate(self%weights_face_)
        if (allocated(self%weights_edge_)) deallocate(self%weights_edge_)

        if (allocated(weights_3d)) self%weights_elem_ = weights_3d
        if (allocated(weights_2d)) self%weights_face_ = weights_2d
        if (allocated(weights_1d)) self%weights_edge_ = weights_1d
        if (allocated(self%weights_elem_) .and. &
            allocated(self%weights_face_) .and. &
            allocated(self%weights_edge_)) self%weights_initialized = .true.

!        if (allocated(weights_3d) .and. allocated(weights_2d) .and. allocated(weights_1d)) then
!            ! Set 3D element weights
!            self%weights_elem_ = weights_3d
!
!            ! Set 2D element weights
!            self%weights_face_ = weights_2d
!
!            self%weights_initialized = .true.
!        end if


    end subroutine init_interpolator_nodes
    !***********************************************************************










    !>  Initialize polynomial interpolation matrices to compute polynomial
    !!  quantities at the interpolation node locations.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/29/2017
    !!
    !-----------------------------------------------------------------------
    subroutine init_interpolator_matrices(self,polynomial,nterms)
        class(reference_element_t), intent(inout)   :: self
        character(*),               intent(in)      :: polynomial
        integer(ik),                intent(in)      :: nterms

        integer(ik)             :: iterm, inode, iedge, iface, spacedim, ierr, nnodes
        real(rk)                :: xi, eta, zeta
        real(rk),   allocatable :: temp(:,:)


        !
        ! Store polynomial
        !
        self%polynomial = trim(polynomial)
        spacedim = 3


        !
        ! Compute element interpolators
        !
        if (allocated(self%val_e)) deallocate(self%val_e,           &
                                              self%ddxi_e,          &
                                              self%ddeta_e,         &
                                              self%ddzeta_e,        &
                                              self%dd_dxidxi_e,     &
                                              self%dd_dxideta_e,    &
                                              self%dd_dxidzeta_e,   &
                                              self%dd_detadeta_e,   &
                                              self%dd_detadzeta_e,  &
                                              self%dd_dzetadzeta_e)
        allocate(self%val_e(          self%nnodes_elem(), nterms),    &
                 self%ddxi_e(         self%nnodes_elem(), nterms),    &
                 self%ddeta_e(        self%nnodes_elem(), nterms),    &
                 self%ddzeta_e(       self%nnodes_elem(), nterms),    &
                 self%dd_dxidxi_e(    self%nnodes_elem(), nterms),    &
                 self%dd_dxideta_e(   self%nnodes_elem(), nterms),    &
                 self%dd_dxidzeta_e(  self%nnodes_elem(), nterms),    &
                 self%dd_detadeta_e(  self%nnodes_elem(), nterms),    &
                 self%dd_detadzeta_e( self%nnodes_elem(), nterms),    &
                 self%dd_dzetadzeta_e(self%nnodes_elem(), nterms), stat=ierr)
        if (ierr /= 0) call AllocationError

        do iterm = 1,nterms
            do inode = 1,self%nnodes_elem()
                xi   = self%nodes_elem_(inode,1)
                eta  = self%nodes_elem_(inode,2)
                zeta = self%nodes_elem_(inode,3)

                ! Value
                self%val_e(   inode,iterm) =  polynomial_val(spacedim,nterms,iterm,[xi,eta,zeta])

                ! First derivatives
                self%ddxi_e(  inode,iterm) = dpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],XI_DIR)
                self%ddeta_e( inode,iterm) = dpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ETA_DIR)
                self%ddzeta_e(inode,iterm) = dpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ZETA_DIR)

                ! Second/mixed derivatives
                self%dd_dxidxi_e(    inode,iterm) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],XI_DIR,XI_DIR  )
                self%dd_dxideta_e(   inode,iterm) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],XI_DIR,ETA_DIR )
                self%dd_dxidzeta_e(  inode,iterm) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],XI_DIR,ZETA_DIR)

                self%dd_detadeta_e(  inode,iterm) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ETA_DIR,ETA_DIR )
                self%dd_detadzeta_e( inode,iterm) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ETA_DIR,ZETA_DIR)

                self%dd_dzetadzeta_e(inode,iterm) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ZETA_DIR,ZETA_DIR)
            end do
        end do



        !
        ! Compute face interpolators
        !
        if (allocated(self%val_f))    deallocate(self%val_f,            &
                                                 self%ddxi_f,           &
                                                 self%ddeta_f,          &
                                                 self%ddzeta_f,         &
                                                 self%dd_dxidxi_f,      &
                                                 self%dd_dxideta_f,     &
                                                 self%dd_dxidzeta_f,    &
                                                 self%dd_detadeta_f,    &
                                                 self%dd_detadzeta_f,   &
                                                 self%dd_dzetadzeta_f)

        allocate(self%val_f(          self%nnodes_face(), nterms, NFACES), &
                 self%ddxi_f(         self%nnodes_face(), nterms, NFACES), &
                 self%ddeta_f(        self%nnodes_face(), nterms, NFACES), &
                 self%ddzeta_f(       self%nnodes_face(), nterms, NFACES), &
                 self%dd_dxidxi_f(    self%nnodes_face(), nterms, NFACES), &
                 self%dd_dxideta_f(   self%nnodes_face(), nterms, NFACES), &
                 self%dd_dxidzeta_f(  self%nnodes_face(), nterms, NFACES), &
                 self%dd_detadeta_f(  self%nnodes_face(), nterms, NFACES), &
                 self%dd_detadzeta_f( self%nnodes_face(), nterms, NFACES), &
                 self%dd_dzetadzeta_f(self%nnodes_face(), nterms, NFACES), stat=ierr)
        if (ierr /= 0) call AllocationError

        do iface = 1,NFACES
            do iterm = 1,nterms
                do inode = 1,self%nnodes_face()
                    xi   = self%nodes_face_(inode,1,iface)
                    eta  = self%nodes_face_(inode,2,iface)
                    zeta = self%nodes_face_(inode,3,iface)

                    ! Value
                    self%val_f(   inode,iterm,iface) =  polynomial_val(spacedim,nterms,iterm,[xi,eta,zeta])

                    ! First-derivatives
                    self%ddxi_f(  inode,iterm,iface) = dpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],XI_DIR)
                    self%ddeta_f( inode,iterm,iface) = dpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ETA_DIR)
                    self%ddzeta_f(inode,iterm,iface) = dpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ZETA_DIR)

                    ! Second/midxed-derivatives
                    self%dd_dxidxi_f(     inode,iterm,iface) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],XI_DIR,XI_DIR)
                    self%dd_dxideta_f(    inode,iterm,iface) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],XI_DIR,ETA_DIR)
                    self%dd_dxidzeta_f(   inode,iterm,iface) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],XI_DIR,ZETA_DIR)
                    self%dd_detadeta_f(   inode,iterm,iface) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ETA_DIR,ETA_DIR)
                    self%dd_detadzeta_f(  inode,iterm,iface) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ETA_DIR,ZETA_DIR)
                    self%dd_dzetadzeta_f( inode,iterm,iface) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ZETA_DIR,ZETA_DIR)

                end do
            end do
        end do


        !
        ! Compute edge interpolators
        !
        if (allocated(self%val_edge)) deallocate(self%val_edge,             &
                                                 self%ddxi_edge,            &
                                                 self%ddeta_edge,           &
                                                 self%ddzeta_edge,          &
                                                 self%dd_dxidxi_edge,       &
                                                 self%dd_dxideta_edge,      &
                                                 self%dd_dxidzeta_edge,     &
                                                 self%dd_detadeta_edge,     &
                                                 self%dd_detadzeta_edge,    &
                                                 self%dd_dzetadzeta_edge)

        allocate(self%val_edge(          self%nnodes_edge(), nterms, NEDGES), &
                 self%ddxi_edge(         self%nnodes_edge(), nterms, NEDGES), &
                 self%ddeta_edge(        self%nnodes_edge(), nterms, NEDGES), &
                 self%ddzeta_edge(       self%nnodes_edge(), nterms, NEDGES), &
                 self%dd_dxidxi_edge(    self%nnodes_edge(), nterms, NEDGES), &
                 self%dd_dxideta_edge(   self%nnodes_edge(), nterms, NEDGES), &
                 self%dd_dxidzeta_edge(  self%nnodes_edge(), nterms, NEDGES), &
                 self%dd_detadeta_edge(  self%nnodes_edge(), nterms, NEDGES), &
                 self%dd_detadzeta_edge( self%nnodes_edge(), nterms, NEDGES), &
                 self%dd_dzetadzeta_edge(self%nnodes_edge(), nterms, NEDGES), stat=ierr)
        if (ierr /= 0) call AllocationError

        do iedge = 1,NEDGES
            do iterm = 1,nterms
                do inode = 1,self%nnodes_edge()
                    xi   = self%nodes_edge_(inode,1,iedge)
                    eta  = self%nodes_edge_(inode,2,iedge)
                    zeta = self%nodes_edge_(inode,3,iedge)

                    ! Value
                    self%val_edge(   inode,iterm,iedge) =  polynomial_val(spacedim,nterms,iterm,[xi,eta,zeta])

                    ! First-derivatives
                    self%ddxi_edge(  inode,iterm,iedge) = dpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],XI_DIR)
                    self%ddeta_edge( inode,iterm,iedge) = dpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ETA_DIR)
                    self%ddzeta_edge(inode,iterm,iedge) = dpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ZETA_DIR)

                    ! Second/midxed-derivatives
                    self%dd_dxidxi_edge(     inode,iterm,iedge) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],XI_DIR,XI_DIR)
                    self%dd_dxideta_edge(    inode,iterm,iedge) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],XI_DIR,ETA_DIR)
                    self%dd_dxidzeta_edge(   inode,iterm,iedge) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],XI_DIR,ZETA_DIR)
                    self%dd_detadeta_edge(   inode,iterm,iedge) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ETA_DIR,ETA_DIR)
                    self%dd_detadzeta_edge(  inode,iterm,iedge) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ETA_DIR,ZETA_DIR)
                    self%dd_dzetadzeta_edge( inode,iterm,iedge) = ddpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ZETA_DIR,ZETA_DIR)

                end do
            end do
        end do










    end subroutine init_interpolator_matrices
    !***********************************************************************












    !>  Return the interpolation nodes for the element.
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------
    function nodes_element(self) result(nodes_)
        class(reference_element_t), intent(in)  :: self

        real(rk),       allocatable :: nodes_(:,:)
        character(:),   allocatable :: user_msg

        user_msg = "reference_element_t: tried to return element nodes without being &
                    initialized. Make sure init_nodes is getting called."
        if (.not. self%nodes_initialized) call chidg_signal(FATAL,user_msg)

        nodes_ = self%nodes_elem_

    end function nodes_element
    !*******************************************************************


    !>  Return the interpolation nodes for a given face.
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------
    function nodes_face(self,iface) result(nodes_)
        class(reference_element_t), intent(in)  :: self
        integer(ik),                intent(in)  :: iface

        real(rk),       allocatable :: nodes_(:,:)
        character(:),   allocatable :: user_msg

        user_msg = "reference_element_t: tried to return face nodes without being &
                    initialized. Make sure init_nodes is getting called."
        if (.not. self%nodes_initialized) call chidg_signal(FATAL,user_msg)

        nodes_ = self%nodes_face_(:,:,iface)

    end function nodes_face
    !*******************************************************************




    !>  Return the interpolation nodes for a given edge.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/6/2017
    !!
    !-------------------------------------------------------------------
    function nodes_edge(self,iedge) result(nodes_)
        class(reference_element_t), intent(in)  :: self
        integer(ik),                intent(in)  :: iedge

        real(rk),       allocatable :: nodes_(:,:)
        character(:),   allocatable :: user_msg

        user_msg = "reference_element_t: tried to return edge nodes without being &
                    initialized. Make sure init_nodes is getting called."
        if (.not. self%nodes_initialized) call chidg_signal(FATAL,user_msg)

        nodes_ = self%nodes_edge_(:,:,iedge)

    end function nodes_edge
    !*******************************************************************



    !>
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------
    function weights_element(self) result(weights_)
        class(reference_element_t), intent(in)  :: self

        real(rk),       allocatable :: weights_(:)
        character(:),   allocatable :: user_msg

        user_msg = "reference_element_t: tried to return element weights without being &
                    initialized. Make sure init_nodes is getting called."
        if (.not. self%weights_initialized) call chidg_signal(FATAL,user_msg)

        weights_ = self%weights_elem_

    end function weights_element
    !*******************************************************************




    !>
    !!
    !!  NOTE: we don't actually use the 'iface' argument in the call.
    !!  it is just there to be consistent with similar procedures.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!
    !!
    !-------------------------------------------------------------------
    function weights_face(self,iface) result(weights_)
        class(reference_element_t), intent(in)  :: self
        integer(ik),                intent(in)  :: iface

        real(rk),       allocatable :: weights_(:)
        character(:),   allocatable :: user_msg

        user_msg = "reference_element_t: tried to return face weights without being &
                    initialized. Make sure init_nodes is getting called."
        if (.not. self%weights_initialized) call chidg_signal(FATAL,user_msg)

        weights_ = self%weights_face_

    end function weights_face
    !*******************************************************************





    !>
    !!
    !!  NOTE: we don't actually use the 'iedge' argument in the call.
    !!  it is just there to be consistent with similar procedures.
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   12/6/2017
    !!
    !-------------------------------------------------------------------
    function weights_edge(self,iedge) result(weights_)
        class(reference_element_t), intent(in)  :: self
        integer(ik),                intent(in)  :: iedge

        real(rk),       allocatable :: weights_(:)
        character(:),   allocatable :: user_msg

        user_msg = "reference_element_t: tried to return edge weighs without being &
                    initialized. Make sure init_nodes is getting called."
        if (.not. self%weights_initialized) call chidg_signal(FATAL,user_msg)

        weights_ = self%weights_edge_

    end function weights_edge
    !*******************************************************************






    !>  Return an interpolator for the element volume node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/28/2017
    !!
    !-------------------------------------------------------------------
    function interpolator_element(self,selector) result(interpolator_)  
        class(reference_element_t), intent(in)  :: self
        character(*),               intent(in)  :: selector

        real(rk),   allocatable :: interpolator_(:,:)


        select case(trim(selector))
            case('Value')
                interpolator_ = self%val_e
            case('ddxi')
                interpolator_ = self%ddxi_e
            case('ddeta')
                interpolator_ = self%ddeta_e
            case('ddzeta')
                interpolator_ = self%ddzeta_e
            case('dxidxi')
                interpolator_ = self%dd_dxidxi_e
            case('dxideta')
                interpolator_ = self%dd_dxideta_e
            case('dxidzeta')
                interpolator_ = self%dd_dxidzeta_e
            case('detadeta')
                interpolator_ = self%dd_detadeta_e
            case('detadzeta')
                interpolator_ = self%dd_detadzeta_e
            case('dzetadzeta')
                interpolator_ = self%dd_dzetadzeta_e
            case default
                call chidg_signal_one(FATAL,"reference_element%interpolator_element: Invalid selector for element interpolator.", trim(selector))
        end select
        
    end function interpolator_element
    !*******************************************************************



    !>  Return an interpolator for the face node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/28/2017
    !!
    !-------------------------------------------------------------------
    function interpolator_face(self,selector,iface) result(interpolator_)  
        class(reference_element_t), intent(in)  :: self
        character(*),               intent(in)  :: selector
        integer(ik),                intent(in)  :: iface

        real(rk),   allocatable :: interpolator_(:,:)


        select case(trim(selector))
            case('Value')
                interpolator_ = self%val_f(:,:,iface)
            case('ddxi')
                interpolator_ = self%ddxi_f(:,:,iface)
            case('ddeta')
                interpolator_ = self%ddeta_f(:,:,iface)
            case('ddzeta')
                interpolator_ = self%ddzeta_f(:,:,iface)
            case('dxidxi')
                interpolator_ = self%dd_dxidxi_f(:,:,iface)
            case('dxideta')
                interpolator_ = self%dd_dxideta_f(:,:,iface)
            case('dxidzeta')
                interpolator_ = self%dd_dxidzeta_f(:,:,iface)
            case('detadeta')
                interpolator_ = self%dd_detadeta_f(:,:,iface)
            case('detadzeta')
                interpolator_ = self%dd_detadzeta_f(:,:,iface)
            case('dzetadzeta')
                interpolator_ = self%dd_dzetadzeta_f(:,:,iface)
            case default
                call chidg_signal_one(FATAL,"reference_element%interpolator_face: Invalid selector for face interpolator.", trim(selector))
        end select
        
    end function interpolator_face
    !*******************************************************************





    !>  Return an interpolator for the edge node set.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/6/2017
    !!
    !-------------------------------------------------------------------
    function interpolator_edge(self,selector,iedge) result(interpolator_)  
        class(reference_element_t), intent(in)  :: self
        character(*),               intent(in)  :: selector
        integer(ik),                intent(in)  :: iedge

        real(rk),   allocatable :: interpolator_(:,:)


        select case(trim(selector))
            case('Value')
                interpolator_ = self%val_edge(:,:,iedge)
            case('ddxi')
                interpolator_ = self%ddxi_edge(:,:,iedge)
            case('ddeta')
                interpolator_ = self%ddeta_edge(:,:,iedge)
            case('ddzeta')
                interpolator_ = self%ddzeta_edge(:,:,iedge)
            case('dxidxi')
                interpolator_ = self%dd_dxidxi_edge(:,:,iedge)
            case('dxideta')
                interpolator_ = self%dd_dxideta_edge(:,:,iedge)
            case('dxidzeta')
                interpolator_ = self%dd_dxidzeta_edge(:,:,iedge)
            case('detadeta')
                interpolator_ = self%dd_detadeta_edge(:,:,iedge)
            case('detadzeta')
                interpolator_ = self%dd_detadzeta_edge(:,:,iedge)
            case('dzetadzeta')
                interpolator_ = self%dd_dzetadzeta_edge(:,:,iedge)
            case default
                call chidg_signal_one(FATAL,"reference_edge%interpolator_edge: Invalid selector for edge interpolator.", trim(selector))
        end select
        
    end function interpolator_edge
    !*******************************************************************






    !>  Return number of nodes in the reference node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2017
    !!
    !-------------------------------------------------------------------
    pure function nnodes_r(self) result(nnodes)
        class(reference_element_t), intent(in)  :: self

        integer(ik) :: nnodes

        !if (.not. self%element_initialized) call chidg_signal(FATAL,"reference_element%nnodes_r: reference element has not been initialized.")

        nnodes = size(self%nodes_r,1)

    end function nnodes_r
    !*******************************************************************




    !>  Return number of nodes in the element interpolation node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2017
    !!
    !-------------------------------------------------------------------
    pure function nnodes_elem(self) result(nnodes)
        class(reference_element_t), intent(in)  :: self

        integer(ik) :: nnodes

        nnodes = size(self%nodes_elem_,1)

    end function nnodes_elem
    !*******************************************************************


    !>  Return number of nodes in the face interpolation node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2017
    !!
    !-------------------------------------------------------------------
    pure function nnodes_face(self) result(nnodes)
        class(reference_element_t), intent(in)  :: self

        integer(ik) :: nnodes

        nnodes = size(self%nodes_face_,1)

    end function nnodes_face
    !*******************************************************************


    !>  Return number of nodes in the edge interpolation node set.
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   12/06/2017
    !!
    !-------------------------------------------------------------------
    pure function nnodes_edge(self) result(nnodes)
        class(reference_element_t), intent(in)  :: self

        integer(ik) :: nnodes

        nnodes = size(self%nodes_edge_,1)

    end function nnodes_edge
    !*******************************************************************





    !>  Return number of terms in the linear projector.
    !!
    !!  NOTE: the number of terms in the linear projector matrix should
    !!        equal the number of nodes defined for the reference element.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2017
    !!
    !-------------------------------------------------------------------
    function nterms_r(self) result(nterms)
        class(reference_element_t), intent(in)  :: self

        integer(ik) :: nterms

        if (.not. self%element_initialized) call chidg_signal(FATAL,"reference_element%nterms_r: reference element has not been initialized.")

        nterms = size(self%nodes_to_modes,1)

    end function nterms_r
    !*******************************************************************




    !>  Return number of terms the interpolators are contstructed for.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2017
    !!
    !-------------------------------------------------------------------
    function nterms_i(self) result(nterms)
        class(reference_element_t), intent(in)  :: self

        integer(ik) :: nterms

        if (.not. self%interpolation_initialized) call chidg_signal(FATAL,"reference_element%nterms_i: reference element interpolators have not been initialized.")

        nterms = size(self%val_e,2)

    end function nterms_i
    !*******************************************************************
















    !>
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/28/2017
    !!
    !-------------------------------------------------------------------------------
    function get_reference_nodes(level, dim) result(nodes)
        integer(ik),    intent(in)  :: level
        integer(ik),    intent(in)  :: dim

        ! Level for quadrature nodes corresponds to the number of nodes
        ! in the 1D tensor construction of the node set.
        integer(ik)                 :: inode, inode_inner, nnodes, ierr, &
                                       nnodes1d, nedge_vertices, nedge_interiors, nface_interiors, &
                                       nelement_interiors, node_start, node_end, line_begin, line_end
        real(rk),       allocatable :: nodes(:,:), line_space(:)
        integer(ik),    allocatable :: nodes_1daccess(:,:), interior_forward(:), interior_reverse(:)
        logical                     :: duplicate_node


        if (dim /= 3) call chidg_signal_one(FATAL,'uniform_nodes: only dim=3 is valid.',dim)


        nnodes1d = level+1
        line_space = linspace(-ONE,ONE,nnodes1d)

        !
        ! Uniform distribution starts with two nodes at end points.
        ! Therefore, minimum level(1) is two points(level+1)
        !
        nnodes1d = level + 1
        nnodes = nnodes1d*nnodes1d*nnodes1d
        allocate(nodes(nnodes,3), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Get tensor product access indices
        nodes_1daccess = get_reference_element_1daccess(level)

        ! Check dimensions are correct
        if (size(nodes_1daccess,1) /= size(nodes,1)) call chidg_signal(FATAL,'uniform_weights: inconsistent node dimensions.')

        ! Assemble nodes
        do inode = 1,size(nodes,1)
            nodes(inode,1) = line_space(nodes_1daccess(inode,1))
            nodes(inode,2) = line_space(nodes_1daccess(inode,2))
            nodes(inode,3) = line_space(nodes_1daccess(inode,3))
        end do !inode


        ! Check to make sure there are no duplicate nodes
        do inode = 1,size(nodes,1)

            duplicate_node = .false.
            do inode_inner = 1,size(nodes,1)
                if ( (inode /= inode_inner) .and. &
                     (nodes(inode,1) == nodes(inode_inner,1)) .and. &
                     (nodes(inode,2) == nodes(inode_inner,2)) .and. &
                     (nodes(inode,3) == nodes(inode_inner,3)) &
                    ) then
                    duplicate_node = .true. 
                end if
            end do !inode_inner

            if (duplicate_node) call chidg_signal(FATAL,'uniform_nodes: found duplicate node. Must be error in node set convention.')

        end do !inode




    end function get_reference_nodes
    !***********************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   8/1/2019
    !!
    !----------------------------------------------------------------------------------
    function get_reference_element_1daccess(order) result(nodes_1daccess)
        integer(ik),    intent(in)  :: order

        integer(ik)                 :: nedge_vertices, nedge_interiors, nface_interiors, nelement_interiors, &
                                       line_begin, line_end, node_start, node_end, nnodes, nnodes1d, ierr, inode, inode_inner
        integer(ik),    allocatable :: nodes_1daccess(:,:), interior_forward(:), interior_reverse(:)
        logical :: duplicate_node


        nnodes1d = order+1
        nnodes = nnodes1d*nnodes1d*nnodes1d
        allocate(nodes_1daccess(nnodes,3), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Fill edge, face, interior nodes
        !
        select case (order)
            ! Linear (HEXA_8)
            case(1)
                nedge_vertices     = 2
                nedge_interiors    = 0
                nface_interiors    = 0
                nelement_interiors = 0
                line_begin = 1
                line_end  = nedge_interiors+nedge_vertices
                !interior_forward = []
                !interior_reverse = []


                ! Hex vertices are same for all orders
                nodes_1daccess(1,1) = line_begin
                nodes_1daccess(1,2) = line_begin
                nodes_1daccess(1,3) = line_begin

                nodes_1daccess(2,1) = line_end
                nodes_1daccess(2,2) = line_begin
                nodes_1daccess(2,3) = line_begin

                nodes_1daccess(3,1) = line_end
                nodes_1daccess(3,2) = line_end
                nodes_1daccess(3,3) = line_begin

                nodes_1daccess(4,1) = line_begin
                nodes_1daccess(4,2) = line_end
                nodes_1daccess(4,3) = line_begin

                nodes_1daccess(5,1) = line_begin
                nodes_1daccess(5,2) = line_begin
                nodes_1daccess(5,3) = line_end

                nodes_1daccess(6,1) = line_end
                nodes_1daccess(6,2) = line_begin
                nodes_1daccess(6,3) = line_end

                nodes_1daccess(7,1) = line_end
                nodes_1daccess(7,2) = line_end
                nodes_1daccess(7,3) = line_end

                nodes_1daccess(8,1) = line_begin
                nodes_1daccess(8,2) = line_end
                nodes_1daccess(8,3) = line_end




            ! Quadratic (HEXA_27)
            case(2)
                nedge_vertices     = 2
                nedge_interiors    = 1
                nface_interiors    = 1
                nelement_interiors = 1
                line_begin = 1
                line_end  = nedge_interiors+nedge_vertices
                interior_forward = [2]
                interior_reverse = [2]




                ! Hex vertices are same for all orders
                nodes_1daccess(1,1) = line_begin
                nodes_1daccess(1,2) = line_begin
                nodes_1daccess(1,3) = line_begin

                nodes_1daccess(2,1) = line_end
                nodes_1daccess(2,2) = line_begin
                nodes_1daccess(2,3) = line_begin

                nodes_1daccess(3,1) = line_end
                nodes_1daccess(3,2) = line_end
                nodes_1daccess(3,3) = line_begin

                nodes_1daccess(4,1) = line_begin
                nodes_1daccess(4,2) = line_end
                nodes_1daccess(4,3) = line_begin

                nodes_1daccess(5,1) = line_begin
                nodes_1daccess(5,2) = line_begin
                nodes_1daccess(5,3) = line_end

                nodes_1daccess(6,1) = line_end
                nodes_1daccess(6,2) = line_begin
                nodes_1daccess(6,3) = line_end

                nodes_1daccess(7,1) = line_end
                nodes_1daccess(7,2) = line_end
                nodes_1daccess(7,3) = line_end

                nodes_1daccess(8,1) = line_begin
                nodes_1daccess(8,2) = line_end
                nodes_1daccess(8,3) = line_end











                ! Edge interiors

                !* Edge 1
                node_start = 9
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = interior_forward
                nodes_1daccess(node_start:node_end, 2) = line_begin
                nodes_1daccess(node_start:node_end, 3) = line_begin

                !* Edge 2
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_end
                nodes_1daccess(node_start:node_end, 2) = interior_forward
                nodes_1daccess(node_start:node_end, 3) = line_begin

                !* Edge 3
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = interior_reverse
                nodes_1daccess(node_start:node_end, 2) = line_end
                nodes_1daccess(node_start:node_end, 3) = line_begin

                !* Edge 4
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_begin
                nodes_1daccess(node_start:node_end, 2) = interior_reverse
                nodes_1daccess(node_start:node_end, 3) = line_begin

                !* Edge 5
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_begin
                nodes_1daccess(node_start:node_end, 2) = line_begin
                nodes_1daccess(node_start:node_end, 3) = interior_forward

                !* Edge 6
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_end
                nodes_1daccess(node_start:node_end, 2) = line_begin
                nodes_1daccess(node_start:node_end, 3) = interior_forward

                !* Edge 7
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_end
                nodes_1daccess(node_start:node_end, 2) = line_end
                nodes_1daccess(node_start:node_end, 3) = interior_forward

                !* Edge 8
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_begin
                nodes_1daccess(node_start:node_end, 2) = line_end
                nodes_1daccess(node_start:node_end, 3) = interior_forward

                !* Edge 9
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = interior_forward
                nodes_1daccess(node_start:node_end, 2) = line_begin
                nodes_1daccess(node_start:node_end, 3) = line_end

                !* Edge 10
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_end
                nodes_1daccess(node_start:node_end, 2) = interior_forward
                nodes_1daccess(node_start:node_end, 3) = line_end

                !* Edge 11
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = interior_reverse
                nodes_1daccess(node_start:node_end, 2) = line_end
                nodes_1daccess(node_start:node_end, 3) = line_end

                !* Edge 12
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_begin
                nodes_1daccess(node_start:node_end, 2) = interior_reverse
                nodes_1daccess(node_start:node_end, 3) = line_end



                ! Face interiors

                ! Face 1
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end,1) = 2
                nodes_1daccess(node_start:node_end,2) = 2
                nodes_1daccess(node_start:node_end,3) = line_begin

                ! Face 2
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end,1) = 2
                nodes_1daccess(node_start:node_end,2) = line_begin
                nodes_1daccess(node_start:node_end,3) = 2

                ! Face 3
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end,1) = line_end
                nodes_1daccess(node_start:node_end,2) = 2
                nodes_1daccess(node_start:node_end,3) = 2

                ! Face 4
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end,1) = 2
                nodes_1daccess(node_start:node_end,2) = line_end
                nodes_1daccess(node_start:node_end,3) = 2

                ! Face 5
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end,1) = line_begin
                nodes_1daccess(node_start:node_end,2) = 2
                nodes_1daccess(node_start:node_end,3) = 2

                ! Face 6
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end,1) = 2
                nodes_1daccess(node_start:node_end,2) = 2
                nodes_1daccess(node_start:node_end,3) = line_end

                !* Interior
                node_start = node_end+1
                node_end   = node_start+nelement_interiors-1
                nodes_1daccess(node_start:node_end,1) = 2
                nodes_1daccess(node_start:node_end,2) = 2
                nodes_1daccess(node_start:node_end,3) = 2



            ! Cubic (HEXA_64)
            case(3)

                nedge_vertices     = 2
                nedge_interiors    = 2
                nface_interiors    = 4
                nelement_interiors = 8
                line_begin = 1
                line_end  = nedge_interiors+nedge_vertices
                interior_forward = [2,3]
                interior_reverse = [3,2]



                ! Hex vertices are same for all orders
                nodes_1daccess(1,1) = line_begin
                nodes_1daccess(1,2) = line_begin
                nodes_1daccess(1,3) = line_begin

                nodes_1daccess(2,1) = line_end
                nodes_1daccess(2,2) = line_begin
                nodes_1daccess(2,3) = line_begin

                nodes_1daccess(3,1) = line_end
                nodes_1daccess(3,2) = line_end
                nodes_1daccess(3,3) = line_begin

                nodes_1daccess(4,1) = line_begin
                nodes_1daccess(4,2) = line_end
                nodes_1daccess(4,3) = line_begin

                nodes_1daccess(5,1) = line_begin
                nodes_1daccess(5,2) = line_begin
                nodes_1daccess(5,3) = line_end

                nodes_1daccess(6,1) = line_end
                nodes_1daccess(6,2) = line_begin
                nodes_1daccess(6,3) = line_end

                nodes_1daccess(7,1) = line_end
                nodes_1daccess(7,2) = line_end
                nodes_1daccess(7,3) = line_end

                nodes_1daccess(8,1) = line_begin
                nodes_1daccess(8,2) = line_end
                nodes_1daccess(8,3) = line_end









                ! Edge interiors
                !* Edge 1
                node_start = 9
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = interior_forward
                nodes_1daccess(node_start:node_end, 2) = line_begin
                nodes_1daccess(node_start:node_end, 3) = line_begin

                !* Edge 2
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_end
                nodes_1daccess(node_start:node_end, 2) = interior_forward
                nodes_1daccess(node_start:node_end, 3) = line_begin

                !* Edge 3
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = interior_reverse
                nodes_1daccess(node_start:node_end, 2) = line_end
                nodes_1daccess(node_start:node_end, 3) = line_begin

                !* Edge 4
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_begin
                nodes_1daccess(node_start:node_end, 2) = interior_reverse
                nodes_1daccess(node_start:node_end, 3) = line_begin

                !* Edge 5
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_begin
                nodes_1daccess(node_start:node_end, 2) = line_begin
                nodes_1daccess(node_start:node_end, 3) = interior_forward

                !* Edge 6
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_end
                nodes_1daccess(node_start:node_end, 2) = line_begin
                nodes_1daccess(node_start:node_end, 3) = interior_forward

                !* Edge 7
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_end
                nodes_1daccess(node_start:node_end, 2) = line_end
                nodes_1daccess(node_start:node_end, 3) = interior_forward

                !* Edge 8
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_begin
                nodes_1daccess(node_start:node_end, 2) = line_end
                nodes_1daccess(node_start:node_end, 3) = interior_forward

                !* Edge 9
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = interior_forward
                nodes_1daccess(node_start:node_end, 2) = line_begin
                nodes_1daccess(node_start:node_end, 3) = line_end

                !* Edge 10
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_end
                nodes_1daccess(node_start:node_end, 2) = interior_forward
                nodes_1daccess(node_start:node_end, 3) = line_end

                !* Edge 11
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = interior_reverse
                nodes_1daccess(node_start:node_end, 2) = line_end
                nodes_1daccess(node_start:node_end, 3) = line_end

                !* Edge 12
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_begin
                nodes_1daccess(node_start:node_end, 2) = interior_reverse
                nodes_1daccess(node_start:node_end, 3) = line_end



                !* Face 1
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end, 1) = [2, 3, 3, 2]
                nodes_1daccess(node_start:node_end, 2) = [2, 2, 3, 3]
                nodes_1daccess(node_start:node_end, 3) = line_begin

                !* Face 2
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end, 1) = [2, 3, 3, 2]
                nodes_1daccess(node_start:node_end, 2) = line_begin
                nodes_1daccess(node_start:node_end, 3) = [2, 2, 3, 3]

                !* Face 3
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_end
                nodes_1daccess(node_start:node_end, 2) = [2, 3, 3, 2]
                nodes_1daccess(node_start:node_end, 3) = [2, 2, 3, 3]

                !* Face 4
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end, 1) = [3, 2, 2, 3]
                nodes_1daccess(node_start:node_end, 2) = line_end
                nodes_1daccess(node_start:node_end, 3) = [2, 2, 3, 3]

                !* Face 5
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_begin
                nodes_1daccess(node_start:node_end, 2) = [3, 2, 2, 3]
                nodes_1daccess(node_start:node_end, 3) = [2, 2, 3, 3]

                !* Face 6
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end, 1) = [2, 3, 3, 2]
                nodes_1daccess(node_start:node_end, 2) = [2, 2, 3, 3]
                nodes_1daccess(node_start:node_end, 3) = line_end

                !* Interior
                node_start = node_end+1
                node_end   = node_start+nelement_interiors-1
                nodes_1daccess(node_start:node_end, 1) = [2, 3, 3, 2, 2, 3, 3, 2]
                nodes_1daccess(node_start:node_end, 2) = [2, 2, 3, 3, 2, 2, 3, 3]
                nodes_1daccess(node_start:node_end, 3) = [2, 2, 2, 2, 3, 3, 3, 3]



            ! Quartic (HEXA_125)
            case(4)

                nedge_vertices     = 2
                nedge_interiors    = 3
                nface_interiors    = 9
                nelement_interiors = 27
                line_begin = 1
                line_end  = nedge_interiors+nedge_vertices
                interior_forward = [2,3,4]
                interior_reverse = [4,3,2]



                ! Hex vertices are same for all orders
                nodes_1daccess(1,1) = line_begin
                nodes_1daccess(1,2) = line_begin
                nodes_1daccess(1,3) = line_begin

                nodes_1daccess(2,1) = line_end
                nodes_1daccess(2,2) = line_begin
                nodes_1daccess(2,3) = line_begin

                nodes_1daccess(3,1) = line_end
                nodes_1daccess(3,2) = line_end
                nodes_1daccess(3,3) = line_begin

                nodes_1daccess(4,1) = line_begin
                nodes_1daccess(4,2) = line_end
                nodes_1daccess(4,3) = line_begin

                nodes_1daccess(5,1) = line_begin
                nodes_1daccess(5,2) = line_begin
                nodes_1daccess(5,3) = line_end

                nodes_1daccess(6,1) = line_end
                nodes_1daccess(6,2) = line_begin
                nodes_1daccess(6,3) = line_end

                nodes_1daccess(7,1) = line_end
                nodes_1daccess(7,2) = line_end
                nodes_1daccess(7,3) = line_end

                nodes_1daccess(8,1) = line_begin
                nodes_1daccess(8,2) = line_end
                nodes_1daccess(8,3) = line_end











                ! Edge interiors
                !* Edge 1
                node_start = 9
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = interior_forward
                nodes_1daccess(node_start:node_end, 2) = line_begin
                nodes_1daccess(node_start:node_end, 3) = line_begin

                !* Edge 2
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_end
                nodes_1daccess(node_start:node_end, 2) = interior_forward
                nodes_1daccess(node_start:node_end, 3) = line_begin

                !* Edge 3
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = interior_reverse
                nodes_1daccess(node_start:node_end, 2) = line_end
                nodes_1daccess(node_start:node_end, 3) = line_begin

                !* Edge 4
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_begin
                nodes_1daccess(node_start:node_end, 2) = interior_reverse
                nodes_1daccess(node_start:node_end, 3) = line_begin

                !* Edge 5
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_begin
                nodes_1daccess(node_start:node_end, 2) = line_begin
                nodes_1daccess(node_start:node_end, 3) = interior_forward

                !* Edge 6
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_end
                nodes_1daccess(node_start:node_end, 2) = line_begin
                nodes_1daccess(node_start:node_end, 3) = interior_forward

                !* Edge 7
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_end
                nodes_1daccess(node_start:node_end, 2) = line_end
                nodes_1daccess(node_start:node_end, 3) = interior_forward

                !* Edge 8
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_begin
                nodes_1daccess(node_start:node_end, 2) = line_end
                nodes_1daccess(node_start:node_end, 3) = interior_forward

                !* Edge 9
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = interior_forward
                nodes_1daccess(node_start:node_end, 2) = line_begin
                nodes_1daccess(node_start:node_end, 3) = line_end

                !* Edge 10
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_end
                nodes_1daccess(node_start:node_end, 2) = interior_forward
                nodes_1daccess(node_start:node_end, 3) = line_end

                !* Edge 11
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = interior_reverse
                nodes_1daccess(node_start:node_end, 2) = line_end
                nodes_1daccess(node_start:node_end, 3) = line_end

                !* Edge 12
                node_start = node_end+1
                node_end   = node_start+nedge_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_begin
                nodes_1daccess(node_start:node_end, 2) = interior_reverse
                nodes_1daccess(node_start:node_end, 3) = line_end



                !* Face 1
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end, 1) = [2, 3, 4, 4, 4, &
                                                          3, 2, 2, 3]
                nodes_1daccess(node_start:node_end, 2) = [2, 2, 2, 3, 4, &
                                                          4, 4, 3, 3]
                nodes_1daccess(node_start:node_end, 3) = line_begin


                !* Face 2
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end, 1) = [2, 3, 4, 4, 4, &
                                                          3, 2, 2, 3]
                nodes_1daccess(node_start:node_end, 2) = line_begin
                nodes_1daccess(node_start:node_end, 3) = [2, 2, 2, 3, 4, &
                                                          4, 4, 3, 3]

                !* Face 3
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_end
                nodes_1daccess(node_start:node_end, 2) = [2, 3, 4, 4, 4, &
                                                          3, 2, 2, 3]
                nodes_1daccess(node_start:node_end, 3) = [2, 2, 2, 3, 4, &
                                                          4, 4, 3, 3]

                !* Face 4
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end, 1) = [4, 3, 2, 2, 2, &
                                                          3, 4, 4, 3]
                nodes_1daccess(node_start:node_end, 2) = line_end
                nodes_1daccess(node_start:node_end, 3) = [2, 2, 2, 3, 4, &
                                                          4, 4, 3, 3]

                !* Face 5
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end, 1) = line_begin
                nodes_1daccess(node_start:node_end, 2) = [4, 3, 2, 2, 2, &
                                                          3, 4, 4, 3]
                nodes_1daccess(node_start:node_end, 3) = [2, 2, 2, 3, 4, &
                                                          4, 4, 3, 3]

                !* Face 6
                node_start = node_end+1
                node_end   = node_start+nface_interiors-1
                nodes_1daccess(node_start:node_end, 1) = [2, 3, 4, 4, 4, &
                                                          3, 2, 2, 3]
                nodes_1daccess(node_start:node_end, 2) = [2, 2, 2, 3, 4, &
                                                          4, 4, 3, 3]
                nodes_1daccess(node_start:node_end, 3) = line_end

                !* Interior
                node_start = node_end+1
                node_end   = node_start+nelement_interiors-1
                nodes_1daccess(node_start:node_end, 1) = [2, 3, 4, 4, 4, &
                                                          3, 2, 2, 3,    &
                                                          2, 3, 4, 4, 4, &
                                                          3, 2, 2, 3,    &
                                                          2, 3, 4, 4, 4, &
                                                          3, 2, 2, 3]
                nodes_1daccess(node_start:node_end, 2) = [2, 2, 2, 3, 4, &
                                                          4, 4, 3, 3,    &
                                                          2, 2, 2, 3, 4, &
                                                          4, 4, 3, 3,    &
                                                          2, 2, 2, 3, 4, &
                                                          4, 4, 3, 3]
                nodes_1daccess(node_start:node_end, 3) = [2, 2, 2, &
                                                          2, 2, 2, & 
                                                          2, 2, 2, & 
                                                          3, 3, 3, & 
                                                          3, 3, 3, & 
                                                          3, 3, 3, & 
                                                          4, 4, 4, & 
                                                          4, 4, 4, & 
                                                          4, 4, 4]




            case default
                call chidg_signal_one(FATAL,"uniform_nodes: invalid value for 'level'", order)

        end select


        ! Check access indices are all unique
        do inode = 1,size(nodes_1daccess,1)
            duplicate_node = .false.
            do inode_inner = 1,size(nodes_1daccess,1)
                if ( (inode /= inode_inner) .and. &
                     (nodes_1daccess(inode,1) == nodes_1daccess(inode_inner,1)) .and. &
                     (nodes_1daccess(inode,2) == nodes_1daccess(inode_inner,2)) .and. &
                     (nodes_1daccess(inode,3) == nodes_1daccess(inode_inner,3)) ) then
                     duplicate_node = .true.
                end if 
            end do
            if (duplicate_node) call chidg_signal(FATAL,'get_reference_element_1daccess: duplicate node detected.')
        end do



    end function get_reference_element_1daccess
    !**********************************************************************************

























































end module type_reference_element
