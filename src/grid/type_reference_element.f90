module type_reference_element
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: NFACES, ZERO, ONE, XI_DIR, ETA_DIR, ZETA_DIR
    use mod_gauss_legendre, only: quadrature_nodes, quadrature_weights
    use mod_nodes_uniform,  only: uniform_nodes
    use mod_polynomial,     only: polynomial_val, dpolynomial_val
    use mod_inv,            only: inv
    implicit none




    !>  A reference element data type.
    !!
    !!  Contains references nodes/weights, interpolators, and projectors.
    !!
    !!
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
        real(rk),       allocatable :: nodes_r(:,:)          ! Nodes defining the parametric reference element.
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
        real(rk),       allocatable :: nodes_ie(:,:)        ! nodes_3d(nnodes, 3). xi = nodes(:,1), eta = nodes(:,2), zeta = nodes(:,3)
        real(rk),       allocatable :: iweights_e(:)        ! weights_3d(nnodes)
        real(rk),       allocatable :: nodes_if(:,:,:)      ! nodes_2d(nnodes, 3, nfaces). xi = nodes(:,1,iface), eta = nodes(:,2,iface), zeta = nodes(:,3,iface)
        real(rk),       allocatable :: iweights_f(:)        ! weights_2d(nnodes)


        ! Interpolators
        real(rk),       allocatable :: val_e(:,:)           ! element interpolator, expansion value  to elem nodes
        real(rk),       allocatable :: ddxi_e(:,:)          ! element interpolator, expansion ddxi   to elem nodes
        real(rk),       allocatable :: ddeta_e(:,:)         ! element interpolator, expansion ddeta  to elem nodes
        real(rk),       allocatable :: ddzeta_e(:,:)        ! element interpolator, expansion ddzeta to elem nodes
        real(rk),       allocatable :: val_f(:,:,:)         ! face    interpolator, expansion value  to face nodes
        real(rk),       allocatable :: ddxi_f(:,:,:)        ! face    interpolator, expansion ddxi   to face nodes
        real(rk),       allocatable :: ddeta_f(:,:,:)       ! face    interpolator, expansion ddeta  to face nodes
        real(rk),       allocatable :: ddzeta_f(:,:,:)      ! face    interpolator, expansion ddzeta to face nodes


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
        generic             :: nodes        => nodes_element,   nodes_face
        procedure, private  :: nodes_element
        procedure, private  :: nodes_face

        generic             :: weights      => weights_element, weights_face
        procedure, private  :: weights_element
        procedure, private  :: weights_face

        generic             :: interpolator => interpolate_element, interpolate_face
        procedure, private  :: interpolate_element
        procedure, private  :: interpolate_face
        
        ! Query
        procedure   :: nnodes_r     ! number of nodes in the reference node set
        procedure   :: nnodes_ie    ! Number of interpolation nodes in the element set
        procedure   :: nnodes_if    ! Number of interpolation nodes in the face set

        procedure   :: nterms_r     ! number of terms in the reference node data
        procedure   :: nterms_i     ! number of terms in the interpolator node data

    end type reference_element_t
    !********************************************************************






contains




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
        self%nodes_r      = uniform_nodes(element_type,dim=3)



        !
        ! Compute nodes-to-modes projection matrix
        !
        !
        ! Compute the values of each mapping term at each mesh point
        !
        !nnodes = size(self%nodes_r,1)
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
    !!  nodes_ie
    !!  iweights_e (maybe)
    !!  nodes_if
    !!  iweights_f (maybe)
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

        integer(ik)             :: iface, ierr
        real(rk),   allocatable :: nodes_3d(:,:), nodes_2d(:,:)
        real(rk),   allocatable :: weights_3d(:), weights_2d(:)



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
                nodes_3d = uniform_nodes(level,dim=3)
                nodes_2d = uniform_nodes(level,dim=2)


            case('Quadrature')
                nodes_3d   = quadrature_nodes(  nterms=nterms_rule, level=level,dim=3)
                nodes_2d   = quadrature_nodes(  nterms=nterms_rule, level=level,dim=2)
                weights_3d = quadrature_weights(nterms=nterms_rule, level=level,dim=3)
                weights_2d = quadrature_weights(nterms=nterms_rule, level=level,dim=2)

                self%weights_initialized = .true.
            case default
                call chidg_signal_one(FATAL,"reference_element%init_nodes: requested node set not implemented.", trim(node_set))
        end select


        !
        ! Set nodes (if provided)
        !
        if (allocated(nodes_3d) .and. allocated(nodes_2d)) then
            ! Set 3D element nodes
            self%nodes_ie = nodes_3d

            ! Set 2D face nodes
            if (allocated(self%nodes_if)) deallocate(self%nodes_if)
            allocate(self%nodes_if(size(nodes_2d,1),3,NFACES), stat=ierr)
            if (ierr /= 0) call AllocationError

            iface = 1
            self%nodes_if(:,1,iface) = -ONE
            self%nodes_if(:,2,iface) = nodes_2d(:,1)
            self%nodes_if(:,3,iface) = nodes_2d(:,2)

            iface = 2
            self%nodes_if(:,1,iface) = ONE
            self%nodes_if(:,2,iface) = nodes_2d(:,1)
            self%nodes_if(:,3,iface) = nodes_2d(:,2)

            iface = 3
            self%nodes_if(:,1,iface) = nodes_2d(:,1)
            self%nodes_if(:,2,iface) = -ONE
            self%nodes_if(:,3,iface) = nodes_2d(:,2)

            iface = 4
            self%nodes_if(:,1,iface) = nodes_2d(:,1)
            self%nodes_if(:,2,iface) = ONE
            self%nodes_if(:,3,iface) = nodes_2d(:,2)

            iface = 5
            self%nodes_if(:,1,iface) = nodes_2d(:,1)
            self%nodes_if(:,2,iface) = nodes_2d(:,2)
            self%nodes_if(:,3,iface) = -ONE

            iface = 6
            self%nodes_if(:,1,iface) = nodes_2d(:,1)
            self%nodes_if(:,2,iface) = nodes_2d(:,2)
            self%nodes_if(:,3,iface) = ONE


            self%nodes_initialized = .true.


        end if !nodes



        !
        ! Set weights (if provided)
        !
        if (allocated(self%iweights_e)) deallocate(self%iweights_e)
        if (allocated(self%iweights_f)) deallocate(self%iweights_f)
        if (allocated(weights_3d) .and. allocated(weights_2d)) then
            ! Set 3D element weights
            self%iweights_e = weights_3d

            ! Set 2D element weights
            self%iweights_f = weights_2d

            self%weights_initialized = .true.
        end if


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

        integer(ik)             :: iterm, inode, iface, spacedim, ierr, nnodes
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
        if (allocated(self%val_e))    deallocate(self%val_e)
        if (allocated(self%ddxi_e))   deallocate(self%ddxi_e)
        if (allocated(self%ddeta_e))  deallocate(self%ddeta_e)
        if (allocated(self%ddzeta_e)) deallocate(self%ddzeta_e)
        allocate(self%val_e(   self%nnodes_ie(), nterms), &
                 self%ddxi_e(  self%nnodes_ie(), nterms), &
                 self%ddeta_e( self%nnodes_ie(), nterms), &
                 self%ddzeta_e(self%nnodes_ie(), nterms), stat=ierr)
        if (ierr /= 0) call AllocationError

        do iterm = 1,nterms
            do inode = 1,self%nnodes_ie()
                xi   = self%nodes_ie(inode,1)
                eta  = self%nodes_ie(inode,2)
                zeta = self%nodes_ie(inode,3)
                self%val_e(   inode,iterm) =  polynomial_val(3,nterms,iterm,[xi,eta,zeta])
                self%ddxi_e(  inode,iterm) = dpolynomial_val(3,nterms,iterm,[xi,eta,zeta],XI_DIR)
                self%ddeta_e( inode,iterm) = dpolynomial_val(3,nterms,iterm,[xi,eta,zeta],ETA_DIR)
                self%ddzeta_e(inode,iterm) = dpolynomial_val(3,nterms,iterm,[xi,eta,zeta],ZETA_DIR)
            end do
        end do



        !
        ! Compute face interpolators
        !
        if (allocated(self%val_f))    deallocate(self%val_f)
        if (allocated(self%ddxi_f))   deallocate(self%ddxi_f)
        if (allocated(self%ddeta_f))  deallocate(self%ddeta_f)
        if (allocated(self%ddzeta_f)) deallocate(self%ddzeta_f)
        allocate(self%val_f(   self%nnodes_if(), nterms, NFACES), &
                 self%ddxi_f(  self%nnodes_if(), nterms, NFACES), &
                 self%ddeta_f( self%nnodes_if(), nterms, NFACES), &
                 self%ddzeta_f(self%nnodes_if(), nterms, NFACES), stat=ierr)
        if (ierr /= 0) call AllocationError

        do iface = 1,NFACES
            do iterm = 1,nterms
                do inode = 1,self%nnodes_if()
                    xi   = self%nodes_if(inode,1,iface)
                    eta  = self%nodes_if(inode,2,iface)
                    zeta = self%nodes_if(inode,3,iface)
                    self%val_f(   inode,iterm,iface) =  polynomial_val(spacedim,nterms,iterm,[xi,eta,zeta])
                    self%ddxi_f(  inode,iterm,iface) = dpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],XI_DIR)
                    self%ddeta_f( inode,iterm,iface) = dpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ETA_DIR)
                    self%ddzeta_f(inode,iterm,iface) = dpolynomial_val(spacedim,nterms,iterm,[xi,eta,zeta],ZETA_DIR)
                end do
            end do
        end do


    end subroutine init_interpolator_matrices
    !***********************************************************************












    !>
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

        nodes_ = self%nodes_ie

    end function nodes_element
    !*******************************************************************


    !>
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

        nodes_ = self%nodes_if(:,:,iface)

    end function nodes_face
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

        weights_ = self%iweights_e

    end function weights_element
    !*******************************************************************




    !>
    !!
    !!  NOTE: we don't actually use the 'iface' argument in the call.
    !!  it is just there to be consistent with the other calls and it
    !!  also helps differentiate this routine from the weights_element
    !!  routine for the generic interface.
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

        weights_ = self%iweights_f

    end function weights_face
    !*******************************************************************








    !>  Return an interpolator for the element volume node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/28/2017
    !!
    !-------------------------------------------------------------------
    function interpolate_element(self,selector) result(interpolator_)  
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
            case default
                call chidg_signal_one(FATAL,"reference_element%interpolate: Invalid selector for element interpolator.", trim(selector))
        end select
        
    end function interpolate_element
    !*******************************************************************



    !>  Return an interpolator for the face node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/28/2017
    !!
    !-------------------------------------------------------------------
    function interpolate_face(self,selector,iface) result(interpolator_)  
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
            case default
                call chidg_signal_one(FATAL,"reference_element%interpolate: Invalid selector for face interpolator.", trim(selector))
        end select
        
    end function interpolate_face
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
    pure function nnodes_ie(self) result(nnodes)
        class(reference_element_t), intent(in)  :: self

        integer(ik) :: nnodes

        !if (.not. self%interpolation_initialized) call chidg_signal(FATAL,"reference_element%nnodes_i: interpolation has not been initialized.")

        !nnodes = size(self%val_e,1)
        nnodes = size(self%nodes_ie,1)

    end function nnodes_ie
    !*******************************************************************


    !>  Return number of nodes in the face interpolation node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2017
    !!
    !-------------------------------------------------------------------
    pure function nnodes_if(self) result(nnodes)
        class(reference_element_t), intent(in)  :: self

        integer(ik) :: nnodes

        !if (.not. self%interpolation_initialized) call chidg_signal(FATAL,"reference_element%nnodes_i: interpolation has not been initialized.")

        !nnodes = size(self%val_f,1)
        nnodes = size(self%nodes_if,1)

    end function nnodes_if
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


end module type_reference_element
