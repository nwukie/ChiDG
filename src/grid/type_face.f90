module type_face
#include  <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX,                 &
                                      ZETA_MIN, ZETA_MAX, XI_DIR, ETA_DIR, ZETA_DIR,    &
                                      X_DIR, Y_DIR, Z_DIR,                              &
                                      TWO_DIM, THREE_DIM,                               &
                                      NFACES, NO_INTERIOR_NEIGHBOR,                     &
                                      ZERO, ONE, TWO

    use type_point,             only: point_t
    use type_element,           only: element_t
    use type_quadrature,        only: quadrature_t
    use type_densevector,       only: densevector_t
    implicit none



    !> Face data type
    !!
    !!  NOTE: could be dangerous to declare static arrays of elements using gfortran because
    !!        the compiler doens't have complete finalization rules implemented. Useing allocatables
    !!        seems to work fine.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------------------------
    type, public :: face_t
        integer(ik)                  :: spacedim            !< Number of spatial dimensions
        integer(ik)                  :: neqns
        integer(ik)                  :: nterms_s
        integer(ik)                  :: ftype               !< interior, or boundary condition, or Chimera interface, or Orphan
        integer(ik)                  :: iface               !< XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, etc

        integer(ik)                  :: idomain             !< Domain index of parent element
        integer(ik)                  :: iparent             !< Element index of parent element

        integer(ik)                  :: ineighbor           !< Block-local index of neighbor element


        ! Chimera identifier
        integer(ik)                  :: ChiID = 0                   !< Identifier for domain-local Chimera interfaces

        ! Chimera face offset. For periodic boundary condition.
        real(rk)                     :: chimera_offset_x = 0._rk
        real(rk)                     :: chimera_offset_y = 0._rk
        real(rk)                     :: chimera_offset_z = 0._rk


        ! Geometry
        !---------------------------------------------------------
        type(point_t),      allocatable :: quad_pts(:)      !< Cartesian coordinates of quadrature nodes
        type(densevector_t)             :: coords           !< Element coordinates

        ! Metric terms
        !---------------------------------------------------------
        real(rk),           allocatable :: jinv(:)          !< array of inverse element jacobians on the face
        real(rk),           allocatable :: metric(:,:,:)    !< Face metric terms
        real(rk),           allocatable :: norm(:,:)        !< Face normals
        real(rk),           allocatable :: unorm(:,:)       !< Unit Face normals in cartesian coordinates


        ! Matrices of cartesian gradients of basis/test functions
        !---------------------------------------------------------
        real(rk), allocatable       :: dtdx(:,:)            !< Derivative of basis functions in x-direction at quadrature nodes
        real(rk), allocatable       :: dtdy(:,:)            !< Derivative of basis functions in y-direction at quadrature nodes
        real(rk), allocatable       :: dtdz(:,:)            !< Derivative of basis functions in z-direction at quadrature nodes


        ! Quadrature matrices
        !---------------------------------------------------------
        type(quadrature_t),  pointer :: gq     => null()    !< Pointer to solution quadrature instance
        type(quadrature_t),  pointer :: gqmesh => null()    !< Pointer to mesh quadrature instance



        ! Logical tests
        !---------------------------------------------------------
        logical :: geomInitialized = .false.
        logical :: numInitialized  = .false.



    contains

        procedure           :: init_geom
        procedure           :: init_sol

        procedure           :: compute_quadrature_metrics   !< Compute metric terms at quadrature nodes
        procedure           :: compute_quadrature_normals   !< Compute normals at quadrature nodes
        procedure           :: compute_quadrature_coords    !< Compute cartesian coordinates at quadrature nodes
        procedure           :: compute_gradients_cartesian  !< Compute gradients in cartesian coordinates

        procedure           :: get_neighbor_element         !< Return neighbor element index
        procedure           :: get_neighbor_face            !< Return neighbor face index

        final               :: destructor

    end type face_t
    !********************************************************************************************************

    private





contains





    !> Face geometry initialization procedure
    !!
    !!  Set integer values for face index, face type, parent element index, neighbor element index and coordinates
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in] iface        Element face integer (XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX)
    !!  @param[in] ftype        Face type (0 - interior face, 1 - boundary face)
    !!  @param[in] elem         Parent element which many face members point to
    !!  @param[in] ineighbor    Index of neighbor element in the block
    !!
    !----------------------------------------------------------------------------------------------------------
    subroutine init_geom(self,iface,ftype,elem,ineighbor)
        class(face_t),      intent(inout)       :: self
        integer(ik),        intent(in)          :: iface
        integer(ik),        intent(in)          :: ftype
        type(element_t),    intent(in), target  :: elem
        integer(ik),        intent(in)          :: ineighbor

        !
        ! Set indices
        !
        self%iface      = iface
        self%ftype      = ftype
        self%spacedim   = elem%spacedim
        self%idomain    = elem%idomain
        self%iparent    = elem%ielem
        self%ineighbor  = ineighbor

        
        !
        ! Set coordinates
        !
        self%coords     = elem%coords


        !
        ! Confirm face grid initialization
        !
        self%geomInitialized = .true.

    end subroutine init_geom
    !**********************************************************************************************************







    !> Face initialization procedure
    !!
    !!  Call procedures to compute metrics, normals, and cartesian face coordinates.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in] elem     Parent element which many face members point to
    !!
    !----------------------------------------------------------------------------------------------------------
    subroutine init_sol(self,elem)
        class(face_t),      intent(inout)       :: self
        type(element_t),    intent(in), target  :: elem

        integer(ik) :: ierr, nnodes

        !
        ! Set indices and associate quadrature instances.
        !
        self%neqns      = elem%neqns
        self%nterms_s   = elem%nterms_s
        self%gq         => elem%gq
        self%gqmesh     => elem%gqmesh

        nnodes = self%gq%face%nnodes


        !
        ! Allocate storage for face data structures.
        !
        allocate(self%quad_pts(nnodes),                     &
                 self%jinv(nnodes),                         &
                 self%metric(3,3,nnodes),                   &
                 self%norm(nnodes,3),                       &
                 self%unorm(nnodes,3),                      &
                 self%dtdx(nnodes,self%nterms_s),           &
                 self%dtdy(nnodes,self%nterms_s),           &
                 self%dtdz(nnodes,self%nterms_s), stat=ierr) 
        if (ierr /= 0) call AllocationError



        !
        ! Compute metrics, normals, node coordinates
        !
        call self%compute_quadrature_metrics()
        call self%compute_quadrature_normals()
        call self%compute_quadrature_coords()
        call self%compute_gradients_cartesian()


        !
        ! Confirm face numerics were initialized
        !
        self%numInitialized  = .true.

    end subroutine init_sol
    !**********************************************************************************************************











    !> Compute metric terms and cell jacobians at face quadrature nodes
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  TODO: Generalize 2D physical coordinates. Currently assumes x-y.
    !!
    !---------------------------------------------------------------------------------------------------------
    subroutine compute_quadrature_metrics(self)
        class(face_t),  intent(inout)   :: self

        integer(ik) :: inode, iface
        integer(ik) :: nnodes

        real(rk)    :: dxdxi(self%gq%face%nnodes), dxdeta(self%gq%face%nnodes), dxdzeta(self%gq%face%nnodes)
        real(rk)    :: dydxi(self%gq%face%nnodes), dydeta(self%gq%face%nnodes), dydzeta(self%gq%face%nnodes)
        real(rk)    :: dzdxi(self%gq%face%nnodes), dzdeta(self%gq%face%nnodes), dzdzeta(self%gq%face%nnodes)
        real(rk)    :: invjac(self%gq%face%nnodes)

        iface  = self%iface
        nnodes = self%gq%face%nnodes

        !
        ! Evaluate directional derivatives of coordinates at quadrature nodes.
        !
        associate (gq_f => self%gqmesh%face)
            dxdxi   = matmul(gq_f%ddxi(  :,:,iface), self%coords%getvar(1))
            dxdeta  = matmul(gq_f%ddeta( :,:,iface), self%coords%getvar(1))
            dxdzeta = matmul(gq_f%ddzeta(:,:,iface), self%coords%getvar(1))

            dydxi   = matmul(gq_f%ddxi(  :,:,iface), self%coords%getvar(2))
            dydeta  = matmul(gq_f%ddeta( :,:,iface), self%coords%getvar(2))
            dydzeta = matmul(gq_f%ddzeta(:,:,iface), self%coords%getvar(2))

            dzdxi   = matmul(gq_f%ddxi(  :,:,iface), self%coords%getvar(3))
            dzdeta  = matmul(gq_f%ddeta( :,:,iface), self%coords%getvar(3))
            dzdzeta = matmul(gq_f%ddzeta(:,:,iface), self%coords%getvar(3))
        end associate

        

        !
        ! TODO: Generalize 2D physical coordinates. Currently assumes x-y.
        !
        if ( self%spacedim == TWO_DIM ) then
            dzdxi   = ZERO
            dzdeta  = ZERO
            dzdzeta = ONE
        end if



        !
        ! At each quadrature node, compute metric terms.
        !
        do inode = 1,nnodes
            self%metric(1,1,inode) = dydeta(inode)*dzdzeta(inode) - dydzeta(inode)*dzdeta(inode)
            self%metric(2,1,inode) = dydzeta(inode)*dzdxi(inode)  - dydxi(inode)*dzdzeta(inode)
            self%metric(3,1,inode) = dydxi(inode)*dzdeta(inode)   - dydeta(inode)*dzdxi(inode)

            self%metric(1,2,inode) = dxdzeta(inode)*dzdeta(inode) - dxdeta(inode)*dzdzeta(inode)
            self%metric(2,2,inode) = dxdxi(inode)*dzdzeta(inode)  - dxdzeta(inode)*dzdxi(inode)
            self%metric(3,2,inode) = dxdeta(inode)*dzdxi(inode)   - dxdxi(inode)*dzdeta(inode)

            self%metric(1,3,inode) = dxdeta(inode)*dydzeta(inode) - dxdzeta(inode)*dydeta(inode)
            self%metric(2,3,inode) = dxdzeta(inode)*dydxi(inode)  - dxdxi(inode)*dydzeta(inode)
            self%metric(3,3,inode) = dxdxi(inode)*dydeta(inode)   - dxdeta(inode)*dydxi(inode)
        end do

        
        !
        ! compute inverse cell mapping jacobian terms
        !
        invjac = dxdxi*dydeta*dzdzeta - dxdeta*dydxi*dzdzeta - &
                 dxdxi*dydzeta*dzdeta + dxdzeta*dydxi*dzdeta + &
                 dxdeta*dydzeta*dzdxi - dxdzeta*dydeta*dzdxi



        self%jinv = invjac

    end subroutine compute_quadrature_metrics
    !*****************************************************************************************************












    !> Compute normal vector components at face quadrature nodes
    !!
    !!  NOTE: be sure to differentiate between normals self%norm and unit-normals self%unorm
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  TODO: Generalize 2D physical coordinates. Currently assumes x-y.
    !!
    !------------------------------------------------------------------------------------------------------
    subroutine compute_quadrature_normals(self)
        class(face_t),  intent(inout)   :: self
        integer(ik)                     :: inode, iface, nnodes

        real(rk)    :: dxdxi(self%gq%face%nnodes), dxdeta(self%gq%face%nnodes), dxdzeta(self%gq%face%nnodes)
        real(rk)    :: dydxi(self%gq%face%nnodes), dydeta(self%gq%face%nnodes), dydzeta(self%gq%face%nnodes)
        real(rk)    :: dzdxi(self%gq%face%nnodes), dzdeta(self%gq%face%nnodes), dzdzeta(self%gq%face%nnodes)
        real(rk)    :: invjac(self%gq%face%nnodes), norm_mag(self%gq%face%nnodes)

        iface = self%iface
        nnodes = self%gq%face%nnodes

        
        !
        ! Evaluate directional derivatives of coordinates at quadrature nodes.
        !
        associate (gq_f => self%gqmesh%face)
            dxdxi   = matmul(gq_f%ddxi(:,:,iface),  self%coords%getvar(1))
            dxdeta  = matmul(gq_f%ddeta(:,:,iface), self%coords%getvar(1))
            dxdzeta = matmul(gq_f%ddzeta(:,:,iface),self%coords%getvar(1))

            dydxi   = matmul(gq_f%ddxi(:,:,iface),  self%coords%getvar(2))
            dydeta  = matmul(gq_f%ddeta(:,:,iface), self%coords%getvar(2))
            dydzeta = matmul(gq_f%ddzeta(:,:,iface),self%coords%getvar(2))

            dzdxi   = matmul(gq_f%ddxi(:,:,iface),  self%coords%getvar(3))
            dzdeta  = matmul(gq_f%ddeta(:,:,iface), self%coords%getvar(3))
            dzdzeta = matmul(gq_f%ddzeta(:,:,iface),self%coords%getvar(3))
        end associate



        !
        ! TODO: Generalize 2D physical coordinates. Currently assumes x-y.
        !
        if ( self%spacedim == TWO_DIM ) then
            dzdxi   = ZERO
            dzdeta  = ZERO
            dzdzeta = ONE
        end if


        !
        ! Compute normal vectors for each face
        !
        select case (self%iface)
            case (XI_MIN, XI_MAX)

                do inode = 1,nnodes
                    self%norm(inode,XI_DIR)   = dydeta(inode)*dzdzeta(inode) - dydzeta(inode)*dzdeta(inode)
                    self%norm(inode,ETA_DIR)  = dxdzeta(inode)*dzdeta(inode) - dxdeta(inode)*dzdzeta(inode)
                    self%norm(inode,ZETA_DIR) = dxdeta(inode)*dydzeta(inode) - dxdzeta(inode)*dydeta(inode)
                end do

            case (ETA_MIN, ETA_MAX)

                do inode = 1,nnodes
                    self%norm(inode,XI_DIR)   = dydzeta(inode)*dzdxi(inode)  - dydxi(inode)*dzdzeta(inode)
                    self%norm(inode,ETA_DIR)  = dxdxi(inode)*dzdzeta(inode)  - dxdzeta(inode)*dzdxi(inode)
                    self%norm(inode,ZETA_DIR) = dxdzeta(inode)*dydxi(inode)  - dxdxi(inode)*dydzeta(inode)
                end do

            case (ZETA_MIN, ZETA_MAX)

                do inode = 1,nnodes
                    self%norm(inode,XI_DIR)   = dydxi(inode)*dzdeta(inode)   - dzdxi(inode)*dydeta(inode)
                    self%norm(inode,ETA_DIR)  = dzdxi(inode)*dxdeta(inode)   - dxdxi(inode)*dzdeta(inode)
                    self%norm(inode,ZETA_DIR) = dxdxi(inode)*dydeta(inode)   - dydxi(inode)*dxdeta(inode)
                end do

            case default
                stop "Error: invalid face index in face initialization"
        end select


        
        !
        ! Reverse normal vectors for faces XI_MIN,ETA_MIN,ZETA_MIN
        !
        if (self%iface == XI_MIN .or. self%iface == ETA_MIN .or. self%iface == ZETA_MIN) then
            self%norm(:,XI_DIR)   = -self%norm(:,XI_DIR)
            self%norm(:,ETA_DIR)  = -self%norm(:,ETA_DIR)
            self%norm(:,ZETA_DIR) = -self%norm(:,ZETA_DIR)
        end if



        !
        ! Compute unit normals
        !
        norm_mag = sqrt(self%norm(:,XI_DIR)**TWO + self%norm(:,ETA_DIR)**TWO + self%norm(:,ZETA_DIR)**TWO)
        self%unorm(:,XI_DIR)   = self%norm(:,XI_DIR)/norm_mag
        self%unorm(:,ETA_DIR)  = self%norm(:,ETA_DIR)/norm_mag
        self%unorm(:,ZETA_DIR) = self%norm(:,ZETA_DIR)/norm_mag





    end subroutine compute_quadrature_normals
    !************************************************************************************************************











    !>  Compute matrices containing cartesian gradients of basis/test function
    !!  at each quadrature node.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/1/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------------
    subroutine compute_gradients_cartesian(self)
        class(face_t),      intent(inout)   :: self

        integer(ik)                         :: iterm,inode,iface,nnodes

        iface  = self%iface
        nnodes = self%gq%face%nnodes



        do iterm = 1,self%nterms_s
            do inode = 1,nnodes
                self%dtdx(inode,iterm) = self%metric(1,1,inode) * self%gq%face%ddxi(inode,iterm,iface)   * (ONE/self%jinv(inode)) + &
                                         self%metric(2,1,inode) * self%gq%face%ddeta(inode,iterm,iface)  * (ONE/self%jinv(inode)) + &
                                         self%metric(3,1,inode) * self%gq%face%ddzeta(inode,iterm,iface) * (ONE/self%jinv(inode))

                self%dtdy(inode,iterm) = self%metric(1,2,inode) * self%gq%face%ddxi(inode,iterm,iface)   * (ONE/self%jinv(inode)) + &
                                         self%metric(2,2,inode) * self%gq%face%ddeta(inode,iterm,iface)  * (ONE/self%jinv(inode)) + &
                                         self%metric(3,2,inode) * self%gq%face%ddzeta(inode,iterm,iface) * (ONE/self%jinv(inode))

                self%dtdz(inode,iterm) = self%metric(1,3,inode) * self%gq%face%ddxi(inode,iterm,iface)   * (ONE/self%jinv(inode)) + &
                                         self%metric(2,3,inode) * self%gq%face%ddeta(inode,iterm,iface)  * (ONE/self%jinv(inode)) + &
                                         self%metric(3,3,inode) * self%gq%face%ddzeta(inode,iterm,iface) * (ONE/self%jinv(inode))
            end do
        end do

    end subroutine compute_gradients_cartesian
    !*********************************************************************************************************













    !> Compute cartesian coordinates at face quadrature nodes
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine compute_quadrature_coords(self)
        class(face_t),  intent(inout)   :: self
        integer(ik)                     :: iface, inode
        real(rk)                        :: x(self%gq%face%nnodes),y(self%gq%face%nnodes),z(self%gq%face%nnodes)

        iface = self%iface

        !
        ! compute real coordinates associated with quadrature points
        !
        associate(gq_f => self%gqmesh%face)
            x = matmul(gq_f%val(:,:,iface),self%coords%getvar(1))
            y = matmul(gq_f%val(:,:,iface),self%coords%getvar(2))
            z = matmul(gq_f%val(:,:,iface),self%coords%getvar(3))
        end associate

        
        !
        ! For each quadrature node, store real coordinates
        !
        do inode = 1,self%gq%face%nnodes
            call self%quad_pts(inode)%set(x(inode),y(inode),z(inode))
        end do

    end subroutine compute_quadrature_coords
    !***********************************************************************************************************







    !> Return neighbor element index
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------------------
    function get_neighbor_element(self) result(neighbor_e)
        class(face_t),  intent(in)   ::  self

        integer(ik) :: neighbor_e


        neighbor_e = self%ineighbor


    end function get_neighbor_element
    !***********************************************************************************************************





    !> Return neighbor face index
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------------------
    function get_neighbor_face(self) result(neighbor_f)
        class(face_t),  intent(in)   ::  self

        integer(ik) :: neighbor_e
        integer(ik) :: neighbor_f



        neighbor_e = self%get_neighbor_element()



        if ( neighbor_e == NO_INTERIOR_NEIGHBOR ) then
            
            neighbor_f = NO_INTERIOR_NEIGHBOR

        else

            !& ASSUMPTION: BLOCK-STRUCTURED
            if ( self%iface == XI_MIN ) then
                neighbor_f = XI_MAX
            else if ( self%iface == XI_MAX ) then
                neighbor_f = XI_MIN
            else if ( self%iface == ETA_MIN ) then
                neighbor_f = ETA_MAX
            else if ( self%iface == ETA_MAX ) then
                neighbor_f = ETA_MIN
            else if ( self%iface == ZETA_MIN ) then
                neighbor_f = ZETA_MAX
            else if ( self%iface == ZETA_MAX ) then
                neighbor_f = ZETA_MIN
            end if

        end if


    end function get_neighbor_face
    !***********************************************************************************************************




















    subroutine destructor(self)
        type(face_t), intent(inout) :: self


    end subroutine




end module type_face
