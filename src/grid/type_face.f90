module type_face
#include  <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX,                 &
                                      ZETA_MIN, ZETA_MAX, XI_DIR, ETA_DIR, ZETA_DIR,    &
                                      X_DIR, Y_DIR, Z_DIR,                              &
                                      TWO_DIM, THREE_DIM,                               &
                                      NFACES, NO_INTERIOR_NEIGHBOR, NO_PROC,            &
                                      ZERO, ONE, TWO, ORPHAN

    use type_point,             only: point_t
    use type_element,           only: element_t
    use type_quadrature,        only: quadrature_t
    use type_densevector,       only: densevector_t
    implicit none



    !>  Face data type
    !!
    !!  ************************************************************************************
    !!  NOTE: could be dangerous to declare static arrays of elements using gfortran because
    !!        the compiler doens't have complete finalization rules implemented. Using 
    !!        allocatables seems to work fine.
    !!  ************************************************************************************
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/23/2016
    !!
    !------------------------------------------------------------------------------------------
    type, public :: face_t
        integer(ik)                 :: spacedim             !< Number of spatial dimensions

        ! Self information
        integer(ik)                 :: ftype                !< INTERIOR, BOUNDARY, CHIMERA, ORPHAN 
        integer(ik)                 :: iface                !< XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, etc
        integer(ik)                 :: ChiID = 0            !< Identifier for domain-local Chimera interfaces

        integer(ik)                 :: BC_ID  = 0           !< Identifier for Boundary Condition index bcs(BC_ID)
        integer(ik)                 :: BC_face = 0          !< Index in bc_patch bcs(BC_ID)%bc_patch%coupled_elements(BC_face)
        integer(ik)                 :: BC_ndepend = 0       !< Number of coupled element if bc face.

        ! Owner-element information
        integer(ik)                 :: idomain_g            !< Global index of the parent domain
        integer(ik)                 :: idomain_l            !< Processor-local index of the parent domain
        integer(ik)                 :: iparent_g            !< Domain-global index of the parent element
        integer(ik)                 :: iparent_l            !< Processor-local index of the parent element
        integer(ik)                 :: neqns                !< Number of equations in equationset_t
        integer(ik)                 :: nterms_s             !< Number of terms in solution polynomial expansion


        ! Neighbor information
        integer(ik)                 :: ineighbor_proc      = NO_PROC    !< MPI processor rank of the neighboring element
        integer(ik)                 :: ineighbor_domain_g  = 0          !< Global index of the neighboring element's domain
        integer(ik)                 :: ineighbor_domain_l  = 0          !< Processor-local index of the neighboring element's domain
        integer(ik)                 :: ineighbor_element_g = 0          !< Domain-global index of the neighboring element
        integer(ik)                 :: ineighbor_element_l = 0          !< Processor-local index of the neighboring element
        integer(ik)                 :: ineighbor_face      = 0
        integer(ik)                 :: ineighbor_neqns     = 0
        integer(ik)                 :: ineighbor_nterms_s  = 0
        integer(ik)                 :: recv_comm           = 0
        integer(ik)                 :: recv_domain         = 0
        integer(ik)                 :: recv_element        = 0

        real(rk),           allocatable :: neighbor_ddx(:,:)            !< Derivative of basis functions in x-direction at quadrature nodes
        real(rk),           allocatable :: neighbor_ddy(:,:)            !< Derivative of basis functions in y-direction at quadrature nodes
        real(rk),           allocatable :: neighbor_ddz(:,:)            !< Derivative of basis functions in z-direction at quadrature nodes
        real(rk),           allocatable :: neighbor_br2_face(:,:)       !< Matrix for computing/obtaining br2 modes at face nodes
        real(rk),           allocatable :: neighbor_br2_vol(:,:)        !< Matrix for computing/obtaining br2 modes at volume nodes
        real(rk),           allocatable :: neighbor_invmass(:,:)    


        ! Chimera face offset. For periodic boundary condition.
        character(len=:),   allocatable :: periodic_type
        real(rk)                        :: chimera_offset_x = 0._rk
        real(rk)                        :: chimera_offset_y = 0._rk
        real(rk)                        :: chimera_offset_z = 0._rk
        real(rk)                        :: chimera_offset_theta = 0._rk


        ! Geometry
        type(point_t),      allocatable :: quad_pts(:)          !< Cartesian coordinates of quadrature nodes
        type(densevector_t)             :: coords               !< Element coordinates

        ! Metric terms
        real(rk),           allocatable :: jinv(:)              !< array of inverse element jacobians on the face
        real(rk),           allocatable :: metric(:,:,:)        !< Face metric terms
        real(rk),           allocatable :: norm(:,:)            !< Face normals
        real(rk),           allocatable :: unorm(:,:)           !< Unit Face normals in cartesian coordinates


        ! Matrices of cartesian gradients of basis/test functions
        real(rk),           allocatable :: ddx(:,:)             !< Derivative of basis functions in x-direction at quadrature nodes
        real(rk),           allocatable :: ddy(:,:)             !< Derivative of basis functions in y-direction at quadrature nodes
        real(rk),           allocatable :: ddz(:,:)             !< Derivative of basis functions in z-direction at quadrature nodes


        ! Quadrature matrices
        type(quadrature_t),  pointer    :: gq     => null()     !< Pointer to solution quadrature instance
        type(quadrature_t),  pointer    :: gqmesh => null()     !< Pointer to mesh quadrature instance


        ! BR2 matrix
        real(rk),           allocatable :: br2_face(:,:)
        real(rk),           allocatable :: br2_vol(:,:)


        ! Logical tests
        logical :: geomInitialized     = .false.
        logical :: neighborInitialized = .false.
        logical :: numInitialized      = .false.



    contains

        procedure           :: init_geom
        procedure           :: init_neighbor
        procedure           :: init_sol

        procedure           :: compute_quadrature_metrics       !< Compute metric terms at quadrature nodes
        procedure           :: compute_quadrature_normals       !< Compute normals at quadrature nodes
        procedure           :: compute_quadrature_coords        !< Compute cartesian coordinates at quadrature nodes
        procedure           :: compute_gradients_cartesian      !< Compute gradients in cartesian coordinates

        procedure           :: get_neighbor_element_g           !< Return neighbor element index
        procedure           :: get_neighbor_element_l           !< Return neighbor element index
        procedure           :: get_neighbor_face                !< Return neighbor face index

        final               :: destructor

    end type face_t
    !******************************************************************************************

    private





contains





    !> Face geometry initialization procedure
    !!
    !!  Set integer values for face index, face type, parent element index, neighbor element
    !!  index and coordinates.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in] iface        Element face integer (XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX)
    !!  @param[in] elem         Parent element which many face members point to
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_geom(self,iface,elem)
        class(face_t),      intent(inout)       :: self
        integer(ik),        intent(in)          :: iface
        type(element_t),    intent(in)          :: elem

        !
        ! Set indices
        !
        self%iface     = iface
        self%ftype     = ORPHAN
        self%spacedim  = elem%spacedim


        !
        ! Set owner element
        !
        self%idomain_g = elem%idomain_g
        self%idomain_l = elem%idomain_l
        self%iparent_g = elem%ielement_g
        self%iparent_l = elem%ielement_l


        !
        ! No neighbor associated at this point
        !
        self%ineighbor_domain_g  = NO_INTERIOR_NEIGHBOR
        self%ineighbor_domain_l  = NO_INTERIOR_NEIGHBOR
        self%ineighbor_element_g = NO_INTERIOR_NEIGHBOR
        self%ineighbor_element_l = NO_INTERIOR_NEIGHBOR
        self%ineighbor_proc      = NO_PROC
        
        !
        ! Set coordinates
        !
        self%coords = elem%coords


        !
        ! Confirm face grid initialization
        !
        self%geomInitialized = .true.

    end subroutine init_geom
    !******************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_neighbor(self,ftype,ineighbor_domain_g,ineighbor_domain_l,              &
                                        ineighbor_element_g,ineighbor_element_l,            &
                                        ineighbor_face,ineighbor_neqns, ineighbor_nterms_s, &
                                        ineighbor_proc)
        class(face_t),  intent(inout)   :: self
        integer(ik),    intent(in)      :: ftype
        integer(ik),    intent(in)      :: ineighbor_domain_g
        integer(ik),    intent(in)      :: ineighbor_domain_l
        integer(ik),    intent(in)      :: ineighbor_element_g
        integer(ik),    intent(in)      :: ineighbor_element_l
        integer(ik),    intent(in)      :: ineighbor_face
        integer(ik),    intent(in)      :: ineighbor_neqns
        integer(ik),    intent(in)      :: ineighbor_nterms_s
        integer(ik),    intent(in)      :: ineighbor_proc


        self%ftype               = ftype
        self%ineighbor_domain_g  = ineighbor_domain_g
        self%ineighbor_domain_l  = ineighbor_domain_l
        self%ineighbor_element_g = ineighbor_element_g
        self%ineighbor_element_l = ineighbor_element_l
        self%ineighbor_face      = ineighbor_face
        self%ineighbor_neqns     = ineighbor_neqns
        self%ineighbor_nterms_s  = ineighbor_nterms_s
        self%ineighbor_proc      = ineighbor_proc


        self%neighborInitialized = .true.

    end subroutine init_neighbor
    !*******************************************************************************************













    !> Face initialization procedure
    !!
    !!  Call procedures to compute metrics, normals, and cartesian face coordinates.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in] elem     Parent element which many face members point to
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_sol(self,elem)
        class(face_t),      intent(inout)       :: self
        type(element_t),    intent(in), target  :: elem

        integer(ik)             :: ierr, nnodes
        real(rk), allocatable   :: tmp(:,:)

        !
        ! Set indices and associate quadrature instances.
        !
        self%neqns      = elem%neqns
        self%nterms_s   = elem%nterms_s
        self%gq         => elem%gq
        self%gqmesh     => elem%gqmesh

        nnodes = self%gq%face%nnodes


        !
        ! (Re)Allocate storage for face data structures.
        !
        if (allocated(self%jinv)) deallocate(self%jinv, self%quad_pts, self%metric, &
                                             self%norm, self%unorm, self%ddx, self%ddy, self%ddz)
        allocate(self%quad_pts(nnodes),                     &
                 self%jinv(nnodes),                         &
                 self%metric(3,3,nnodes),                   &
                 self%norm(nnodes,3),                       &
                 self%unorm(nnodes,3),                      &
                 self%ddx(nnodes,self%nterms_s),            &
                 self%ddy(nnodes,self%nterms_s),            &
                 self%ddz(nnodes,self%nterms_s), stat=ierr) 
        if (ierr /= 0) call AllocationError



        !
        ! Compute metrics, normals, node coordinates
        !
        call self%compute_quadrature_metrics()
        call self%compute_quadrature_normals()
        call self%compute_quadrature_coords()
        call self%compute_gradients_cartesian()


        !
        ! Compute BR2 matrix
        !
        tmp = matmul(elem%invmass,self%gq%face%val_trans(:,:,self%iface))
        self%br2_face = matmul(self%gq%face%val(:,:,self%iface),tmp)
        self%br2_vol  = matmul(self%gq%vol%val(:,:),tmp)


        !
        ! Confirm face numerics were initialized
        !
        self%numInitialized  = .true.

    end subroutine init_sol
    !******************************************************************************************











    !> Compute metric terms and cell jacobians at face quadrature nodes
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  TODO: Generalize 2D physical coordinates. Currently assumes x-y.
    !!
    !------------------------------------------------------------------------------------------
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
    !******************************************************************************************












    !> Compute normal vector components at face quadrature nodes
    !!
    !!  NOTE: be sure to differentiate between normals self%norm and unit-normals self%unorm
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  TODO: Generalize 2D physical coordinates. Currently assumes x-y.
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_quadrature_normals(self)
        class(face_t),  intent(inout)   :: self
        integer(ik)                     :: inode, iface, nnodes

        real(rk)    :: dxdxi(self%gq%face%nnodes), dxdeta(self%gq%face%nnodes), dxdzeta(self%gq%face%nnodes)
        real(rk)    :: dydxi(self%gq%face%nnodes), dydeta(self%gq%face%nnodes), dydzeta(self%gq%face%nnodes)
        real(rk)    :: dzdxi(self%gq%face%nnodes), dzdeta(self%gq%face%nnodes), dzdzeta(self%gq%face%nnodes)
        real(rk)    :: norm_mag(self%gq%face%nnodes)

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
                call chidg_signal(FATAL,"face%compute_quadrature_normals: Invalid face index in face initialization.")
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
    !******************************************************************************************











    !>  Compute matrices containing cartesian gradients of basis/test function
    !!  at each quadrature node.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_gradients_cartesian(self)
        class(face_t),      intent(inout)   :: self

        integer(ik)                         :: iterm,inode,iface,nnodes

        iface  = self%iface
        nnodes = self%gq%face%nnodes



        do iterm = 1,self%nterms_s
            do inode = 1,nnodes
                self%ddx(inode,iterm) = self%metric(1,1,inode) * self%gq%face%ddxi(inode,iterm,iface)   * (ONE/self%jinv(inode)) + &
                                        self%metric(2,1,inode) * self%gq%face%ddeta(inode,iterm,iface)  * (ONE/self%jinv(inode)) + &
                                        self%metric(3,1,inode) * self%gq%face%ddzeta(inode,iterm,iface) * (ONE/self%jinv(inode))

                self%ddy(inode,iterm) = self%metric(1,2,inode) * self%gq%face%ddxi(inode,iterm,iface)   * (ONE/self%jinv(inode)) + &
                                        self%metric(2,2,inode) * self%gq%face%ddeta(inode,iterm,iface)  * (ONE/self%jinv(inode)) + &
                                        self%metric(3,2,inode) * self%gq%face%ddzeta(inode,iterm,iface) * (ONE/self%jinv(inode))

                self%ddz(inode,iterm) = self%metric(1,3,inode) * self%gq%face%ddxi(inode,iterm,iface)   * (ONE/self%jinv(inode)) + &
                                        self%metric(2,3,inode) * self%gq%face%ddeta(inode,iterm,iface)  * (ONE/self%jinv(inode)) + &
                                        self%metric(3,3,inode) * self%gq%face%ddzeta(inode,iterm,iface) * (ONE/self%jinv(inode))
            end do
        end do

    end subroutine compute_gradients_cartesian
    !*******************************************************************************************













    !> Compute cartesian coordinates at face quadrature nodes
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_quadrature_coords(self)
        class(face_t),  intent(inout)   :: self
        integer(ik)                     :: iface, inode
        real(rk)                        :: x(self%gq%face%nnodes), &
                                           y(self%gq%face%nnodes), &
                                           z(self%gq%face%nnodes)

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
    !******************************************************************************************







    !> Return neighbor element index
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    function get_neighbor_element_l(self) result(neighbor_e)
        class(face_t),  intent(in)   ::  self

        integer(ik) :: neighbor_e


        neighbor_e = self%ineighbor_element_l


    end function get_neighbor_element_l
    !******************************************************************************************






    !> Return neighbor element index
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    function get_neighbor_element_g(self) result(neighbor_e)
        class(face_t),  intent(in)   ::  self

        integer(ik) :: neighbor_e


        neighbor_e = self%ineighbor_element_g


    end function get_neighbor_element_g
    !******************************************************************************************











    !> Return neighbor face index
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    function get_neighbor_face(self) result(neighbor_f)
        class(face_t),  intent(in)   ::  self

        integer(ik) :: neighbor_e
        integer(ik) :: neighbor_f



        neighbor_e = self%get_neighbor_element_l()



        if ( neighbor_e == NO_INTERIOR_NEIGHBOR ) then
            
            neighbor_f = NO_INTERIOR_NEIGHBOR

        else

            !& ASSUMPTION: All elements have same orientation.
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
    !******************************************************************************************




















    subroutine destructor(self)
        type(face_t), intent(inout) :: self


    end subroutine




end module type_face
