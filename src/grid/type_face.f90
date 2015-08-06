module type_face
#include  <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, &
                                      ZETA_MIN, ZETA_MAX, XI_DIR, ETA_DIR, ZETA_DIR, &
                                      SPACEDIM, NFACES
    use type_point,             only: point_t
    use type_element,           only: element_t
    use type_quadrature,        only: quadrature_t
    use type_expansion,         only: expansion_t

    implicit none



    !! Face data type
    !!
    !!  NOTE: could be dangerous to declare static arrays of elements using gfortran because
    !!        the compiler doens't have complete finalization rules implemented. Useing allocatables
    !!        seems to work fine.
    !!
    !!
    !------------------------------
    type, public :: face_t
        integer(ik), pointer         :: neqns    => null()
        integer(ik), pointer         :: nterms_s => null()
        integer(ik)                  :: ftype               !> interior (0) or boundary face (1)
        integer(ik)                  :: iface               !> XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, etc
        integer(ik)                  :: iparent             !> Block-local index of parent element
        integer(ik)                  :: ineighbor           !> Block-local index of neighbor element

        !> Geometry
        !---------------------------------------------------------
        type(point_t),  allocatable  :: quad_pts(:)         !> Cartesian coordinates of quadrature nodes
        type(expansion_t), pointer   :: coords => null()    !> Pointer to element coordinates

        !> Metric terms
        !---------------------------------------------------------
        real(rk),       allocatable  :: jinv(:)                     !> array of inverse element jacobians on the face
        real(rk),       allocatable  :: metric(:,:,:)               !> Face metric terms
        real(rk),       allocatable  :: norm(:,:)                   !> Face normals
        real(rk),       pointer      :: invmass(:,:) => null()      !> Pointer to element inverse mass matrix


        !> Quadrature matrices
        !---------------------------------------------------------
        type(quadrature_t),  pointer :: gq     => null()            !> Pointer to solution quadrature instance
        type(quadrature_t),  pointer :: gqmesh => null()            !> Pointer to mesh quadrature instance



        !> Logical tests
        !---------------------------------------------------------
        logical :: geomInitialized = .false.
        logical :: numInitialized  = .false.
    contains
        procedure           :: init_geom
        procedure           :: init_sol
        procedure, public   :: integrate_flux               !> Integrate face flux
        procedure, public   :: integrate_scalar             !> Integrate face scalar

        procedure           :: compute_quadrature_metrics   !> Compute metric terms at quadrature nodes
        procedure           :: compute_quadrature_normals   !> Compute normals at quadrature nodes
        procedure           :: compute_quadrature_coords    !> Compute cartesian coordinates at quadrature nodes

        final               :: destructor
    end type face_t
    !------------------------------

    private
contains


    !> Face geometry initialization procedure
    !!
    !!  Set integer values for face index, face type, parent element index, neighbor element index and coordinates
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in] iface        Element face integer (XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX)
    !!  @param[in] ftype        Face type (0 - interior face, 1 - boundary face)
    !!  @param[in] elem         Parent element which many face members point to
    !!  @param[in] ineighbor    Index of neighbor element in the block

    !!
    !--------------------------------------------------------------------------
    subroutine init_geom(self,iface,ftype,elem,ineighbor)
        class(face_t),      intent(inout)       :: self
        integer(ik),        intent(in)          :: iface
        integer(ik),        intent(in)          :: ftype
        type(element_t),    intent(in), target  :: elem
        integer(ik),        intent(in)          :: ineighbor

        self%iface      = iface
        self%ftype      = ftype
        self%iparent    = elem%ielem
        self%ineighbor  = ineighbor

        self%coords     => elem%coords

        self%geomInitialized = .true.       !> Confirm face grid initialization
    end subroutine



    !> Face initialization procedure
    !!
    !!  Call procedures to compute metrics, normals, and cartesian face coordinates.
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in] elem     Parent element which many face members point to
    !---------------------------------------------------------------------
    subroutine init_sol(self,elem)
        class(face_t),      intent(inout)       :: self
        type(element_t),    intent(in), target  :: elem

        integer(ik) :: ierr, nnodes

        self%neqns      => elem%neqns
        self%nterms_s   => elem%nterms_s
        self%gq         => elem%gq
        self%gqmesh     => elem%gqmesh
        self%invmass    => elem%invmass


        nnodes = self%gq%face%nnodes

        ! Allocate storage for face data structures
        allocate(self%quad_pts(nnodes),                    &
                 self%jinv(nnodes),                        &
                 self%metric(SPACEDIM,SPACEDIM,nnodes),    &
                 self%norm(nnodes,SPACEDIM), stat=ierr)
        if (ierr /= 0) call AllocationError

        call self%compute_quadrature_metrics()
        call self%compute_quadrature_normals()
        call self%compute_quadrature_coords()

        self%numInitialized  = .true.            !> Confirm face numerics were initialized
    end subroutine



    !> Compute metric terms and cell jacobians at face quadrature nodes
    !!
    !!  @author Nathan A. Wukie
    !------------------------------------------------------------------------
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

        associate (gq_f => self%gqmesh%face)
            dxdxi   = matmul(gq_f%ddxi(  :,:,iface), self%coords%var(1))
            dxdeta  = matmul(gq_f%ddeta( :,:,iface), self%coords%var(1))
            dxdzeta = matmul(gq_f%ddzeta(:,:,iface), self%coords%var(1))

            dydxi   = matmul(gq_f%ddxi(  :,:,iface), self%coords%var(2))
            dydeta  = matmul(gq_f%ddeta( :,:,iface), self%coords%var(2))
            dydzeta = matmul(gq_f%ddzeta(:,:,iface), self%coords%var(2))

            dzdxi   = matmul(gq_f%ddxi(  :,:,iface), self%coords%var(3))
            dzdeta  = matmul(gq_f%ddeta( :,:,iface), self%coords%var(3))
            dzdzeta = matmul(gq_f%ddzeta(:,:,iface), self%coords%var(3))
        end associate

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

        ! compute inverse cell mapping jacobian terms
        invjac = dxdxi*dydeta*dzdzeta - dxdeta*dydxi*dzdzeta - &
                 dxdxi*dydzeta*dzdeta + dxdzeta*dydxi*dzdeta + &
                 dxdeta*dydzeta*dzdxi - dxdzeta*dydeta*dzdxi
        self%jinv = invjac
    end subroutine





    !> Compute normal vector components at face quadrature nodes
    !!
    !!  NOTE: These are not unit normals
    !!
    !!  @author Nathan A. Wukie
    !-------------------------------------------------------------------------
    subroutine compute_quadrature_normals(self)
        class(face_t),  intent(inout)   :: self
        integer(ik)                     :: inode, iface, nnodes

        real(rk)    :: dxdxi(self%gq%face%nnodes), dxdeta(self%gq%face%nnodes), dxdzeta(self%gq%face%nnodes)
        real(rk)    :: dydxi(self%gq%face%nnodes), dydeta(self%gq%face%nnodes), dydzeta(self%gq%face%nnodes)
        real(rk)    :: dzdxi(self%gq%face%nnodes), dzdeta(self%gq%face%nnodes), dzdzeta(self%gq%face%nnodes)
        real(rk)    :: invjac(self%gq%face%nnodes)

        iface = self%iface
        nnodes = self%gq%face%nnodes

        associate (gq_f => self%gqmesh%face)
            dxdxi   = matmul(gq_f%ddxi(:,:,iface),  self%coords%var(1))
            dxdeta  = matmul(gq_f%ddeta(:,:,iface), self%coords%var(1))
            dxdzeta = matmul(gq_f%ddzeta(:,:,iface),self%coords%var(1))

            dydxi   = matmul(gq_f%ddxi(:,:,iface),  self%coords%var(2))
            dydeta  = matmul(gq_f%ddeta(:,:,iface), self%coords%var(2))
            dydzeta = matmul(gq_f%ddzeta(:,:,iface),self%coords%var(2))

            dzdxi   = matmul(gq_f%ddxi(:,:,iface),  self%coords%var(3))
            dzdeta  = matmul(gq_f%ddeta(:,:,iface), self%coords%var(3))
            dzdzeta = matmul(gq_f%ddzeta(:,:,iface),self%coords%var(3))
        end associate

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

        ! Reverse normal vectors for faces XI_MIN,ETA_MIN,ZETA_MIN
        if (self%iface == XI_MIN .or. self%iface == ETA_MIN .or. self%iface == ZETA_MIN) then
            self%norm(:,XI_DIR)   = -self%norm(:,XI_DIR)
            self%norm(:,ETA_DIR)  = -self%norm(:,ETA_DIR)
            self%norm(:,ZETA_DIR) = -self%norm(:,ZETA_DIR)
        end if


    end subroutine





    !> Compute cartesian coordinates at face quadrature nodes
    !!
    !!  @author Nathan A. Wukie
    !----------------------------------------------------------------------
    subroutine compute_quadrature_coords(self)
        class(face_t),  intent(inout)   :: self
        integer(ik)                     :: iface, inode
        real(rk)                        :: x(self%gq%face%nnodes),y(self%gq%face%nnodes),z(self%gq%face%nnodes)

        iface = self%iface
        ! compute cartesian coordinates associated with quadrature points
        associate(gq_f => self%gqmesh%face)
            x = matmul(gq_f%val(:,:,iface),self%coords%var(1))
            y = matmul(gq_f%val(:,:,iface),self%coords%var(2))
            z = matmul(gq_f%val(:,:,iface),self%coords%var(3))
        end associate

        do inode = 1,self%gq%face%nnodes
            call self%quad_pts(inode)%set(x(inode),y(inode),z(inode))
        end do

    end subroutine





    !=============================================================================
    !
    !   Integrate boundary flux.
    !
    !   Multiplies fluxes by test functions and integrates.
    !
    !=============================================================================
    subroutine integrate_flux(self,flux_x,flux_y,flux_z,ivar)
        class(face_t), intent(inout)  :: self
        real(rk),      intent(inout)  :: flux_x(:), flux_y(:), flux_z(:)
        integer(ik),   intent(in)     :: ivar

        real(rk),  dimension(self%nterms_s)  :: integral
        integer(4) :: iface

        iface = self%iface

        ! Multiply by quadrature weights
        flux_x = flux_x * self%gq%face%weights(:,iface)
        flux_y = flux_y * self%gq%face%weights(:,iface)
        flux_z = flux_z * self%gq%face%weights(:,iface)

        ! Multiply by column of test functions, integrate, and add to RHS
        integral = matmul(transpose(self%gq%face%val(:,:,iface)),flux_x)
!        self%rhs%vals(:,ivar) = self%rhs%vals(:,varindex) - integral

        integral = matmul(transpose(self%gq%face%val(:,:,iface)),flux_y)
!        self%rhs%vals(:,ivar) = self%rhs%vals(:,varindex) - integral

        integral = matmul(transpose(self%gq%face%val(:,:,iface)),flux_z)
!        self%rhs%vals(:,ivar) = self%rhs%vals(:,varindex) - integral


        ! ADD SECTION FOR AD LINEARIZATION AND ADDING THAT TO THE SYSTEM JACOBIAN MATRIX

    end subroutine


    !=============================================================================
    !
    !   Integrate scalar flux.
    !
    !
    !=============================================================================
    subroutine integrate_scalar(self,scalar,ivar)
        class(face_t),      intent(inout)  :: self
        real(rk),           intent(inout)  :: scalar(:)
        integer(ik),        intent(in)     :: ivar

        real(rk),       dimension(self%nterms_s)    :: integral
        integer(4)                                  :: iface

        iface = self%iface

        ! Multiply by quadrature weights
        scalar = (scalar)  *  (self%gq%face%weights(:,iface))

        ! Multiply by column of test functions, integrate, and add to RHS
        integral = matmul(transpose(self%gq%face%val(:,:,iface)),scalar)
!        self%rhs%vals(:,ivar) = self%rhs%vals(:,ivar) - integral

    end subroutine
















    !=============================================================================
    !
    !
    !   Compute variable at quadrature nodes.
    !
    !
    !=============================================================================
!    subroutine compute_var(self,varindex,vargq)
!        class(face_t),      intent(in)      :: self
!        integer(kind=ik),   intent(in)      :: varindex
!        real(kind=rk),      intent(inout)   :: vargq(:)
!
!        integer(kind=ik)                    :: iface
!        integer(4)                          :: numrows, numcols
!
!        iface = self%iface
!
!        ! Compute variables at volume GQ nodes
!        vargq = matmul(self%gq%face%val(:,:,iface),self%q%vals(:,varindex))
!
!    end subroutine










    subroutine destructor(self)
        type(face_t), intent(inout) :: self


    end subroutine

end module type_face
