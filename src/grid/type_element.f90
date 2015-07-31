module type_element
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: SPACEDIM,NFACES,XI_MIN,XI_MAX,ETA_MIN, &
                                      ETA_MAX,ZETA_MIN,ZETA_MAX,ONE,ZERO
    use type_point,             only: point_t
    use type_expansion,         only: expansion_t
    use type_quadrature,        only: quadrature_t
    use DNAD_D
    use mod_quadrature,         only: GQ, get_quadrature
    use mod_grid,               only: ELEM_MAP
    use mod_polynomial,         only: polynomialVal
    use mod_grid_tools,         only: compute_modal_coordinates
    use mod_inv,                only: inv
    implicit none


    !> Element data type
    !!
    !!  NOTE: could be dangerous to declare static arrays of elements using gfortran because
    !!        the compiler doens't have complete finalization rules implemented. Useing allocatables
    !!        seems to work fine.
    !!
    !!
    !------------------------------------------------------------
    type, public :: element_t
        integer(ik)      :: neqns            !> Number of equations being solved
        integer(ik)      :: nterms_s         !> Number of terms in solution expansion.   (polynomial_order + 1) or (Ns + 1)
        integer(ik)      :: nterms_c         !> Number of terms in coordinate expansion. (polynomial_order + 1) or (Nc + 1)
        integer(ik)      :: ielem            !> Block-local element index

        !> Element quadrature points, mesh points and modes
        !---------------------------------------------------------
        type(point_t), allocatable  :: quad_pts(:)  !> Cartesian coordinates of discrete quadrature points
        type(point_t), allocatable  :: elem_pts(:)  !> Cartesian coordinates of discrete points defining element
        type(expansion_t)           :: coords       !> Modal representation of cartesian coordinates (nterms_var,(x,y,z))

        !> Element metric terms
        !---------------------------------------------------------
        real(rk), allocatable       :: metric(:,:,:)    !> metric matrix for each quadrature node    (mat_i,mat_j,quad_pt)
        real(rk), allocatable       :: jinv(:)          !> jacobian terms at quadrature nodes

        !> Matrices of cartesian gradients of basis/test functions
        !---------------------------------------------------------
        real(rk), allocatable       :: dtdx(:,:)
        real(rk), allocatable       :: dtdy(:,:)
        real(rk), allocatable       :: dtdz(:,:)

        !> Quadrature matrices
        !---------------------------------------------------------
        type(quadrature_t), pointer  :: gq     => null()
        type(quadrature_t), pointer  :: gqmesh => null()

        !> Element-local mass matrices
        !---------------------------------------------------------
        real(rk), allocatable   :: mass(:,:)
        real(rk), allocatable   :: invmass(:,:)

        !> Logical tests
        logical :: geomInitialized = .false.
        logical :: numInitialized  = .false.
    contains
        ! Initialization procedures
        procedure, public   :: init_geom
        procedure, public   :: init_sol

        ! Solution procedures
        procedure, public   :: integrate_volume_flux
        procedure, public   :: integrate_volume_source

        ! Public utility procedures
        procedure, public   :: x
        procedure, public   :: y
        procedure, public   :: z

        ! Private utility procedure
        procedure           :: compute_element_matrices
        procedure           :: compute_mass_matrix
        procedure           :: compute_gradients_cartesian
        procedure           :: compute_quadrature_metrics
        procedure           :: compute_quadrature_coords
        procedure           :: assign_quadrature

        final               :: destructor
    end type element_t
    !-----------------------
    private
contains
    


    !>  Initialize element geometry
    !!      - Set element points
    !!      - Compute modal representation of element cartesian coordinates
    !!      - Compute element metric terms
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in] nterms_c     Number of terms in the modal representation of the cartesian coordinates
    !!  @param[in] points       Array of cartesian points defining the element
    !---------------------------------------------------------------------------------------
    subroutine init_geom(self,mapping,points,ielem)
        class(element_t),   intent(inout) :: self
        integer(ik),        intent(in)    :: mapping
        type(point_t),      intent(in)    :: points(:)
        integer(ik),        intent(in)    :: ielem

        integer(ik)                         :: ierr, nterms_c
        real(rk), dimension(:,:), pointer   :: imap => null()

        if (self%geomInitialized) call signal(FATAL,'element%init_geom -- element already initialized')



        if (allocated(elem_map(mapping)%mat)) then
            imap            => elem_map(mapping)%mat
            nterms_c        = size(elem_map(mapping)%mat,1) !> Get number of terms if coordinate expansion from size of mapping matrix
            self%nterms_c   = nterms_c                      !> Set number of terms in coordinate expansion


            if (nterms_c /= size(points)) call signal(FATAL,'element%init_geom -- mapping and points do not match')
        else
            call signal(FATAL,'element%init_geom -- element mapping not initialized')
        end if

        ! Allocate and compute mesh x,y,z modes
        allocate(self%elem_pts(nterms_c),stat=ierr)
        call self%coords%init(nterms_c,SPACEDIM)
        self%ielem    = ielem
        self%elem_pts = points

        call compute_modal_coordinates(self%elem_pts,mapping,self%coords)

        self%geomInitialized = .true.                       !> Confirm element grid was initialized
    end subroutine






    !>  Initialize element numerics
    !!      - Allocate storage for solution and supporting matrices
    !!      - Compute element-local matrices (cartesian gradients, mass matrices)
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]  nterms_s    Number of terms in the modal representation of the solution
    !!  @param[in]  neqns       Number of equations contained in the element solution
    !--------------------------------------------------------------------------------------
    subroutine init_sol(self,neqns,nterms_s)
        class(element_t),   intent(inout) :: self
        integer(ik),        intent(in)    :: neqns
        integer(ik),        intent(in)    :: nterms_s

        integer(ik) :: ierr
        integer(ik) :: nnodes,nnodes_face,nnodes_vol

        if (self%numInitialized) call signal(FATAL,'element%init_sol -- element already initialized')


        self%nterms_s    = nterms_s                 !> Set number of terms in modal expansion of solution
        self%neqns       = neqns                    !> Set number of equations being solved

        call self%assign_quadrature()               !> With nterms_s and nterms_c defined, we can assign a quadrature instance
        nnodes = self%gq%vol%nnodes                 !> With a quadrature instance assigned, we have the number of quadrature nodes


        allocate(self%jinv(nnodes),stat=ierr)
        allocate(self%metric(SPACEDIM,SPACEDIM,nnodes),stat=ierr)
        if (ierr /= 0) call AllocationError

        allocate(self%quad_pts(nnodes))
        allocate(self%dtdx(nnodes,nterms_s))
        allocate(self%dtdy(nnodes,nterms_s))
        allocate(self%dtdz(nnodes,nterms_s),stat=ierr)
        if (ierr /= 0) call AllocationError

        allocate(self%mass(nterms_s,nterms_s),stat=ierr)
        allocate(self%invmass(nterms_s,nterms_s),stat=ierr)
        if (ierr /= 0) call AllocationError

        call self%compute_quadrature_metrics()                  !> Compute element metrics
        call self%compute_element_matrices()                    !> Compute mass matrices and derivative matrices

        self%numInitialized = .true.                            !> Confirm element numerics were initialized
    end subroutine







    !>  Assign quadrature instances for solution modes (GQ) and mesh modes (GQMESH)
    !!
    !!  @author Nathan A. Wukie
    !---------------------------------------
    subroutine assign_quadrature(self)
        use mod_quadrature,     only: compute_nnodes_gq
        class(element_t),   intent(inout)   :: self

        integer(ik) :: nterms_s,nterms_c
        integer(ik) :: nnodes_face, nnodes_vol, igq, igq_s, igq_f
        logical     :: has_correct_nodes_terms

        nterms_s = self%nterms_s
        nterms_c = self%nterms_c

        if (nterms_c == 0) call signal(FATAL,'element%assign_quadrature -- coordinate expansion not defined')

        call compute_nnodes_gq(nterms_s,nterms_c,nnodes_face,nnodes_vol)    !> Get number of quadrature nodes


        ! Get solution quadrature instance
        call get_quadrature(nterms_s,nnodes_vol,nnodes_face,igq_s)
        self%gq => GQ(igq_s)

        ! Get coordinate quadrature instance
        call get_quadrature(nterms_c,nnodes_vol,nnodes_face,igq_f)
        self%gqmesh => GQ(igq_f)


    end subroutine







    !> Compute element metric and jacobian terms
    !!
    !!  @author Nathan A. Wukie
    !-------------------------------------------------------------------------
    subroutine compute_quadrature_metrics(self)
        class(element_t),    intent(inout)   :: self

        integer(ik)             :: inode
        integer(ik)             :: nnodes

        real(rk)    :: dxdxi(self%gq%vol%nnodes), dxdeta(self%gq%vol%nnodes), dxdzeta(self%gq%vol%nnodes)
        real(rk)    :: dydxi(self%gq%vol%nnodes), dydeta(self%gq%vol%nnodes), dydzeta(self%gq%vol%nnodes)
        real(rk)    :: dzdxi(self%gq%vol%nnodes), dzdeta(self%gq%vol%nnodes), dzdzeta(self%gq%vol%nnodes)

        nnodes = self%gq%vol%nnodes
        ! compute element metric terms
        dxdxi   = matmul(self%gqmesh%vol%ddxi,  self%coords%mat(:,1))
        dxdeta  = matmul(self%gqmesh%vol%ddeta, self%coords%mat(:,1))
        dxdzeta = matmul(self%gqmesh%vol%ddzeta,self%coords%mat(:,1))

        dydxi   = matmul(self%gqmesh%vol%ddxi,  self%coords%mat(:,2))
        dydeta  = matmul(self%gqmesh%vol%ddeta, self%coords%mat(:,2))
        dydzeta = matmul(self%gqmesh%vol%ddzeta,self%coords%mat(:,2))

        dzdxi   = matmul(self%gqmesh%vol%ddxi,  self%coords%mat(:,3))
        dzdeta  = matmul(self%gqmesh%vol%ddeta, self%coords%mat(:,3))
        dzdzeta = matmul(self%gqmesh%vol%ddzeta,self%coords%mat(:,3))

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

        !> compute inverse cell mapping jacobian
        self%jinv = dxdxi*dydeta*dzdzeta - dxdeta*dydxi*dzdzeta - &
                    dxdxi*dydzeta*dzdeta + dxdzeta*dydxi*dzdeta + &
                    dxdeta*dydzeta*dzdxi - dxdzeta*dydeta*dzdxi
    end subroutine





    !> Subroutine computes element-specific matrices
    !!      - Mass matrix   (mass, invmass)
    !!      - Matrices of cartesian gradients of basis/test functions (dtdx, dtdy, dtdz)
    !!      - Cartesian coordinates of quadrature points (quad_pts)
    !! @author Nathan A. Wukie
    !---------------------------------------------------------------------
    subroutine compute_element_matrices(self)
        class(element_t),   intent(inout)   :: self

        !> Call to compute matrices of cartesian gradients at each quadrature node
        call self%compute_gradients_cartesian()

        !> Call to compute mass matrix
        call self%compute_mass_matrix()

        !> Call to compute cartesian coordinates at each quadrature node
        call self%compute_quadrature_coords()
    end subroutine






    !> Compute matrices containing cartesian gradients of basis/test function
    !! at each quadrature node.
    !!
    !! @author Nathan A. Wukie
    !-----------------------------------------------------------------------------------
    subroutine compute_gradients_cartesian(self)
        class(element_t),   intent(inout)   :: self
        integer(ik)                         :: iterm,inode

        do iterm = 1,self%nterms_s
            do inode = 1,self%gq%vol%nnodes
                self%dtdx(inode,iterm) = self%metric(1,1,inode) * self%gq%vol%ddxi(inode,iterm)   * (ONE/self%jinv(inode)) + &
                                         self%metric(2,1,inode) * self%gq%vol%ddeta(inode,iterm)  * (ONE/self%jinv(inode)) + &
                                         self%metric(3,1,inode) * self%gq%vol%ddzeta(inode,iterm) * (ONE/self%jinv(inode))

                self%dtdy(inode,iterm) = self%metric(1,2,inode) * self%gq%vol%ddxi(inode,iterm)   * (ONE/self%jinv(inode)) + &
                                         self%metric(2,2,inode) * self%gq%vol%ddeta(inode,iterm)  * (ONE/self%jinv(inode)) + &
                                         self%metric(3,2,inode) * self%gq%vol%ddzeta(inode,iterm) * (ONE/self%jinv(inode))

                self%dtdz(inode,iterm) = self%metric(1,3,inode) * self%gq%vol%ddxi(inode,iterm)   * (ONE/self%jinv(inode)) + &
                                         self%metric(2,3,inode) * self%gq%vol%ddeta(inode,iterm)  * (ONE/self%jinv(inode)) + &
                                         self%metric(3,3,inode) * self%gq%vol%ddzeta(inode,iterm) * (ONE/self%jinv(inode))
            end do
        end do

    end subroutine






    !> Compute cartesian coordinates at each quadrature point
    !!
    !! @author Nathan A. Wukie
    !-----------------------------------------------------------------------
    subroutine compute_quadrature_coords(self)
        class(element_t),   intent(inout)   :: self
        integer(ik)                         :: nnodes
        real(rk)                            :: x(self%gq%vol%nnodes),y(self%gq%vol%nnodes),z(self%gq%vol%nnodes)
        integer(ik)                         :: inode

        nnodes = self%gq%vol%nnodes
        !> compute cartesian coordinates associated with quadrature points
        x = matmul(self%gqmesh%vol%val,self%coords%mat(:,1))
        y = matmul(self%gqmesh%vol%val,self%coords%mat(:,2))
        z = matmul(self%gqmesh%vol%val,self%coords%mat(:,3))

        !> Initialize each point with cartesian coordinates
        do inode = 1,nnodes
            call self%quad_pts(inode)%set(x(inode),y(inode),z(inode))
        end do
    end subroutine






    !>  Compute element-local mass matrix
    !!
    !!  @author Nathan A. Wukie
    !---------------------------------------------------------------------
    subroutine compute_mass_matrix(self)
        class(element_t), intent(inout) :: self
        integer(ik)  :: iterm,i,j
        real(rk)     :: temp(self%nterms_s,self%gq%vol%nnodes)

        self%invmass = 0._rk
        self%mass = 0._rk
        temp = transpose(self%gq%vol%val)

        ! Multiply rows by quadrature weights and cell jacobians
        do iterm = 1,self%nterms_s
            temp(iterm,:) = temp(iterm,:)*(self%gq%vol%weights)*(self%jinv)
        end do

        ! Perform the matrix multiplication of the transpose val matrix by
        ! the standard matrix. This produces the mass matrix. I think...
        self%mass = matmul(temp,self%gq%vol%val)

        ! Compute and store the inverted mass matrix
        self%invmass = inv(self%mass)

    end subroutine ! compute_mass_matrix











    !=============================================================================
    !
    !   Integrate volume flux.
    !
    !   Multiplies volume fluxes by gradient of test functions and integrates.
    !
    !=============================================================================
    subroutine integrate_volume_flux(self,flux_x,flux_y,flux_z,varindex)
        class(element_t),   intent(inout)  :: self
        type(AD_D),         intent(inout)  :: flux_x(:), flux_y(:), flux_z(:)
        integer(ik),        intent(in)     :: varindex

        type(AD_D),  dimension(self%nterms_s) :: integral

        ! Multiply each component by quadrature weights and jacobian
        flux_x = (flux_x) * (self%gq%vol%weights) * (self%jinv)
        flux_y = (flux_y) * (self%gq%vol%weights) * (self%jinv)
        flux_z = (flux_z) * (self%gq%vol%weights) * (self%jinv)

        ! Multiply by column of test function gradients, integrate, and add to RHS
        integral = matmul(transpose(self%dtdx),flux_x)
!        self%rhs%vals(:,varindex) = self%rhs%vals(:,varindex) + integral

        integral = matmul(transpose(self%dtdy),flux_y)
!        self%rhs%vals(:,varindex) = self%rhs%vals(:,varindex) + integral

        integral = matmul(transpose(self%dtdz),flux_z)
!        self%rhs%vals(:,varindex) = self%rhs%vals(:,varindex) + integral

    end subroutine


    !=============================================================================
    !
    !   Integrate volume source.
    !
    !   Multiplies volume source by test functions and integrates.
    !
    !=============================================================================
    subroutine integrate_volume_source(self,source,varindex)
        class(element_t),   intent(inout)  :: self
        type(AD_D),         intent(inout)  :: source(:)
        integer(ik),        intent(in)     :: varindex

        type(AD_D),  dimension(self%nterms_s)  :: integral

        ! Multiply by weights
        source = (source) * (self%gq%vol%weights) * (self%jinv)

        ! Multiply by column of test functions, integrate, and contribute to RHS
        integral = matmul(transpose(self%gq%vol%val),source)
!        self%rhs%vars(:,varindex) = self%rhs%vars(:,varindex) - integral


    end subroutine








    !> Convert local(xi,eta,zeta) coordinates to global coordinates(x,y,z)
    !!
    !----------------------------------------------------------------
    function x(self,xi,eta,zeta) result(xval)
        class(element_t),   intent(in)  :: self
        real(rk),      intent(in)  :: xi,eta,zeta
        real(rk)                   :: xval

        type(point_t)              :: node
        real(rk)                   :: polyvals(self%nterms_c)
        integer(ik)                :: iterm

        call node%set(xi,eta,zeta)

        ! Evaluate polynomial modes at node location
        do iterm = 1,self%nterms_c
            polyvals(iterm)  = polynomialVal(3,self%nterms_c,iterm,node)
        end do

        ! Evaluate x from dot product of modes and polynomial values
        xval = dot_product(self%coords%mat(:,1),polyvals)

    end function

    function y(self,xi,eta,zeta) result(yval)
        class(element_t),   intent(in)  :: self
        real(rk),           intent(in)  :: xi,eta,zeta
        real(rk)                        :: yval

        type(point_t)              :: node
        real(rk)                   :: polyvals(self%nterms_c)
        integer(ik)                :: iterm

        call node%set(xi,eta,zeta)

        ! Evaluate polynomial modes at node location
        do iterm = 1,self%nterms_c
            polyvals(iterm)  = polynomialVal(3,self%nterms_c,iterm,node)
        end do

        ! Evaluate x from dot product of modes and polynomial values
        yval = dot_product(self%coords%mat(:,2),polyvals)

    end function

    function z(self,xi,eta,zeta) result(zval)
        class(element_t),   intent(in)  :: self
        real(rk),           intent(in)  :: xi,eta,zeta
        real(rk)                        :: zval

        type(point_t)              :: node
        real(rk)                   :: polyvals(self%nterms_c)
        integer(ik)                :: iterm

        call node%set(xi,eta,zeta)

        ! Evaluate polynomial modes at node location
        do iterm = 1,self%nterms_c
            polyvals(iterm)  = polynomialVal(3,self%nterms_c,iterm,node)
        end do

        ! Evaluate x from dot product of modes and polynomial values
        zval = dot_product(self%coords%mat(:,3),polyvals)

    end function










    subroutine destructor(self)
        type(element_t), intent(inout) :: self

!> Shouldn't need to deallocate an 'allocatable'. The compiler is supposed to do that for you
!        if (allocated(self%quad_pts))   deallocate(self%quad_pts)
!        if (allocated(self%elem_pts))   deallocate(self%elem_pts)
!        if (allocated(self%metric))     deallocate(self%metric)
!        if (allocated(self%jinv))       deallocate(self%jinv)
!        if (allocated(self%dtdx))       deallocate(self%dtdx)
!        if (allocated(self%dtdy))       deallocate(self%dtdy)
!        if (allocated(self%dtdz))       deallocate(self%dtdz)
!        if (allocated(self%mass))       deallocate(self%mass)
!        if (allocated(self%invmass))    deallocate(self%invmass)

    end subroutine

end module type_element
