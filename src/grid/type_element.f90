module type_element
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,XI_MIN,XI_MAX,ETA_MIN, &
                                      ETA_MAX,ZETA_MIN,ZETA_MAX,ONE,ZERO, &
                                      X_DIR, Y_DIR, Z_DIR, XI_DIR, ETA_DIR, ZETA_DIR, &
                                      TWO_DIM, THREE_DIM, RKTOL, VALID_POINT, INVALID_POINT
    use mod_quadrature,         only: GQ, get_quadrature
    use mod_grid,               only: get_element_mapping, face_corners
    use mod_polynomial,         only: polynomialVal, dpolynomialVal
    use mod_inv,                only: inv


    use type_point,                 only: point_t
    use type_densevector,           only: densevector_t
    use type_quadrature,            only: quadrature_t
    use type_function,              only: function_t
    use type_element_connectivity,  only: element_connectivity_t
    use DNAD_D
    implicit none




    !>  Element data type
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
    !!  @author Mayank Sharma
    !!  @date   11/12/2016
    !!
    !-----------------------------------------------------------------------------------------
    type, public :: element_t

        ! Element info
        integer(ik)     :: idomain_g                        !< Global index of parent domain
        integer(ik)     :: idomain_l                        !< Proc-local index of parent domain
        integer(ik)     :: ielement_g                       !< Domain-global index of element
        integer(ik)     :: ielement_l                       !< Proc-local index of the element

        integer(ik)     :: spacedim                         !< Number of spatial dimensions for the element
        integer(ik)     :: neqns                            !< Number of equations being solved
        integer(ik)     :: nterms_s                         !< Number of terms in solution expansion.  
        integer(ik)     :: nterms_c                         !< Number of terms in coordinate expansion. 
        integer(ik)     :: ntime                            !< Number of time levels in solution

        ! Element quadrature points, mesh points and modes
        type(element_connectivity_t)    :: connectivity         !< Connectivity list. Integer indices of the associated nodes in block node list
        type(point_t), allocatable      :: quad_pts(:)          !< Cartesian coordinates of discrete quadrature points
        type(point_t), allocatable      :: elem_pts(:)          !< Cartesian coordinates of discrete points defining element
        type(densevector_t)             :: coords               !< Modal representation of cartesian coordinates (nterms_var,(x,y,z))

        ! Element metric terms
        real(rk), allocatable           :: metric(:,:,:)        !< metric matrix for each quadrature node    (mat_i,mat_j,quad_pt)
        real(rk), allocatable           :: jinv(:)              !< jacobian terms at quadrature nodes

        ! Matrices of cartesian gradients of basis/test functions
        real(rk), allocatable           :: ddx(:,:)             !< Derivative of basis functions in x-direction at quadrature nodes
        real(rk), allocatable           :: ddy(:,:)             !< Derivative of basis functions in y-direction at quadrature nodes
        real(rk), allocatable           :: ddz(:,:)             !< Derivative of basis functions in z-direction at quadrature nodes
        real(rk), allocatable           :: ddx_trans(:,:)       !< Derivative of basis functions in x-direction at quadrature nodes - transposed
        real(rk), allocatable           :: ddy_trans(:,:)       !< Derivative of basis functions in y-direction at quadrature nodes - transposed
        real(rk), allocatable           :: ddz_trans(:,:)       !< Derivative of basis functions in z-direction at quadrature nodes - transposed

        ! Quadrature matrices
        type(quadrature_t), pointer     :: gq     => null()     !< Pointer to instance for solution expansion
        type(quadrature_t), pointer     :: gqmesh => null()     !< Pointer to instance for coordinate expansion

        ! Element-local mass, inverse mass matrices
        real(rk), allocatable           :: mass(:,:)        
        real(rk), allocatable           :: invmass(:,:)

        ! Element volume
        real(rk)                        :: vol


        ! A psudo-timestep for each equation in the element. Used in the nonlinear solver. 
        ! Quasi-Newton, for example.
        real(rk),   allocatable         :: dtau(:)

        ! Logical tests
        logical :: geomInitialized = .false.
        logical :: numInitialized  = .false.


    contains

        ! Initialization procedures
        procedure, public   :: init_geom
        procedure, public   :: init_sol


        ! Compute discrete value for x/y/z-coordinate at a given xi,eta,zeta.
        procedure, public   :: x                      
        procedure, public   :: y                      
        procedure, public   :: z                      

        procedure, public   :: grid_point             !< Compute a discrete value for a physical coordinate at a given xi, eta, zeta.
        procedure, public   :: computational_point    !< Compute a discrete value for a computational coordinate at a given x, y, z.
        procedure, public   :: metric_point           !< Compute a discrete value for a metric term at a given xi, eta, zeta.
        procedure, public   :: solution_point         !< Compute a discrete value for the solution at a given xi,eta, zeta.
        procedure, public   :: derivative_point
        procedure, public   :: project                !< Compute a projection of a function onto the solution basis


        ! Get connected face
        procedure, public   :: get_face_from_corners

        ! Private utility procedure
        procedure           :: compute_element_matrices
        procedure           :: compute_mass_matrix
        procedure           :: compute_gradients_cartesian
        procedure           :: compute_quadrature_metrics
        procedure           :: compute_quadrature_coords
        procedure           :: assign_quadrature



        final               :: destructor

    end type element_t
    !******************************************************************************************

    private




contains
    





    !>  Initialize element geometry
    !!      - Set element points
    !!      - Compute modal representation of element cartesian coordinates
    !!      - Compute element metric terms
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in] nterms_c     Number of terms in the modal representation of the 
    !!                          cartesian coordinates.
    !!  @param[in] points       Array of cartesian points defining the element
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init_geom(self,spacedim,nodes,connectivity,idomain_l,ielem_l)
        class(element_t),               intent(inout)   :: self
        integer(ik),                    intent(in)      :: spacedim
        type(point_t),                  intent(in)      :: nodes(:)
        type(element_connectivity_t),   intent(in)      :: connectivity
        integer(ik),                    intent(in)      :: idomain_l
        integer(ik),                    intent(in)      :: ielem_l

        character(:),   allocatable :: user_msg
        type(point_t),  allocatable :: points(:)
        real(rk),       allocatable :: element_mapping(:,:)
        real(rk),       allocatable :: xmodes(:), ymodes(:), zmodes(:)
        integer(ik)                 :: ierr, nterms_c, ipt, npts_1d, npts, &
                                       mapping, inode, idomain_g, ielem_g
        integer(ik)                 :: ntime = 1


        user_msg = "element%init_geom: element already initialized."
        if (self%geomInitialized) call chidg_signal(FATAL,user_msg)


        !
        ! Get connectivity info
        !
        idomain_g = connectivity%get_domain_index()
        ielem_g   = connectivity%get_element_index()
        mapping   = connectivity%get_element_mapping()



        !
        ! Accumulate coordinates for current element from node list.
        !
        npts_1d = mapping+1
        npts    = npts_1d * npts_1d * npts_1d
        allocate(points(npts), stat=ierr)
        if (ierr /= 0) call AllocationError

        do ipt = 1,npts
            !
            ! Get node index
            !
            inode = connectivity%get_element_node(ipt)

            !
            ! Add node to element points list
            !
            points(ipt) = nodes(inode)
        end do !ipt





        !
        ! Get element mapping
        !
        element_mapping = get_element_mapping(spacedim,mapping)
        nterms_c = size(element_mapping,1)
        self%nterms_c = nterms_c

        user_msg = "element%init_geom: mapping and points do not match."
        if (nterms_c /= size(points)) call chidg_signal(FATAL,user_msg)


        !
        ! Allocate storage
        !
        allocate(self%elem_pts(nterms_c),stat=ierr)
        call self%coords%init(nterms_c,3,ntime,idomain_g,idomain_l,ielem_g,ielem_l)
        self%spacedim       = spacedim
        self%idomain_g      = idomain_g
        self%idomain_l      = idomain_l
        self%ielement_g     = ielem_g
        self%ielement_l     = ielem_l
        self%elem_pts       = points
        self%connectivity   = connectivity

        
        !
        ! Compute mesh x,y,z modes
        !
        xmodes = matmul(element_mapping,self%elem_pts(:)%c1_)
        ymodes = matmul(element_mapping,self%elem_pts(:)%c2_)
        zmodes = matmul(element_mapping,self%elem_pts(:)%c3_)

        call self%coords%setvar(1,itime = 1,vals = xmodes)
        call self%coords%setvar(2,itime = 1,vals = ymodes)
        call self%coords%setvar(3,itime = 1,vals = zmodes)


        !
        ! Confirm element geometry was initialized
        !
        self%geomInitialized = .true.   


    end subroutine init_geom
    !******************************************************************************************









    !>  Initialize element numerics
    !!      - Allocate storage for solution and supporting matrices
    !!      - Compute element-local matrices (cartesian gradients, mass matrices)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  nterms_s    Number of terms in the modal representation of the solution
    !!  @param[in]  neqns       Number of equations contained in the element solution
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/12/2016
    !!
    !!
    subroutine init_sol(self,neqns,nterms_s,ntime)
        class(element_t),   intent(inout) :: self
        integer(ik),        intent(in)    :: neqns
        integer(ik),        intent(in)    :: nterms_s
        integer(ik),        intent(in)    :: ntime

        integer(ik) :: ierr
        integer(ik) :: nnodes


        self%nterms_s    = nterms_s     ! Set number of terms in modal expansion of solution
        self%neqns       = neqns        ! Set number of equations being solved
        self%ntime       = ntime        ! Set number of time steps in solution


        call self%assign_quadrature()   ! With nterms_s and nterms_c defined, we can assign a quadrature instance
        nnodes = self%gq%vol%nnodes     ! With a quadrature instance assigned, we have the number of quadrature nodes


        !
        ! (Re)Allocate storage for element data structures
        !
        if (allocated(self%jinv)) deallocate( self%jinv, self%metric, self%quad_pts,            &
                                              self%ddx, self%ddy, self%ddz,                     &
                                              self%ddx_trans, self%ddy_trans, self%ddz_trans,   &
                                              self%mass, self%invmass, self%dtau)
        allocate(self%jinv(nnodes),                 &
                 self%metric(3,3,nnodes),           &
                 self%quad_pts(nnodes),             &
                 self%ddx(nnodes,nterms_s),         &
                 self%ddy(nnodes,nterms_s),         &
                 self%ddz(nnodes,nterms_s),         &
                 self%ddx_trans(nterms_s,nnodes),   &
                 self%ddy_trans(nterms_s,nnodes),   &
                 self%ddz_trans(nterms_s,nnodes),   &
                 self%mass(nterms_s,nterms_s),      &
                 self%invmass(nterms_s,nterms_s),   &
                 self%dtau(neqns), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Call element metric and matrix calculation routines
        !
        call self%compute_quadrature_metrics()          ! Compute element metrics
        call self%compute_element_matrices()            ! Compute mass matrices and derivative matrices



        !
        ! Confirm element numerics were initialized
        !
        self%numInitialized = .true.    

    end subroutine init_sol
    !*****************************************************************************************







    !>  Assign quadrature instances for solution modes (GQ) and mesh modes (GQMESH)
    !!      self%gq
    !!      self%gqmesh
    !!
    !!  TODO: would be good to eliminate pointers in the element data type and just 
    !!        use integer indices to a global array of quadrature instances.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine assign_quadrature(self)
        use mod_quadrature,     only: compute_nnodes_gq
        class(element_t),   intent(inout)   :: self

        character(:), allocatable   :: user_msg
        integer(ik)                 :: nterms_s,nterms_c,spacedim
        integer(ik)                 :: nnodes_face, nnodes_vol, igq_s, igq_f

        spacedim = self%spacedim
        nterms_s = self%nterms_s
        nterms_c = self%nterms_c

        user_msg = "element%assign_quadrature: coordinate expansion not defined."
        if (nterms_c == 0) call chidg_signal(FATAL,user_msg)



        !
        ! Get number of quadrature nodes
        !
        call compute_nnodes_gq(spacedim,nterms_s,nterms_c,nnodes_face,nnodes_vol)


        !
        ! Get solution quadrature instance
        !
        call get_quadrature(spacedim,nterms_s,nnodes_vol,nnodes_face,igq_s)
        self%gq => GQ(igq_s)


        !
        ! Get coordinate quadrature instance
        !
        call get_quadrature(spacedim,nterms_c,nnodes_vol,nnodes_face,igq_f)
        self%gqmesh => GQ(igq_f)


    end subroutine assign_quadrature
    !*****************************************************************************************











    !> Compute element metric and jacobian terms
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !! TODO: Generalized 2D physical coordinates. Currently assumes x-y
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_quadrature_metrics(self)
        class(element_t),    intent(inout)   :: self

        integer(ik)             :: inode
        integer(ik)             :: nnodes

        real(rk)    :: dxdxi(self%gq%vol%nnodes), dxdeta(self%gq%vol%nnodes), dxdzeta(self%gq%vol%nnodes)
        real(rk)    :: dydxi(self%gq%vol%nnodes), dydeta(self%gq%vol%nnodes), dydzeta(self%gq%vol%nnodes)
        real(rk)    :: dzdxi(self%gq%vol%nnodes), dzdeta(self%gq%vol%nnodes), dzdzeta(self%gq%vol%nnodes)


        nnodes = self%gq%vol%nnodes

        !
        ! Compute element metric terms
        !
        dxdxi   = matmul(self%gqmesh%vol%ddxi,  self%coords%getvar(1,itime = 1))
        dxdeta  = matmul(self%gqmesh%vol%ddeta, self%coords%getvar(1,itime = 1))
        dxdzeta = matmul(self%gqmesh%vol%ddzeta,self%coords%getvar(1,itime = 1))

        dydxi   = matmul(self%gqmesh%vol%ddxi,  self%coords%getvar(2,itime = 1))
        dydeta  = matmul(self%gqmesh%vol%ddeta, self%coords%getvar(2,itime = 1))
        dydzeta = matmul(self%gqmesh%vol%ddzeta,self%coords%getvar(2,itime = 1))

        dzdxi   = matmul(self%gqmesh%vol%ddxi,  self%coords%getvar(3,itime = 1))
        dzdeta  = matmul(self%gqmesh%vol%ddeta, self%coords%getvar(3,itime = 1))
        dzdzeta = matmul(self%gqmesh%vol%ddzeta,self%coords%getvar(3,itime = 1))




        !
        ! TODO: Generalized 2D physical coordinates. Currently assumes x-y
        !
        if ( self%spacedim == TWO_DIM ) then
            dzdxi   = ZERO
            dzdeta  = ZERO
            dzdzeta = ONE
        end if





        !
        ! Loop through quadrature nodes and compute metric terms
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
        ! Compute inverse cell mapping jacobian
        !
        self%jinv = dxdxi*dydeta*dzdzeta - dxdeta*dydxi*dzdzeta - &
                    dxdxi*dydzeta*dzdeta + dxdzeta*dydxi*dzdeta + &
                    dxdeta*dydzeta*dzdxi - dxdzeta*dydeta*dzdxi



        !
        ! Compute element volume
        !
        self%vol = abs(sum(self%jinv * self%gq%vol%weights))

    end subroutine compute_quadrature_metrics
    !******************************************************************************************











    !>  Subroutine computes element-specific matrices
    !!      - Mass matrix   (mass, invmass)
    !!      - Matrices of cartesian gradients of basis/test functions (ddx, ddy, ddz)
    !!      - Cartesian coordinates of quadrature points (quad_pts)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_element_matrices(self)
        class(element_t),   intent(inout)   :: self

        !
        ! Call to compute matrices of cartesian gradients at each quadrature node
        !
        call self%compute_gradients_cartesian()

        !
        ! Call to compute mass matrix
        !
        call self%compute_mass_matrix()

        !
        ! Call to compute cartesian coordinates at each quadrature node
        !
        call self%compute_quadrature_coords()

    end subroutine compute_element_matrices
    !******************************************************************************************











    !>  Compute matrices containing cartesian gradients of basis/test function
    !!  at each quadrature node.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_gradients_cartesian(self)
        class(element_t),   intent(inout)   :: self
        integer(ik)                         :: iterm,inode

        do iterm = 1,self%nterms_s
            do inode = 1,self%gq%vol%nnodes
                self%ddx(inode,iterm) = self%metric(1,1,inode) * self%gq%vol%ddxi(inode,iterm)   * (ONE/self%jinv(inode)) + &
                                        self%metric(2,1,inode) * self%gq%vol%ddeta(inode,iterm)  * (ONE/self%jinv(inode)) + &
                                        self%metric(3,1,inode) * self%gq%vol%ddzeta(inode,iterm) * (ONE/self%jinv(inode))

                self%ddy(inode,iterm) = self%metric(1,2,inode) * self%gq%vol%ddxi(inode,iterm)   * (ONE/self%jinv(inode)) + &
                                        self%metric(2,2,inode) * self%gq%vol%ddeta(inode,iterm)  * (ONE/self%jinv(inode)) + &
                                        self%metric(3,2,inode) * self%gq%vol%ddzeta(inode,iterm) * (ONE/self%jinv(inode))

                self%ddz(inode,iterm) = self%metric(1,3,inode) * self%gq%vol%ddxi(inode,iterm)   * (ONE/self%jinv(inode)) + &
                                        self%metric(2,3,inode) * self%gq%vol%ddeta(inode,iterm)  * (ONE/self%jinv(inode)) + &
                                        self%metric(3,3,inode) * self%gq%vol%ddzeta(inode,iterm) * (ONE/self%jinv(inode))
            end do
        end do


        self%ddx_trans = transpose(self%ddx)
        self%ddy_trans = transpose(self%ddy)
        self%ddz_trans = transpose(self%ddz)

    end subroutine compute_gradients_cartesian
    !******************************************************************************************











    !>  Compute cartesian coordinates at each quadrature point
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_quadrature_coords(self)
        class(element_t),   intent(inout)   :: self
        integer(ik)                         :: nnodes
        real(rk)                            :: x(self%gq%vol%nnodes),y(self%gq%vol%nnodes),z(self%gq%vol%nnodes)
        integer(ik)                         :: inode

        nnodes = self%gq%vol%nnodes

        !
        ! compute cartesian coordinates associated with quadrature points
        !
        x = matmul(self%gqmesh%vol%val,self%coords%getvar(1,itime = 1))
        y = matmul(self%gqmesh%vol%val,self%coords%getvar(2,itime = 1))
        z = matmul(self%gqmesh%vol%val,self%coords%getvar(3,itime = 1))


        !
        ! Initialize each point with cartesian coordinates
        !
        do inode = 1,nnodes
            call self%quad_pts(inode)%set(x(inode),y(inode),z(inode))
        end do

    end subroutine compute_quadrature_coords
    !*****************************************************************************************













    !>  Compute element-local mass matrix
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_mass_matrix(self)
        class(element_t), intent(inout) :: self
        integer(ik)  :: iterm
        real(rk)     :: temp(self%nterms_s,self%gq%vol%nnodes)

        self%invmass = ZERO
        self%mass    = ZERO
        temp = transpose(self%gq%vol%val)



        !
        ! Multiply rows by quadrature weights and cell jacobians
        !
        do iterm = 1,self%nterms_s
            temp(iterm,:) = temp(iterm,:)*(self%gq%vol%weights)*(self%jinv)
        end do


        !
        ! Perform the matrix multiplication of the transpose val matrix by
        ! the standard matrix. This produces the mass matrix. I think...
        !
        self%mass = matmul(temp,self%gq%vol%val)



        !
        ! Compute and store the inverted mass matrix
        !
        self%invmass = inv(self%mass)



    end subroutine compute_mass_matrix
    !******************************************************************************************
















    !>  Convert local(xi,eta,zeta) coordinates to global coordinates(x,y,z)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !------------------------------------------------------------------------------------------
    function x(self,xi,eta,zeta) result(xval)
        class(element_t),   intent(in)  :: self
        real(rk),      intent(in)  :: xi,eta,zeta
        real(rk)                   :: xval

        type(point_t)              :: node
        real(rk)                   :: polyvals(self%nterms_c)
        integer(ik)                :: iterm, spacedim

        call node%set(xi,eta,zeta)

        spacedim = self%spacedim

        !
        ! Evaluate polynomial modes at node location
        !
        do iterm = 1,self%nterms_c

            polyvals(iterm)  = polynomialVal(spacedim,self%nterms_c,iterm,node)

        end do

        
        !
        ! Evaluate x from dot product of modes and polynomial values
        !
        xval = dot_product(self%coords%getvar(1,itime = 1),polyvals)

    end function x
    !******************************************************************************************







    !>  Convert local(xi,eta,zeta) coordinates to global coordinates(x,y,z)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !-----------------------------------------------------------------------------------------
    function y(self,xi,eta,zeta) result(yval)
        class(element_t),   intent(in)  :: self
        real(rk),           intent(in)  :: xi,eta,zeta
        real(rk)                        :: yval

        type(point_t)              :: node
        real(rk)                   :: polyvals(self%nterms_c)
        integer(ik)                :: iterm, spacedim

        call node%set(xi,eta,zeta)

        spacedim = self%spacedim

        !
        ! Evaluate polynomial modes at node location
        !
        do iterm = 1,self%nterms_c

            polyvals(iterm)  = polynomialVal(spacedim,self%nterms_c,iterm,node)

        end do


        !
        ! Evaluate x from dot product of modes and polynomial values
        !
        yval = dot_product(self%coords%getvar(2,itime = 1),polyvals)

    end function y
    !******************************************************************************************






    !>  Convert local(xi,eta,zeta) coordinates to global coordinates(x,y,z)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !-----------------------------------------------------------------------------------------
    function z(self,xi,eta,zeta) result(zval)
        class(element_t),   intent(in)  :: self
        real(rk),           intent(in)  :: xi,eta,zeta
        real(rk)                        :: zval

        type(point_t)              :: node
        real(rk)                   :: polyvals(self%nterms_c)
        integer(ik)                :: iterm, spacedim

        call node%set(xi,eta,zeta)

        spacedim = self%spacedim

        !
        ! Evaluate polynomial modes at node location
        !
        do iterm = 1,self%nterms_c

            polyvals(iterm)  = polynomialVal(spacedim,self%nterms_c,iterm,node)

        end do


        !
        ! Evaluate x from dot product of modes and polynomial values
        !
        zval = dot_product(self%coords%getvar(3,itime = 1),polyvals)

    end function z
    !******************************************************************************************












    !>  Compute a coordinate value, based on the location in reference space (xi, eta, zeta)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  elem    Element containing coordinate expansion
    !!  @param[in]  icoord  Integer corresponding to coordinate index 1(x), 2(y), 3(z)
    !!  @param[in]  xi      Real value for xi-coordinate
    !!  @param[in]  eta     Real value for eta-coordinate
    !!  @param[in]  zeta    Real value for zeta-coordinate
    !!
    !!  @author Mayank Sharma + MAtteo Ugolotti
    !!  @date   11/5/2016
    !!
    !-----------------------------------------------------------------------------------------
    function grid_point(self,icoord,xi,eta,zeta) result(val)
        class(element_t),   intent(in)  :: self
        integer(ik),        intent(in)  :: icoord
        real(rk),           intent(in)  :: xi, eta, zeta

        real(rk)        :: val
        type(point_t)   :: node
        real(rk)        :: polyvals(self%nterms_c)
        integer(ik)     :: iterm, spacedim

        if (icoord > 3)                 call chidg_signal(FATAL,"element%grid_point -- icoord exceeded 3 physical coordinates")
        if (.not. self%geomInitialized) call chidg_signal(FATAL,"element%grid_point: geometry not initialized")


        call node%set(xi,eta,zeta)

        spacedim = self%spacedim


        !
        ! Evaluate polynomial modes at node location
        !
        do iterm = 1,self%nterms_c
                polyvals(iterm) = polynomialVal(spacedim,self%nterms_c,iterm,node)
        end do


        !
        ! Evaluate mesh point from dot product of modes and polynomial values
        !
        val = dot_product(self%coords%getvar(icoord,itime = 1), polyvals)

    end function grid_point
    !******************************************************************************************










    !> Compute coordinate metric term at a given point in computational space
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  elem        element_t containing the geometry definition and data
    !!  @param[in]  cart_dir    Cartesian coordinate being differentiated
    !!  @param[in]  comp_dir    Computational coordinate being differentiated with respect to
    !!  @param[in]  xi          Computational coordinate - xi
    !!  @param[in]  eta         Computational coordinate - eta
    !!  @param[in]  zeta        Computational coordinate - zeta
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !-----------------------------------------------------------------------------------------
    function metric_point(self,cart_dir,comp_dir,xi,eta,zeta) result(val)
        class(element_t),   intent(in)  :: self
        integer(ik),        intent(in)  :: cart_dir
        integer(ik),        intent(in)  :: comp_dir
        real(rk),           intent(in)  :: xi, eta, zeta
        
        real(rk)        :: val
        type(point_t)   :: node
        real(rk)        :: polyvals(self%nterms_c)
        integer(ik)     :: iterm, spacedim


        if (cart_dir > 3) call chidg_signal(FATAL,"Error: metric_point -- card_dir exceeded 3 physical coordinates")
        if (comp_dir > 3) call chidg_signal(FATAL,"Error: metric_point -- comp_dir exceeded 3 physical coordinates")

        call node%set(xi,eta,zeta)

        spacedim = self%spacedim

        !
        ! Evaluate polynomial modes at node location
        !
        do iterm = 1,self%nterms_c
            polyvals(iterm) = dpolynomialVal(spacedim,self%nterms_c,iterm,node,comp_dir)
        end do


        !
        ! Evaluate mesh point from dot product of modes and polynomial values
        !
        val = dot_product(self%coords%getvar(cart_dir,itime = 1), polyvals)



        !
        ! 2D/3D. For metric terms, unlike solution derivatives, dzdzeta is 1 for 2D, 0 else.
        !
        if ( spacedim == TWO_DIM ) then
            if      ( (cart_dir == X_DIR) .and. (comp_dir == ZETA_DIR) ) then
                val = ZERO
            else if ( (cart_dir == Y_DIR) .and. (comp_dir == ZETA_DIR) ) then
                val = ZERO
            else if ( (cart_dir == Z_DIR) .and. (comp_dir == ZETA_DIR) ) then
                val = ONE
            end if
        end if



    end function metric_point
    !*****************************************************************************************









    !>  Compute a variable value, based on the location in reference space (xi, eta, zeta)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  elem    Element that the solution expansion is associated with.
    !!  @param[in]  q       Solution expansion for a given element.
    !!  @param[in]  ivar    Integer corresponding to variable index.
    !!  @param[in]  xi      Real value for xi-coordinate.
    !!  @param[in]  eta     Real value for eta-coordinate.
    !!  @param[in]  zeta    Real value for zeta-coordinate.
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    function solution_point(self,q,ivar,itime,xi,eta,zeta) result(val)
        class(element_t),       intent(in)      :: self
        class(densevector_t),   intent(in)      :: q
        integer(ik),            intent(in)      :: ivar
        integer(ik),            intent(in)      :: itime
        real(rk),               intent(in)      :: xi,eta,zeta

        real(rk)                   :: val
        type(point_t)              :: node
        real(rk)                   :: polyvals(q%nterms())
        integer(ik)                :: iterm, spacedim


        call node%set(xi,eta,zeta)

        spacedim = self%spacedim


        !
        ! Evaluate polynomial modes at node location
        !
        do iterm = 1,q%nterms()
            polyvals(iterm)  = polynomialVal(spacedim,q%nterms(),iterm,node)
        end do


        !
        ! Evaluate x from dot product of modes and polynomial values
        !
        val = dot_product(q%getvar(ivar,itime),polyvals)

    end function solution_point
    !******************************************************************************************




    !>  Compute a variable value, based on the location in reference space (xi, eta, zeta)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  elem    Element that the solution expansion is associated with.
    !!  @param[in]  q       Solution expansion for a given element.
    !!  @param[in]  ivar    Integer corresponding to variable index.
    !!  @param[in]  xi      Real value for xi-coordinate.
    !!  @param[in]  eta     Real value for eta-coordinate.
    !!  @param[in]  zeta    Real value for zeta-coordinate.
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !----------------------------------------------------------------------------------------
    function derivative_point(self,q,ivar,itime,xi,eta,zeta,dir) result(val)
        class(element_t),       intent(in)      :: self
        class(densevector_t),   intent(in)      :: q
        integer(ik),            intent(in)      :: ivar
        integer(ik),            intent(in)      :: itime
        real(rk),               intent(in)      :: xi,eta,zeta
        integer(ik),            intent(in)      :: dir

        real(rk)        :: val
        type(point_t)   :: node
        real(rk)        :: ddxi(q%nterms()), ddeta(q%nterms()), ddzeta(q%nterms()), &
                           deriv(q%nterms())
        real(rk)        :: metric(3,3), jinv, dxi_dx, dxi_dy, dxi_dz, &
                           deta_dx, deta_dy, deta_dz, dzeta_dx, dzeta_dy, dzeta_dz
        integer(ik)     :: iterm, spacedim


        call node%set(xi,eta,zeta)

        spacedim = self%spacedim



        !
        ! Evaluate polynomial mode derivatives at node location
        !
        do iterm = 1,q%nterms()
            ddxi(iterm)   = DpolynomialVal(spacedim,q%nterms(),iterm,node,XI_DIR)
            ddeta(iterm)  = DpolynomialVal(spacedim,q%nterms(),iterm,node,ETA_DIR)
            ddzeta(iterm) = DpolynomialVal(spacedim,q%nterms(),iterm,node,ZETA_DIR)
        end do


        !
        ! Compute metrics at node
        !
        metric(1,1) = self%metric_point(X_DIR,XI_DIR,  xi,eta,zeta)
        metric(2,1) = self%metric_point(Y_DIR,XI_DIR,  xi,eta,zeta)
        metric(3,1) = self%metric_point(Z_DIR,XI_DIR,  xi,eta,zeta)
        metric(1,2) = self%metric_point(X_DIR,ETA_DIR, xi,eta,zeta)
        metric(2,2) = self%metric_point(Y_DIR,ETA_DIR, xi,eta,zeta)
        metric(3,2) = self%metric_point(Z_DIR,ETA_DIR, xi,eta,zeta)
        metric(1,3) = self%metric_point(X_DIR,ZETA_DIR,xi,eta,zeta)
        metric(2,3) = self%metric_point(Y_DIR,ZETA_DIR,xi,eta,zeta)
        metric(3,3) = self%metric_point(Z_DIR,ZETA_DIR,xi,eta,zeta)


        !
        ! Compute inverse cell mapping jacobian
        !
        jinv = metric(1,1)*metric(2,2)*metric(3,3) - metric(1,2)*metric(2,1)*metric(3,3) - &
               metric(1,1)*metric(2,3)*metric(3,2) + metric(1,3)*metric(2,1)*metric(3,2) + &
               metric(1,2)*metric(2,3)*metric(3,1) - metric(1,3)*metric(2,2)*metric(3,1)





        do iterm = 1,self%nterms_s
            if (dir == X_DIR) then
                dxi_dx   = metric(2,2)*metric(3,3) - metric(2,3)*metric(3,2)
                deta_dx  = metric(2,3)*metric(3,1) - metric(2,1)*metric(3,3)
                dzeta_dx = metric(2,1)*metric(3,2) - metric(2,2)*metric(3,1)
                deriv(iterm) = dxi_dx   * ddxi(iterm)   * (ONE/jinv) + &
                               deta_dx  * ddeta(iterm)  * (ONE/jinv) + &
                               dzeta_dx * ddzeta(iterm) * (ONE/jinv)
            else if (dir == Y_DIR) then
                dxi_dy   = metric(1,3)*metric(3,2) - metric(1,2)*metric(3,3)
                deta_dy  = metric(1,1)*metric(3,3) - metric(1,3)*metric(3,1)
                dzeta_dy = metric(1,2)*metric(3,1) - metric(1,1)*metric(3,2)
                deriv(iterm) = dxi_dy   * ddxi(iterm)   * (ONE/jinv) + &
                               deta_dy  * ddeta(iterm)  * (ONE/jinv) + &
                               dzeta_dy * ddzeta(iterm) * (ONE/jinv)
            else if (dir == Z_DIR) then
                dxi_dz   = metric(1,2)*metric(2,3) - metric(1,3)*metric(2,2)
                deta_dz  = metric(1,3)*metric(2,1) - metric(1,1)*metric(2,3)
                dzeta_dz = metric(1,1)*metric(2,2) - metric(1,2)*metric(2,1)
                deriv(iterm) = dxi_dz   * ddxi(iterm)   * (ONE/jinv) + &
                               deta_dz  * ddeta(iterm)  * (ONE/jinv) + &
                               dzeta_dz * ddzeta(iterm) * (ONE/jinv)
            else
                call chidg_signal(FATAL,"element%derivative_point: Invalid value for 'dir' parameter. 'ddx', 'ddy', 'ddz'.")
            end if
        end do



        !
        ! Evaluate x from dot product of modes and polynomial values
        !
        val = dot_product(q%getvar(ivar,itime),deriv)

    end function derivative_point
    !*****************************************************************************************







    
    !>  Compute a computational location(xi,eta,zeta), based on the location in cartesian space (x,y,z)
    !!
    !!  NOTE: Will return a location, even if the newton solve did not converge. So make
    !!        sure to check the 'status' component of the returned point_t to check if 
    !!        the point is valid.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/23/2016
    !!
    !!  @param[in]  coord   Computational coordinate that is getting computed (xi, eta, zeta)
    !!  @param[in]  x       Real value for x-coordinate.
    !!  @param[in]  y       Real value for y-coordinate.
    !!  @param[in]  z       Real value for z-coordinate.
    !!
    !-----------------------------------------------------------------------------------------
    function computational_point(self,x,y,z) result(loc)
        class(element_t),   intent(in)  :: self
        real(rk),           intent(in)  :: x
        real(rk),           intent(in)  :: y
        real(rk),           intent(in)  :: z

        type(point_t)       :: loc
        real(rk)            :: xi,  eta, zeta,   &
                               xn,  yn,  zn
        integer(ik)         :: inewton

        real(rk)    :: mat(3,3), minv(3,3)
        real(rk)    :: R(3)
        real(rk)    :: dcoord(3)
        real(rk)    :: res, tol


        tol = 1000._rk*RKTOL
        !tol = RKTOL


        !
        ! Newton iteration to find the donor local coordinates
        !
        xi   = ZERO
        eta  = ZERO
        zeta = ZERO
        do inewton = 1,20

            !
            ! Compute local cartesian coordinates as a function of xi,eta,zeta
            !
            xn = self%x(xi,eta,zeta)
            yn = self%y(xi,eta,zeta)
            zn = self%z(xi,eta,zeta)


            !
            ! Assemble residual vector
            !
            R(1) = -(xn - x)
            R(2) = -(yn - y)
            R(3) = -(zn - z)


            !
            ! Assemble coordinate jacobian matrix
            !
            mat(1,1) = self%metric_point(X_DIR,XI_DIR,  xi,eta,zeta)
            mat(2,1) = self%metric_point(Y_DIR,XI_DIR,  xi,eta,zeta)
            mat(3,1) = self%metric_point(Z_DIR,XI_DIR,  xi,eta,zeta)
            mat(1,2) = self%metric_point(X_DIR,ETA_DIR, xi,eta,zeta)
            mat(2,2) = self%metric_point(Y_DIR,ETA_DIR, xi,eta,zeta)
            mat(3,2) = self%metric_point(Z_DIR,ETA_DIR, xi,eta,zeta)
            mat(1,3) = self%metric_point(X_DIR,ZETA_DIR,xi,eta,zeta)
            mat(2,3) = self%metric_point(Y_DIR,ZETA_DIR,xi,eta,zeta)
            mat(3,3) = self%metric_point(Z_DIR,ZETA_DIR,xi,eta,zeta)


            !
            ! Invert jacobian matrix
            !
            minv = inv(mat)


            !
            ! Compute coordinate update
            !
            dcoord = matmul(minv,R)


            !
            ! Update coordinates
            !
            xi   = xi   + dcoord(1)
            eta  = eta  + dcoord(2)
            zeta = zeta + dcoord(3)


            !
            ! Compute residual coordinate norm
            !
            res = norm2(R)


            !
            ! Exit if converged
            !
            if ( res < tol ) then
                loc%status = VALID_POINT  ! point found
                call loc%set(xi,eta,zeta)
                exit
            end if


            !
            ! Limit computational coordinates, in case they go out of bounds.
            !
            if ( xi   >  ONE ) xi   =  ONE
            if ( xi   < -ONE ) xi   = -ONE
            if ( eta  >  ONE ) eta  =  ONE
            if ( eta  < -ONE ) eta  = -ONE
            if ( zeta >  ONE ) zeta =  ONE
            if ( zeta < -ONE ) zeta = -ONE


            if ( inewton == 20 ) then
                loc%status = INVALID_POINT  ! point not found
            end if

        end do ! inewton




    end function computational_point
    !*****************************************************************************************










    !>  Project a function to the solution basis. Return modal coefficients.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/25/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function project(self,fcn) result(fmodes)
        class(element_t),       intent(in)      :: self
        class(function_t),      intent(inout)   :: fcn

        real(rk),       allocatable     :: fmodes(:)

        type(point_t),  allocatable     :: pts(:)
        real(rk),       allocatable     :: fvals(:), temp(:)
        real(rk)                        :: time


        !
        ! Call function for evaluation at quadrature nodes and multiply by quadrature weights
        !
        time  = 0._rk
        fvals = fcn%compute(time,self%quad_pts) * self%gq%vol%weights * self%jinv


        !
        ! Project
        !
        temp = matmul(transpose(self%gq%vol%val),fvals)
        fmodes = matmul(self%invmass,temp)


    end function project
    !*****************************************************************************************






    !>  Given a set of node indices, determine and return the associated face index
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/20/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_face_from_corners(self,corner_indices) result(face_index)
        class(element_t),   intent(in)  :: self
        integer(ik),        intent(in)  :: corner_indices(:)

        integer(ik), dimension(size(corner_indices))   :: corner_position

        integer(ik), allocatable    :: element_indices(:)

        integer(ik)                 :: nterms_1d, face_index, cindex, eindex, iface_test
        logical                     :: node_matches, face_match, &
                                       corner_one_in_face, corner_two_in_face, corner_three_in_face, corner_four_in_face

        !
        ! Get nodes from connectivity
        !
        element_indices = self%connectivity%get_element_nodes()


        do cindex = 1,size(corner_indices)
            do eindex = 1,size(element_indices)


                node_matches = (corner_indices(cindex) == element_indices(eindex))

                if (node_matches) then
                    corner_position(cindex) = eindex
                    exit
                end if


            end do
        end do


        !
        ! Determine element mapping index
        !
        nterms_1d = 0
        do while (nterms_1d*nterms_1d*nterms_1d < self%nterms_c)
            nterms_1d = nterms_1d + 1
        end do


        !
        ! Test corner positions against known face configurations
        !
        do iface_test = 1,NFACES
            corner_one_in_face   = any(face_corners(iface_test,:,nterms_1d - 1) == corner_position(1))
            corner_two_in_face   = any(face_corners(iface_test,:,nterms_1d - 1) == corner_position(2))
            corner_three_in_face = any(face_corners(iface_test,:,nterms_1d - 1) == corner_position(3))
            corner_four_in_face  = any(face_corners(iface_test,:,nterms_1d - 1) == corner_position(4))

            face_match = (corner_one_in_face .and. corner_two_in_face .and. corner_three_in_face .and. corner_four_in_face )

            if (face_match) then
                face_index = iface_test
                exit
            end if

        end do


        if (.not. face_match) call chidg_signal(FATAL,"element%get_face_from_corners: couldn't find a face index that matched the provided corners")


    end function get_face_from_corners
    !******************************************************************************************

















    subroutine destructor(self)
        type(element_t), intent(inout) :: self


    end subroutine

end module type_element
