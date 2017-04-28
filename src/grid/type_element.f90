module type_element
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,XI_MIN,XI_MAX,ETA_MIN, &
                                      ETA_MAX,ZETA_MIN,ZETA_MAX,ONE,ZERO, &
                                      DIR_1, DIR_2, DIR_3, DIR_THETA, XI_DIR, ETA_DIR, ZETA_DIR, &
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
    use ieee_arithmetic,            only: ieee_is_nan
    use DNAD_D
    implicit none




    !>  Element data type
    !!
    !!  ************************************************************************************
    !!  NOTE: could be dangerous to declare static arrays of elements using gfortran because
    !!        the compiler doesn't have complete finalization rules implemented. Using 
    !!        allocatables seems to work fine.
    !!  ************************************************************************************
    !!
    !!  Coordinate systems:
    !!      Coordinates could be in either 'Cartesian' or 'Cylindrical' systems.
    !!      As such, coordinate indices are marked by (1,2,3):
    !! 
    !!      'Cartesian'   system: 1 = x      ;  2 = y   ;  3 = z
    !!      'Cylindrical' system: 1 = theta  ;  2 = r   ;  3 = z
    !!
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
        integer(ik)     :: idomain_g                        ! Global index of parent domain
        integer(ik)     :: idomain_l                        ! Proc-local index of parent domain
        integer(ik)     :: ielement_g                       ! Domain-global index of element
        integer(ik)     :: ielement_l                       ! Proc-local index of the element

        integer(ik)     :: spacedim                         ! Number of spatial dimensions for the element
        integer(ik)     :: neqns                            ! Number of equations being solved
        integer(ik)     :: nterms_s                         ! Number of terms in solution expansion.  
        integer(ik)     :: nterms_c                         ! Number of terms in coordinate expansion. 
        integer(ik)     :: nterms_c_1d                      ! N-terms in 1d coordinate expansion.
        integer(ik)     :: ntime                            ! Number of time levels in solution

        ! Element quadrature points, mesh points and modes
        type(element_connectivity_t)    :: connectivity         ! Integer indices of the associated nodes in block node list
        type(point_t),  allocatable     :: quad_pts(:)          ! Coordinates of discrete quadrature points
        type(point_t),  allocatable     :: elem_pts(:)          ! Coordinates of discrete points defining element
        type(densevector_t)             :: coords               ! Modal expansion of coordinates (nterms_var,(x,y,z))
        character(:),   allocatable     :: coordinate_system    ! 'Cartesian', 'Cylindrical'

        ! Element metric terms
        real(rk), allocatable           :: metric(:,:,:)        ! metric matrix for each quadrature node    (mat_i,mat_j,quad_pt)
        real(rk), allocatable           :: jinv(:)              ! volume jacobian at quadrature nodes

        ! Matrices of cartesian gradients of basis/test functions
        real(rk), allocatable           :: grad1(:,:)           ! Grad of basis functions in at quadrature nodes
        real(rk), allocatable           :: grad2(:,:)           ! Grad of basis functions in at quadrature nodes
        real(rk), allocatable           :: grad3(:,:)           ! Grad of basis functions in at quadrature nodes
        real(rk), allocatable           :: grad1_trans(:,:)     ! transpose grad1
        real(rk), allocatable           :: grad2_trans(:,:)     ! transpose grad2
        real(rk), allocatable           :: grad3_trans(:,:)     ! transpose grad3

        ! Quadrature matrices
        type(quadrature_t), pointer     :: gq     => null()     ! Pointer to instance for solution expansion
        type(quadrature_t), pointer     :: gqmesh => null()     ! Pointer to instance for coordinate expansion

        ! Element-local mass, inverse mass matrices
        real(rk), allocatable           :: mass(:,:)        
        real(rk), allocatable           :: invmass(:,:)

        ! Element volume, approx. size of bounding box
        real(rk)                        :: vol
        real(rk)                        :: h(3)     


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


        ! Compute discrete value for at a given xi,eta,zeta.
        procedure, public   :: x                      
        procedure, public   :: y                      
        procedure, public   :: z                      
        procedure, public   :: grid_point           
        procedure, public   :: physical_point
        procedure, public   :: computational_point
        procedure, public   :: metric_point 
        procedure, public   :: solution_point   
        procedure, public   :: derivative_point


        ! Compute a projection of a function onto the solution basis
        procedure, public   :: project


        ! Get connected face
        procedure, public   :: get_face_from_corners

        ! Private utility procedure
        procedure           :: compute_element_matrices
        procedure           :: compute_mass_matrix
        procedure           :: compute_quadrature_gradients
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
    subroutine init_geom(self,spacedim,nodes,connectivity,idomain_l,ielem_l,coord_system)
        class(element_t),               intent(inout)   :: self
        integer(ik),                    intent(in)      :: spacedim
        type(point_t),                  intent(in)      :: nodes(:)
        type(element_connectivity_t),   intent(in)      :: connectivity
        integer(ik),                    intent(in)      :: idomain_l
        integer(ik),                    intent(in)      :: ielem_l
        character(*),                   intent(in)      :: coord_system

        character(:),   allocatable :: user_msg
        type(point_t),  allocatable :: points(:)
        real(rk),       allocatable :: element_mapping(:,:)
        real(rk),       allocatable :: modes1(:), modes2(:), modes3(:)
        real(rk)                    :: xmin, xmax, xwidth,  &
                                       ymin, ymax, ywidth,  &
                                       zmin, zmax, zwidth
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
        npts_1d          = mapping+1
        npts             = npts_1d * npts_1d * npts_1d
        self%nterms_c_1d = npts_1d
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
        ! Compute modal expansion of element coordinates
        !
        modes1 = matmul(element_mapping,self%elem_pts(:)%c1_)
        modes2 = matmul(element_mapping,self%elem_pts(:)%c2_)
        modes3 = matmul(element_mapping,self%elem_pts(:)%c3_)

        call self%coords%setvar(1,itime = 1,vals = modes1)
        call self%coords%setvar(2,itime = 1,vals = modes2)
        call self%coords%setvar(3,itime = 1,vals = modes3)



        !
        ! Compute approximate size of bounding box
        !
        xmax = maxval(points(:)%c1_)
        xmin = minval(points(:)%c1_)
        xwidth = abs(xmax - xmin)

        ymax = maxval(points(:)%c2_)
        ymin = minval(points(:)%c2_)
        ywidth = abs(ymax - ymin)

        zmax = maxval(points(:)%c3_)
        zmin = minval(points(:)%c3_)
        zwidth = abs(zmax - zmin)

        self%h(1) = xwidth
        self%h(2) = ywidth
        self%h(3) = zwidth



        !
        ! Set coordinate system and confirm initialization 
        !
        self%coordinate_system = coord_system
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
    !-----------------------------------------------------------------------------------------
    subroutine init_sol(self,neqns,nterms_s,ntime)
        class(element_t),   intent(inout) :: self
        integer(ik),        intent(in)    :: neqns
        integer(ik),        intent(in)    :: nterms_s
        integer(ik),        intent(in)    :: ntime

        integer(ik) :: ierr
        integer(ik) :: nnodes


        self%nterms_s    = nterms_s     ! number of terms in solution expansion
        self%neqns       = neqns        ! number of equations being solved
        self%ntime       = ntime        ! number of time steps in solution


        !
        ! With nterms_s and nterms_c defined:
        !   - assign quadrature instance
        !   - get number of quadrature nodes
        !
        call self%assign_quadrature()
        nnodes = self%gq%vol%nnodes


        !
        ! (Re)Allocate storage for element data structures
        !
        if (allocated(self%jinv)) &
            deallocate( self%jinv, self%metric, self%quad_pts,                &
                        self%grad1, self%grad2, self%grad3,                   &
                        self%grad1_trans, self%grad2_trans, self%grad3_trans, &
                        self%mass, self%invmass, self%dtau)

        allocate(self%jinv(nnodes),                 &
                 self%metric(3,3,nnodes),           &
                 self%quad_pts(nnodes),             &
                 self%grad1(nnodes,nterms_s),         &
                 self%grad2(nnodes,nterms_s),         &
                 self%grad3(nnodes,nterms_s),         &
                 self%grad1_trans(nterms_s,nnodes),   &
                 self%grad2_trans(nterms_s,nnodes),   &
                 self%grad3_trans(nterms_s,nnodes),   &
                 self%mass(nterms_s,nterms_s),      &
                 self%invmass(nterms_s,nterms_s),   &
                 self%dtau(neqns), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Call element metric and matrix calculation routines
        !
        call self%compute_element_matrices()



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









    !>  Subroutine computes element-specific matrices
    !!      - Mass matrix   (mass, invmass)
    !!      - Matrices of gradients of basis/test functions (grad1, grad2, grad3)
    !!      - Coordinates of quadrature points (quad_pts)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_element_matrices(self)
        class(element_t),   intent(inout)   :: self

        !
        ! Call to compute coordinates at each quadrature node
        !
        call self%compute_quadrature_coords()

        !
        ! Compute quadrature metrics
        !
        call self%compute_quadrature_metrics()

        !
        ! Call to compute mass matrix
        !
        call self%compute_mass_matrix()

        !
        ! Call to compute matrices of gradients at each quadrature node
        !
        call self%compute_quadrature_gradients()


    end subroutine compute_element_matrices
    !******************************************************************************************









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

        integer(ik)                 :: inode
        integer(ik)                 :: nnodes
        character(:),   allocatable :: coordinate_system

        real(rk),   dimension(self%gq%vol%nnodes)   ::  &
            d1dxi, d1deta, d1dzeta,                     &
            d2dxi, d2deta, d2dzeta,                     &
            d3dxi, d3deta, d3dzeta,                     &
            scaling_12, scaling_13, scaling_23, scaling_123


        nnodes = self%gq%vol%nnodes

        !
        ! Compute element metric terms
        !
        d1dxi   = matmul(self%gqmesh%vol%ddxi,  self%coords%getvar(1,itime = 1))
        d1deta  = matmul(self%gqmesh%vol%ddeta, self%coords%getvar(1,itime = 1))
        d1dzeta = matmul(self%gqmesh%vol%ddzeta,self%coords%getvar(1,itime = 1))

        d2dxi   = matmul(self%gqmesh%vol%ddxi,  self%coords%getvar(2,itime = 1))
        d2deta  = matmul(self%gqmesh%vol%ddeta, self%coords%getvar(2,itime = 1))
        d2dzeta = matmul(self%gqmesh%vol%ddzeta,self%coords%getvar(2,itime = 1))

        d3dxi   = matmul(self%gqmesh%vol%ddxi,  self%coords%getvar(3,itime = 1))
        d3deta  = matmul(self%gqmesh%vol%ddeta, self%coords%getvar(3,itime = 1))
        d3dzeta = matmul(self%gqmesh%vol%ddzeta,self%coords%getvar(3,itime = 1))




        !
        ! Define area/volume scaling for coordinate system
        !   Cartesian:
        !       12 = x-y  ;  13 = x-z  ;  23 = y-z
        !
        !   Cylindrical
        !       12 = r-theta  ;  13 = r-z      ;  23 = theta-z
        !
        select case (self%coordinate_system)
            case ('Cartesian')
                scaling_12  = ONE
                scaling_13  = ONE
                scaling_23  = ONE
                scaling_123 = ONE
            case ('Cylindrical')
                scaling_12  = self%quad_pts(:)%c1_
                scaling_13  = ONE
                scaling_23  = self%quad_pts(:)%c1_
                scaling_123 = self%quad_pts(:)%c1_
            case default
                call chidg_signal(FATAL,"element%compute_quadrature_metrics: Invalid coordinate system. Choose 'Cartesian' or 'Cylindrical'.")
        end select


        !
        ! Compute inverse cell mapping jacobian
        !
        self%jinv = scaling_123*(d1dxi*d2deta*d3dzeta  -  d1deta*d2dxi*d3dzeta - &
                                 d1dxi*d2dzeta*d3deta  +  d1dzeta*d2dxi*d3deta + &
                                 d1deta*d2dzeta*d3dxi  -  d1dzeta*d2deta*d3dxi)

        !
        ! Check for negative jacobians
        !
        if (any(self%jinv < ZERO)) call chidg_signal(FATAL,"element%compute_quadrature_metrics: Negative element jacobians. Check element quality and orientation.")


        !
        ! Compute element volume
        !
        self%vol = abs(sum(self%jinv * self%gq%vol%weights))


        !
        ! Loop through quadrature nodes and compute metric terms. This is the explicit formula
        ! for inverting a 3x3 matrix.
        !
        !   See: http://mathworld.wolfram.com/MatrixInverse.html 
        !
        do inode = 1,nnodes
            self%metric(1,1,inode) = ONE/self%jinv(inode) * scaling_23(inode) * (d2deta(inode)*d3dzeta(inode) - d2dzeta(inode)*d3deta(inode))
            self%metric(2,1,inode) = ONE/self%jinv(inode) * scaling_23(inode) * (d2dzeta(inode)*d3dxi(inode)  - d2dxi(inode)*d3dzeta(inode) )
            self%metric(3,1,inode) = ONE/self%jinv(inode) * scaling_23(inode) * (d2dxi(inode)*d3deta(inode)   - d2deta(inode)*d3dxi(inode)  )

            self%metric(1,2,inode) = ONE/self%jinv(inode) * scaling_13(inode) * (d1dzeta(inode)*d3deta(inode) - d1deta(inode)*d3dzeta(inode))
            self%metric(2,2,inode) = ONE/self%jinv(inode) * scaling_13(inode) * (d1dxi(inode)*d3dzeta(inode)  - d1dzeta(inode)*d3dxi(inode) )
            self%metric(3,2,inode) = ONE/self%jinv(inode) * scaling_13(inode) * (d1deta(inode)*d3dxi(inode)   - d1dxi(inode)*d3deta(inode)  )

            self%metric(1,3,inode) = ONE/self%jinv(inode) * scaling_12(inode) * (d1deta(inode)*d2dzeta(inode) - d1dzeta(inode)*d2deta(inode))
            self%metric(2,3,inode) = ONE/self%jinv(inode) * scaling_12(inode) * (d1dzeta(inode)*d2dxi(inode)  - d1dxi(inode)*d2dzeta(inode) )
            self%metric(3,3,inode) = ONE/self%jinv(inode) * scaling_12(inode) * (d1dxi(inode)*d2deta(inode)   - d1deta(inode)*d2dxi(inode)  )
        end do




    end subroutine compute_quadrature_metrics
    !******************************************************************************************











    !>  Compute matrices containing gradients of basis/test function
    !!  at each quadrature node.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_quadrature_gradients(self)
        class(element_t),   intent(inout)   :: self

        character(:),   allocatable :: user_msg
        integer(ik)                 :: iterm,inode

        do iterm = 1,self%nterms_s
            do inode = 1,self%gq%vol%nnodes
                !self%grad1(inode,iterm) = self%metric(1,1,inode) * self%gq%vol%ddxi(inode,iterm)   * (ONE/self%jinv(inode)) + &
                !                          self%metric(2,1,inode) * self%gq%vol%ddeta(inode,iterm)  * (ONE/self%jinv(inode)) + &
                !                          self%metric(3,1,inode) * self%gq%vol%ddzeta(inode,iterm) * (ONE/self%jinv(inode))

                !self%grad2(inode,iterm) = self%metric(1,2,inode) * self%gq%vol%ddxi(inode,iterm)   * (ONE/self%jinv(inode)) + &
                !                          self%metric(2,2,inode) * self%gq%vol%ddeta(inode,iterm)  * (ONE/self%jinv(inode)) + &
                !                          self%metric(3,2,inode) * self%gq%vol%ddzeta(inode,iterm) * (ONE/self%jinv(inode))

                !self%grad3(inode,iterm) = self%metric(1,3,inode) * self%gq%vol%ddxi(inode,iterm)   * (ONE/self%jinv(inode)) + &
                !                          self%metric(2,3,inode) * self%gq%vol%ddeta(inode,iterm)  * (ONE/self%jinv(inode)) + &
                !                          self%metric(3,3,inode) * self%gq%vol%ddzeta(inode,iterm) * (ONE/self%jinv(inode))

                self%grad1(inode,iterm) = self%metric(1,1,inode) * self%gq%vol%ddxi(inode,iterm)   + &
                                          self%metric(2,1,inode) * self%gq%vol%ddeta(inode,iterm)  + &
                                          self%metric(3,1,inode) * self%gq%vol%ddzeta(inode,iterm) 

                self%grad2(inode,iterm) = self%metric(1,2,inode) * self%gq%vol%ddxi(inode,iterm)   + &
                                          self%metric(2,2,inode) * self%gq%vol%ddeta(inode,iterm)  + &
                                          self%metric(3,2,inode) * self%gq%vol%ddzeta(inode,iterm) 

                self%grad3(inode,iterm) = self%metric(1,3,inode) * self%gq%vol%ddxi(inode,iterm)   + &
                                          self%metric(2,3,inode) * self%gq%vol%ddeta(inode,iterm)  + &
                                          self%metric(3,3,inode) * self%gq%vol%ddzeta(inode,iterm) 
            end do
        end do

        !
        ! Check for acceptable element
        !
        if (any(ieee_is_nan(self%grad1)) .or. &
            any(ieee_is_nan(self%grad2)) .or. &
            any(ieee_is_nan(self%grad3)) ) then
            user_msg = "element%compute_quadrature_gradients: Element failed to produce valid gradient information. &
                        Element quality is likely not reasonable."
            call chidg_signal(FATAL,"element%compute_quadrature_gradients: Element failed to produce valid gradient information. Element quality is likely not reasonable.")
        end if

        self%grad1_trans = transpose(self%grad1)
        self%grad2_trans = transpose(self%grad2)
        self%grad3_trans = transpose(self%grad3)

    end subroutine compute_quadrature_gradients
    !******************************************************************************************











    !>  Compute coordinates at each quadrature point
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

        integer(ik)                                 :: nnodes
        real(rk),   dimension(self%gq%vol%nnodes)   :: coord1, coord2, coord3
        integer(ik)                                 :: inode

        nnodes = self%gq%vol%nnodes

        !
        ! compute coordinates associated with quadrature points
        !
        coord1 = matmul(self%gqmesh%vol%val,self%coords%getvar(1,itime = 1))
        coord2 = matmul(self%gqmesh%vol%val,self%coords%getvar(2,itime = 1))
        coord3 = matmul(self%gqmesh%vol%val,self%coords%getvar(3,itime = 1))


        !
        ! Initialize each point with coordinates
        !
        do inode = 1,nnodes
            call self%quad_pts(inode)%set(coord1(inode),coord2(inode),coord3(inode))
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

        !print*, 'jinv: ', self%jinv


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




    !>  Convert local(xi,eta,zeta) coordinates to global coordinates(x,y,z)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !------------------------------------------------------------------------------------------
    function physical_point(self,xi,eta,zeta) result(phys_point)
        class(element_t),   intent(in)  :: self
        real(rk),           intent(in)  :: xi,eta,zeta

        real(rk)                   :: val1, val2, val3
        type(point_t)              :: point, phys_point
        real(rk)                   :: polyvals(self%nterms_c)
        integer(ik)                :: iterm, spacedim

        call point%set(xi,eta,zeta)

        spacedim = self%spacedim

        !
        ! Evaluate polynomial modes at node location
        !
        do iterm = 1,self%nterms_c

            polyvals(iterm)  = polynomialVal(spacedim,self%nterms_c,iterm,point)

        end do

        
        !
        ! Evaluate x from dot product of modes and polynomial values
        !
        val1 = dot_product(self%coords%getvar(1, itime=1),polyvals)
        val2 = dot_product(self%coords%getvar(2, itime=1),polyvals)
        val3 = dot_product(self%coords%getvar(3, itime=1),polyvals)


        !
        ! Set physical coordinates
        !
        call phys_point%set(val1,val2,val3) 


    end function physical_point
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
    !!  @param[in]  phys_dir    physical coordinate being differentiated
    !!  @param[in]  comp_dir    Computational coordinate being differentiated with respect to
    !!  @param[in]  xi          Computational coordinate - xi
    !!  @param[in]  eta         Computational coordinate - eta
    !!  @param[in]  zeta        Computational coordinate - zeta
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !-----------------------------------------------------------------------------------------
    function metric_point(self,phys_dir,comp_dir,xi,eta,zeta,scale) result(val)
        class(element_t),   intent(in)              :: self
        integer(ik),        intent(in)              :: phys_dir
        integer(ik),        intent(in)              :: comp_dir
        real(rk),           intent(in)              :: xi, eta, zeta
        logical,            intent(in), optional    :: scale
        
        real(rk)        :: val, r
        type(point_t)   :: node
        real(rk)        :: polyvals(self%nterms_c)
        integer(ik)     :: iterm, spacedim


        if (phys_dir > 3) call chidg_signal(FATAL,"element%metric_point: phys_dir exceeded 3 physical coordinates")
        if (comp_dir > 3) call chidg_signal(FATAL,"element%metric_point: comp_dir exceeded 3 physical coordinates")

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
        val = dot_product(self%coords%getvar(phys_dir, itime=1), polyvals)



        !
        ! Apply scaling due to coordinate system.
        !
        if (present(scale)) then
            if (scale) then
                if (self%coordinate_system == 'Cartesian') then

                else if (self%coordinate_system == 'Cylindrical') then
                    if (phys_dir == DIR_THETA) then
                        r = self%grid_point(1,xi,eta,zeta)
                        val = val * r
                    end if
                end if
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
    !-----------------------------------------------------------------------------------------
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
        metric(1,1) = self%metric_point(DIR_1,XI_DIR,  xi,eta,zeta)
        metric(2,1) = self%metric_point(DIR_2,XI_DIR,  xi,eta,zeta)
        metric(3,1) = self%metric_point(DIR_3,XI_DIR,  xi,eta,zeta)
        metric(1,2) = self%metric_point(DIR_1,ETA_DIR, xi,eta,zeta)
        metric(2,2) = self%metric_point(DIR_2,ETA_DIR, xi,eta,zeta)
        metric(3,2) = self%metric_point(DIR_3,ETA_DIR, xi,eta,zeta)
        metric(1,3) = self%metric_point(DIR_1,ZETA_DIR,xi,eta,zeta)
        metric(2,3) = self%metric_point(DIR_2,ZETA_DIR,xi,eta,zeta)
        metric(3,3) = self%metric_point(DIR_3,ZETA_DIR,xi,eta,zeta)


        !
        ! Compute inverse cell mapping jacobian
        !
        jinv = metric(1,1)*metric(2,2)*metric(3,3) - metric(1,2)*metric(2,1)*metric(3,3) - &
               metric(1,1)*metric(2,3)*metric(3,2) + metric(1,3)*metric(2,1)*metric(3,2) + &
               metric(1,2)*metric(2,3)*metric(3,1) - metric(1,3)*metric(2,2)*metric(3,1)





        do iterm = 1,self%nterms_s
            if (dir == DIR_1) then
                dxi_dx   = metric(2,2)*metric(3,3) - metric(2,3)*metric(3,2)
                deta_dx  = metric(2,3)*metric(3,1) - metric(2,1)*metric(3,3)
                dzeta_dx = metric(2,1)*metric(3,2) - metric(2,2)*metric(3,1)
                deriv(iterm) = dxi_dx   * ddxi(iterm)   * (ONE/jinv) + &
                               deta_dx  * ddeta(iterm)  * (ONE/jinv) + &
                               dzeta_dx * ddzeta(iterm) * (ONE/jinv)
            else if (dir == DIR_2) then
                dxi_dy   = metric(1,3)*metric(3,2) - metric(1,2)*metric(3,3)
                deta_dy  = metric(1,1)*metric(3,3) - metric(1,3)*metric(3,1)
                dzeta_dy = metric(1,2)*metric(3,1) - metric(1,1)*metric(3,2)
                deriv(iterm) = dxi_dy   * ddxi(iterm)   * (ONE/jinv) + &
                               deta_dy  * ddeta(iterm)  * (ONE/jinv) + &
                               dzeta_dy * ddzeta(iterm) * (ONE/jinv)
            else if (dir == DIR_3) then
                dxi_dz   = metric(1,2)*metric(2,3) - metric(1,3)*metric(2,2)
                deta_dz  = metric(1,3)*metric(2,1) - metric(1,1)*metric(2,3)
                dzeta_dz = metric(1,1)*metric(2,2) - metric(1,2)*metric(2,1)
                deriv(iterm) = dxi_dz   * ddxi(iterm)   * (ONE/jinv) + &
                               deta_dz  * ddeta(iterm)  * (ONE/jinv) + &
                               dzeta_dz * ddzeta(iterm) * (ONE/jinv)
            else
                call chidg_signal(FATAL,"element%derivative_point: Invalid value for 'dir' parameter. (1,2,3).")
            end if
        end do



        !
        ! Evaluate x from dot product of modes and polynomial values
        !
        val = dot_product(q%getvar(ivar,itime),deriv)

    end function derivative_point
    !*****************************************************************************************







    
    !>  Compute a computational location(xi,eta,zeta), based on the location in 
    !!  physical space (x,y,z), (r,theta,z)
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
    function computational_point(self,coord1,coord2,coord3) result(loc)
        class(element_t),   intent(in)  :: self
        real(rk),           intent(in)  :: coord1
        real(rk),           intent(in)  :: coord2
        real(rk),           intent(in)  :: coord3

        type(point_t)       :: loc, point_n
        real(rk)            :: xi, eta, zeta
        integer(ik)         :: inewton

        real(rk)    :: mat(3,3), minv(3,3)
        real(rk)    :: R(3)
        real(rk)    :: dcoord(3)
        real(rk)    :: res, tol


        !tol = 2._rk*RKTOL
        !tol = 10._rk*RKTOL
        !tol = RKTOL

        tol = 1000._rk*RKTOL

        !
        ! Newton iteration to find the donor local coordinates
        !
        xi   = ZERO
        eta  = ZERO
        zeta = ZERO
        do inewton = 1,20


            !
            ! Compute local physical coordinates as a function of xi,eta,zeta
            !
            point_n  = self%physical_point(xi,eta,zeta)


            !
            ! Assemble residual vector
            !
            R(1) = -(point_n%c1_ - coord1)
            R(2) = -(point_n%c2_ - coord2)
            R(3) = -(point_n%c3_ - coord3)


            !
            ! Assemble coordinate jacobian matrix
            !
            mat(1,1) = self%metric_point(DIR_1,XI_DIR,  xi,eta,zeta)
            mat(2,1) = self%metric_point(DIR_2,XI_DIR,  xi,eta,zeta)
            mat(3,1) = self%metric_point(DIR_3,XI_DIR,  xi,eta,zeta)
            mat(1,2) = self%metric_point(DIR_1,ETA_DIR, xi,eta,zeta)
            mat(2,2) = self%metric_point(DIR_2,ETA_DIR, xi,eta,zeta)
            mat(3,2) = self%metric_point(DIR_3,ETA_DIR, xi,eta,zeta)
            mat(1,3) = self%metric_point(DIR_1,ZETA_DIR,xi,eta,zeta)
            mat(2,3) = self%metric_point(DIR_2,ZETA_DIR,xi,eta,zeta)
            mat(3,3) = self%metric_point(DIR_3,ZETA_DIR,xi,eta,zeta)


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

        character(:),   allocatable :: user_msg
        integer(ik),    allocatable :: element_indices(:), face_indices(:)
        integer(ik)                 :: face_index, cindex, eindex, iface_test
        logical                     :: node_matches, face_match, &
                                       corner_one_in_face, corner_two_in_face, &
                                       corner_three_in_face, corner_four_in_face

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
        ! Test corner positions against known face configurations 
        ! to determine face index:
        !
        do iface_test = 1,NFACES

            face_indices = face_corners(iface_test,:,self%nterms_c_1d - 1)
            corner_one_in_face   = any(face_indices == corner_position(1))
            corner_two_in_face   = any(face_indices == corner_position(2))
            corner_three_in_face = any(face_indices == corner_position(3))
            corner_four_in_face  = any(face_indices == corner_position(4))

            face_match = (corner_one_in_face   .and. &
                          corner_two_in_face   .and. &
                          corner_three_in_face .and. &
                          corner_four_in_face )

            if (face_match) then
                face_index = iface_test
                exit
            end if

        end do


        user_msg = "element%get_face_from_cornders: Couldn't find a face index that matched &
                    the provided corner indices."
        if (.not. face_match) call chidg_signal(FATAL,user_msg)


    end function get_face_from_corners
    !******************************************************************************************

















    subroutine destructor(self)
        type(element_t), intent(inout) :: self


    end subroutine

end module type_element
