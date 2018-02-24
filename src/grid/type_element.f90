module type_element
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,NEDGES,XI_MIN,XI_MAX,ETA_MIN, &
                                      ETA_MAX,ZETA_MIN,ZETA_MAX,ONE,ZERO,THIRD, &
                                      DIR_1, DIR_2, DIR_3, DIR_THETA, XI_DIR,   &
                                      ETA_DIR, ZETA_DIR, TWO_DIM, THREE_DIM,    &
                                      RKTOL, VALID_POINT, INVALID_POINT, NO_PMM_ASSIGNED, &
                                      ZERO, TWO, CARTESIAN, CYLINDRICAL, DIR_R, NO_ID
    use mod_grid,               only: get_element_mapping, face_corners
    use mod_reference_elements, only: get_reference_element, ref_elems
    use mod_polynomial,         only: polynomial_val, dpolynomial_val
    use mod_inv,                only: inv, inv_3x3
    use mod_determinant,        only: det_3x3
    use mod_io,                 only: gq_rule


    use type_point,                 only: point_t
    use type_densevector,           only: densevector_t
    use type_function,              only: function_t
    use type_element_connectivity,  only: element_connectivity_t
    use type_reference_element,     only: reference_element_t
    use DNAD_D
    use ieee_arithmetic,            only: ieee_value, ieee_quiet_nan, ieee_is_nan
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
    !----------------------------------------------------------------------------------
    type, public :: element_t

        ! Element location
        integer(ik)                 :: element_location(5)  ! [idomain_g, idomain_l, ielement_g, ielement_l, iproc], useful for nonblocking send
        integer(ik)                 :: idomain_g            ! Global index of parent domain
        integer(ik)                 :: idomain_l            ! Proc-local index of parent domain
        integer(ik)                 :: ielement_g           ! Domain-global index of element
        integer(ik)                 :: ielement_l           ! Proc-local index of the element
        integer(ik)                 :: iproc                ! Processor the element is associated with.

        ! Element data
        integer(ik)                 :: element_data(8)      ! [element_type, spacedim, coordinate_system, neqns, nterms_s, nterms_c, ntime, interpolation_level]
        integer(ik)                 :: element_type         ! 1=linear, 2=quadratic, 3=cubic, 4=quartic, etc.
        integer(ik)                 :: spacedim             ! Number of spatial dimensions for the element.
        integer(ik)                 :: coordinate_system    ! CARTESIAN, CYLINDRICAL. parameters from mod_constants.
        integer(ik)                 :: eqn_ID = NO_ID       ! Equation set identifier the element is associated with.
        integer(ik)                 :: neqns                ! Number of equations being solved.
        integer(ik)                 :: nterms_s             ! Number of terms in solution expansion.  
        integer(ik)                 :: nterms_c             ! Number of terms in coordinate expansion. 
        integer(ik)                 :: ntime                ! Number of time levels in solution.
        integer(ik)                 :: interpolation_level  ! 1=lowest, 2-> are higher.
        integer(ik)                 :: recv_comm    = NO_ID ! chidg_vector access if element is initialized on another processor
        integer(ik)                 :: recv_domain  = NO_ID ! chidg_vector access if element is initialized on another processor
        integer(ik)                 :: recv_element = NO_ID ! chidg_vector access if element is initialized on another processor


        ! Connectivty and linear transformation martrix for 
        ! converting between nodal and modal representations
        integer(ik),    allocatable :: connectivity(:)      ! Integer indices of the associated nodes in block node list
        real(rk),       allocatable :: nodes_to_modes(:,:)  ! Transformation matrix for converting nodal values to modal coefficients


        ! Modal representations of element coordinates/velocity
        type(densevector_t)         :: coords               ! Modal representation of undeformed element coordinates (nterms_var,(x,y,z))
        type(densevector_t)         :: coords_def           ! Modal representation of deformed element coordinates (nterms, (x,y,z))
        type(densevector_t)         :: coords_vel           ! Modal representation of coordinate velocities (nterms,(v_x,v_y,v_z))


        ! Element geometry at support nodes
        real(rk),       allocatable :: node_coords(:,:)     ! Node coordinates, undeformed element, local element ordering
        real(rk),       allocatable :: node_coords_def(:,:) ! Node coordinates, deformed element,   local element ordering
        real(rk),       allocatable :: node_coords_vel(:,:) ! Node velocities,  deformed element,   local element ordering

        ! Element geometry at interpolation nodes
        real(rk),       allocatable :: interp_coords(:,:)       ! Undeformed coordinates at element interpolation nodes
        real(rk),       allocatable :: interp_coords_def(:,:)   ! Deformed coordinates at element interpolation nodes
        real(rk),       allocatable :: interp_coords_vel(:,:)   ! Coordinate velocities at element interpolation nodes
        real(rk),       allocatable :: metric(:,:,:)            ! inverted jacobian matrix for each volume node (mat_i,mat_j,volume_pt)
        real(rk),       allocatable :: edge_metric(:,:,:,:)     ! inverted jacobian matrix for each edge node (mat_i,mat_j,edge_pt)
        real(rk),       allocatable :: jinv(:)                  ! Differential volume ratio: Undeformed Volume/Reference Volume
        real(rk),       allocatable :: jinv_def(:)              ! Differential volume ratio: Deformed Volume/Reference Volume


        ! Arbitrary Lagrangian Eulerian data
        !   : This defines a mapping from some deformed element back to the original
        !   : undeformed element with the idea that the governing equations are transformed
        !   : and solved on the undeformed element.
        integer(ik)                 :: pmm_ID = NO_PMM_ASSIGNED
        real(rk),       allocatable :: ale_Dinv(:,:,:)          ! Deformation gradient: deformed element/undeformed element
        real(rk),       allocatable :: ale_g(:)                 ! Differential volume ratio: Deformed Volume/Undeformed Volume
        real(rk),       allocatable :: ale_g_grad1(:)
        real(rk),       allocatable :: ale_g_grad2(:)
        real(rk),       allocatable :: ale_g_grad3(:)
        real(rk),       allocatable :: ale_g_modes(:)           ! A modal expansion of the ale differential volume ratio in the solution basis


        ! Matrices of physical gradients of basis/test functions
        real(rk),       allocatable :: grad1(:,:)           ! Grad of basis functions in at quadrature nodes
        real(rk),       allocatable :: grad2(:,:)           ! Grad of basis functions in at quadrature nodes
        real(rk),       allocatable :: grad3(:,:)           ! Grad of basis functions in at quadrature nodes
        real(rk),       allocatable :: grad1_trans(:,:)     ! transpose grad1
        real(rk),       allocatable :: grad2_trans(:,:)     ! transpose grad2
        real(rk),       allocatable :: grad3_trans(:,:)     ! transpose grad3

        real(rk),       allocatable :: edge_grad1(:,:,:)
        real(rk),       allocatable :: edge_grad2(:,:,:)
        real(rk),       allocatable :: edge_grad3(:,:,:)

        ! Element-local mass, inverse mass matrices
        real(rk),       allocatable :: mass(:,:)        
        real(rk),       allocatable :: invmass(:,:)

        ! Element volume, approx. size of bounding box
        real(rk)                    :: vol
        real(rk)                    :: vol_ale
        real(rk)                    :: h(3)     
        real(rk),       allocatable :: dtau(:)              ! a pseudo-timestep for each equation. Used in the nonlinear solver.

        ! Reference element and interpolators
        type(reference_element_t),  pointer :: basis_s => null()  ! Pointer to solution basis and interpolator
        type(reference_element_t),  pointer :: basis_c => null()  ! Pointer to coordinate basis and interpolator

        ! Logical tests
        logical :: geom_initialized = .false.
        logical :: numInitialized   = .false.


    contains

        ! Initialization procedures
        procedure, public   :: init_geom
        procedure, public   :: init_sol
        procedure, public   :: init_eqn


        ! Undeformed element procedures
        procedure, public   :: update_interpolations
        procedure, private  :: interpolate_coords
        procedure, private  :: interpolate_metrics
        procedure, private  :: interpolate_gradients
        procedure, private  :: compute_mass_matrix

        ! Deformed element/ALE procedures
        procedure, public   :: set_displacements_velocities
        procedure, public   :: update_interpolations_ale
        procedure, private  :: interpolate_coords_ale
        procedure, private  :: interpolate_metrics_ale

        ! Edge procedures
        procedure, private  :: interpolate_metrics_edge
        procedure, private  :: interpolate_gradients_edge

        ! Compute discrete value for a given xi,eta,zeta.
        procedure, public   :: x                      
        procedure, public   :: y                      
        procedure, public   :: z                      
        procedure, public   :: grid_point           
        procedure, public   :: physical_point
        procedure, public   :: computational_point
        procedure, public   :: metric_point 
        procedure, public   :: solution_point   
        procedure, public   :: ale_point
        !procedure, public   :: derivative_point


        ! Compute a projection of a function onto the solution basis
        procedure, public   :: project

        ! Get connected face
        procedure, public   :: get_face_from_corners




        final               :: destructor

    end type element_t
    !*********************************************************************************

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
    !----------------------------------------------------------------------------------
    subroutine init_geom(self,nodes,connectivity,etype,location,coord_system)
        class(element_t),               intent(inout)   :: self
        real(rk),                       intent(in)      :: nodes(:,:)
        integer(ik),                    intent(in)      :: connectivity(:)
        integer(ik),                    intent(in)      :: etype
        integer(ik),                    intent(in)      :: location(5)
        character(*),                   intent(in)      :: coord_system

        character(:),   allocatable :: user_msg
        real(rk),       allocatable :: nodes_l(:,:), dnodes(:,:), vnodes(:,:)
        real(rk),       allocatable :: modes1(:), modes2(:), modes3(:)
        real(rk)                    :: xmin, xmax, xwidth,  &
                                       ymin, ymax, ywidth,  &
                                       zmin, zmax, zwidth
        integer(ik)                 :: ierr, ipt, npts_1d, npts, &
                                       mapping, inode, ref_ID_c
        integer(ik)                 :: ntime = 1


        user_msg = "element%init_geom: element already initialized."
        if (self%geom_initialized) call chidg_signal(FATAL,user_msg)


        !
        ! Get connectivity info
        !
        mapping               = etype
        self%idomain_g        = location(1)
        self%idomain_l        = location(2)
        self%ielement_g       = location(3)
        self%ielement_l       = location(4)
        self%iproc            = location(5)
        self%element_location = location
        self%connectivity     = connectivity



        !
        ! Get reference element (reference nodes only, no interpolators)
        !
        self%element_type = mapping
        ref_ID_c          = get_reference_element(element_type = mapping)
        self%basis_c => ref_elems(ref_ID_c)



        !
        ! Accumulate coordinates for current element from node list.
        !
        npts = self%basis_c%nnodes_r()
        allocate(nodes_l(npts,3), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Accumulate local nodes
        !
        do ipt = 1,npts
            ! Get node index
            inode = connectivity(ipt)

            ! Assemble local node list from global
            ! Default node coordinate delta's = zero
            nodes_l(ipt,:) = nodes(inode,:)
        end do !ipt





        !
        ! Get element mapping
        !
        self%spacedim = 3
        self%ntime    = 1
        self%nodes_to_modes = self%basis_c%nodes_to_modes
        self%nterms_c = size(self%nodes_to_modes,1)
        user_msg = "element%init_geom: mapping and points do not match."
        if (self%nterms_c /= size(nodes_l,1)) call chidg_signal(FATAL,user_msg)



        !
        ! Allocate storage
        !
        allocate(self%node_coords(self%nterms_c,3),stat=ierr)
        call self%coords%init(self%nterms_c,self%spacedim,self%ntime,self%idomain_g,self%idomain_l,self%ielement_g,self%ielement_l)
        self%node_coords = nodes_l

        
        !
        ! Compute modal expansion of element coordinates
        !
        modes1 = matmul(self%nodes_to_modes,self%node_coords(:,1))
        modes2 = matmul(self%nodes_to_modes,self%node_coords(:,2))
        modes3 = matmul(self%nodes_to_modes,self%node_coords(:,3))

        call self%coords%setvar(1,itime = 1,vals = modes1)
        call self%coords%setvar(2,itime = 1,vals = modes2)
        call self%coords%setvar(3,itime = 1,vals = modes3)



        !
        ! Compute approximate size of bounding box
        !
        xmax = maxval(nodes_l(:,1))
        xmin = minval(nodes_l(:,1))
        xwidth = abs(xmax - xmin)

        ymax = maxval(nodes_l(:,2))
        ymin = minval(nodes_l(:,2))
        ywidth = abs(ymax - ymin)

        zmax = maxval(nodes_l(:,3))
        zmin = minval(nodes_l(:,3))
        zwidth = abs(zmax - zmin)

        self%h(1) = xwidth
        self%h(2) = ywidth
        self%h(3) = zwidth



        !
        ! Set coordinate system and confirm initialization 
        !
        select case(trim(coord_system))
            case('Cartesian')
                self%coordinate_system = CARTESIAN
            case('Cylindrical')
                self%coordinate_system = CYLINDRICAL
            case default
                call chidg_signal_one(FATAL,"element%init_geom: Invalid coordinate system.",trim(coord_system))
        end select
        self%geom_initialized = .true.   


        !
        ! Store element_data(1-2)
        !
        self%element_data(1) = self%element_type
        self%element_data(2) = self%spacedim
        self%element_data(3) = self%coordinate_system


        !
        ! ALE initialization
        !   Default: zero displacements/velocities
        !
        dnodes = nodes
        vnodes = nodes
        dnodes = ZERO
        vnodes = ZERO
        call self%set_displacements_velocities(dnodes,vnodes)


    end subroutine init_geom
    !***********************************************************************************







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
    !----------------------------------------------------------------------------------
    subroutine init_sol(self,interpolation,level,nterms_s,nfields,ntime)
        class(element_t),   intent(inout) :: self
        character(*),       intent(in)    :: interpolation
        integer(ik),        intent(in)    :: level
        integer(ik),        intent(in)    :: nterms_s
        integer(ik),        intent(in)    :: nfields
        integer(ik),        intent(in)    :: ntime

        integer(ik) :: ierr
        integer(ik) :: nnodes, nnodes_edge
        integer(ik) :: ref_ID_s, ref_ID_c
        

        self%nterms_s    = nterms_s     ! number of terms in solution expansion
        self%neqns       = nfields      ! number of equations being solved
        self%ntime       = ntime        ! number of time steps in solution


        !
        ! With nterms_s and nterms_c defined:
        !   - assign quadrature instance
        !   - get number of quadrature nodes
        !
        ref_ID_s = get_reference_element(element_type = self%element_type,  &
                                         polynomial   = 'Legendre',         &
                                         nterms       = nterms_s,           &
                                         node_set     = interpolation,      &
                                         level        = level,              &
                                         nterms_rule  = nterms_s)

        ref_ID_c = get_reference_element(element_type = self%element_type,  &
                                         polynomial   = 'Legendre',         &
                                         nterms       = self%nterms_c,      &
                                         node_set     = interpolation,      &
                                         level        = level,              &
                                         nterms_rule  = nterms_s)
        self%basis_s => ref_elems(ref_ID_s)
        self%basis_c => ref_elems(ref_ID_c)



        !
        ! (Re)Allocate storage for element data structures
        !
        if (allocated(self%jinv))                       &
            deallocate( self%jinv,                      &
                        self%jinv_def,                  &
                        self%metric,                    &
                        self%edge_metric,               &
                        self%interp_coords,             &
                        self%interp_coords_def,         &
                        self%grad1,                     &
                        self%grad2,                     &
                        self%grad3,                     &
                        self%grad1_trans,               &
                        self%grad2_trans,               &
                        self%grad3_trans,               &
                        self%edge_grad1,                &
                        self%edge_grad2,                &
                        self%edge_grad3,                &
                        self%mass,                      &
                        self%invmass,                   &
                        self%interp_coords_vel,         &
                        self%ale_Dinv,                  &
                        self%ale_g,                     &
                        self%ale_g_grad1,               &
                        self%ale_g_grad2,               &
                        self%ale_g_grad3,               &
                        self%ale_g_modes,               &
                        self%dtau                       &
                        )
            

        nnodes = ref_elems(ref_ID_s)%nnodes_elem()
        nnodes_edge = ref_elems(ref_ID_s)%nnodes_edge()
        allocate(self%jinv(nnodes),                         &
                 self%jinv_def(nnodes),                     &
                 self%grad1(nnodes,nterms_s),               &
                 self%grad2(nnodes,nterms_s),               &
                 self%grad3(nnodes,nterms_s),               &
                 self%grad1_trans(nterms_s,nnodes),         &
                 self%grad2_trans(nterms_s,nnodes),         &
                 self%grad3_trans(nterms_s,nnodes),         &
                 self%edge_grad1(nnodes_edge,nterms_s,NEDGES),   &
                 self%edge_grad2(nnodes_edge,nterms_s,NEDGES),   &
                 self%edge_grad3(nnodes_edge,nterms_s,NEDGES),   &
                 self%mass(nterms_s,nterms_s),              &
                 self%invmass(nterms_s,nterms_s),           &
                 self%interp_coords(nnodes,3),              & 
                 self%interp_coords_def(nnodes,3),          &
                 self%interp_coords_vel(nnodes,3),          &
                 self%metric(3,3,nnodes),                   &
                 self%edge_metric(3,3,nnodes_edge,NEDGES),  &
                 self%ale_Dinv(3,3,nnodes),                 &
                 self%ale_g(nnodes),                        &
                 self%ale_g_grad1(nnodes),                  &
                 self%ale_g_grad2(nnodes),                  &
                 self%ale_g_grad3(nnodes),                  &
                 self%ale_g_modes(nterms_s),                &
                 self%dtau(nfields), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Call element metric and matrix calculation routines
        !
        call self%update_interpolations()
        call self%update_interpolations_ale()

        !
        ! Confirm element numerics were initialized
        !
        self%numInitialized = .true.    


        !
        ! Store element_data(4-8)
        !
        self%element_data(4) = self%neqns
        self%element_data(5) = self%nterms_s
        self%element_data(6) = self%nterms_c
        self%element_data(7) = self%ntime
        self%element_data(8) = level




    end subroutine init_sol
    !**********************************************************************************






    !>  Assign equation_set_t instance association.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2018
    !!
    !!  @param[in]  eqn_ID  Equation set index associated with the element.
    !!                      chidg%data%eqnset(eqn_ID)
    !!
    !--------------------------------------------------------------------------------
    subroutine init_eqn(self,eqn_ID)
        class(element_t),   intent(inout)   :: self
        integer(ik),        intent(in)      :: eqn_ID

        self%eqn_ID = eqn_ID

    end subroutine init_eqn
    !********************************************************************************








    !>  Subroutine computes element-specific matrices
    !!      - Coordinates of interpolation node set
    !!      - Metric terms at interpolation node set
    !!      - Mass matrix   (mass, invmass)
    !!      - Matrices of gradients of basis/test functions (grad1, grad2, grad3)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-----------------------------------------------------------------------------------
    subroutine update_interpolations(self)
        class(element_t),   intent(inout)   :: self

        !
        ! Call to compute coordinates at each interpolation node
        !
        call self%interpolate_coords()

        !
        ! Compute interpolation metrics
        !
        call self%interpolate_metrics()
        call self%interpolate_metrics_edge()

        !
        ! Call to compute mass matrix
        !
        call self%compute_mass_matrix()

        !
        ! Call to compute matrices of gradients at each interpolation node
        !
        call self%interpolate_gradients()
        call self%interpolate_gradients_edge()


    end subroutine update_interpolations
    !***********************************************************************************









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
    !----------------------------------------------------------------------------------
    subroutine interpolate_metrics(self)
        class(element_t),    intent(inout)   :: self

        integer(ik)                 :: inode, nnodes, ierr
        character(:),   allocatable :: coordinate_system, user_msg

        real(rk),   dimension(:),       allocatable :: scaling_row2, weights
        real(rk),   dimension(:,:),     allocatable :: val, ddxi, ddeta, ddzeta
        real(rk),   dimension(:,:,:),   allocatable :: jacobian

        nnodes  = self%basis_c%nnodes_elem()
        weights = self%basis_c%weights_element()
        val     = self%basis_c%interpolator_element('Value')
        ddxi    = self%basis_c%interpolator_element('ddxi')
        ddeta   = self%basis_c%interpolator_element('ddeta')
        ddzeta  = self%basis_c%interpolator_element('ddzeta')

        !
        ! Compute coordinate jacobian matrix at interpolation nodes
        !
        allocate(jacobian(3,3,nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError
        jacobian(1,1,:) = matmul(ddxi,   self%coords%getvar(1,itime = 1))
        jacobian(1,2,:) = matmul(ddeta,  self%coords%getvar(1,itime = 1))
        jacobian(1,3,:) = matmul(ddzeta, self%coords%getvar(1,itime = 1))

        jacobian(2,1,:) = matmul(ddxi,   self%coords%getvar(2,itime = 1))
        jacobian(2,2,:) = matmul(ddeta,  self%coords%getvar(2,itime = 1))
        jacobian(2,3,:) = matmul(ddzeta, self%coords%getvar(2,itime = 1))

        jacobian(3,1,:) = matmul(ddxi,   self%coords%getvar(3,itime = 1))
        jacobian(3,2,:) = matmul(ddeta,  self%coords%getvar(3,itime = 1))
        jacobian(3,3,:) = matmul(ddzeta, self%coords%getvar(3,itime = 1))


        !
        ! Add coordinate system scaling to jacobian matrix
        !
        allocate(scaling_row2(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError

        select case (self%coordinate_system)
            case (CARTESIAN)
                scaling_row2 = ONE
            case (CYLINDRICAL)
                scaling_row2 = self%interp_coords(:,1)
            case default
                user_msg = "element%interpolate_metrics: Invalid coordinate system."
                call chidg_signal(FATAL,user_msg)
        end select



        !
        ! Apply coorindate system scaling
        !
        jacobian(2,1,:) = jacobian(2,1,:)*scaling_row2
        jacobian(2,2,:) = jacobian(2,2,:)*scaling_row2
        jacobian(2,3,:) = jacobian(2,3,:)*scaling_row2




        !
        ! Compute inverse cell mapping jacobian
        !
        do inode = 1,nnodes
            self%jinv(inode) = det_3x3(jacobian(:,:,inode))
        end do


        !
        ! Check for negative jacobians
        !
        user_msg = "element%interpolate_metrics: Negative element &
                    volume detected. Check element quality and orientation."
        if (any(self%jinv < ZERO)) call chidg_signal(FATAL,user_msg)


        !
        ! Compute total volume by integrating differential volumes
        !
        self%vol = abs(sum(self%jinv * weights))


        !
        ! Invert jacobian matrix at each interpolation node
        !
        do inode = 1,nnodes
            self%metric(:,:,inode) = inv_3x3(jacobian(:,:,inode))
        end do



    end subroutine interpolate_metrics
    !************************************************************************************




    !> Compute edge metric and jacobian terms
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/6/2017
    !!
    !! TODO: Generalized 2D physical coordinates. Currently assumes x-y
    !!
    !------------------------------------------------------------------------------------
    subroutine interpolate_metrics_edge(self)
        class(element_t),    intent(inout)   :: self

        integer(ik)                 :: inode, iedge, nnodes_edge, ierr
        character(:),   allocatable :: coordinate_system, user_msg

        real(rk),   dimension(:),       allocatable :: scaling_row2
        real(rk),   dimension(:,:),     allocatable :: val, ddxi, ddeta, ddzeta
        real(rk),   dimension(:,:,:),   allocatable :: jacobian

        nnodes_edge  = self%basis_c%nnodes_edge()

        ! Element jacobian matrix
        allocate(jacobian(3,3,nnodes_edge), stat=ierr)
        if (ierr /= 0) call AllocationError
        
        ! Coordinate system scaling
        allocate(scaling_row2(nnodes_edge), stat=ierr)
        if (ierr /= 0) call AllocationError

        do iedge = 1,NEDGES

            val     = self%basis_c%interpolator_edge('Value',  iedge)
            ddxi    = self%basis_c%interpolator_edge('ddxi',   iedge)
            ddeta   = self%basis_c%interpolator_edge('ddeta',  iedge)
            ddzeta  = self%basis_c%interpolator_edge('ddzeta', iedge)

            !
            ! Compute coordinate jacobian matrix at interpolation nodes
            !
            jacobian(1,1,:) = matmul(ddxi,   self%coords%getvar(1,itime = 1))
            jacobian(1,2,:) = matmul(ddeta,  self%coords%getvar(1,itime = 1))
            jacobian(1,3,:) = matmul(ddzeta, self%coords%getvar(1,itime = 1))

            jacobian(2,1,:) = matmul(ddxi,   self%coords%getvar(2,itime = 1))
            jacobian(2,2,:) = matmul(ddeta,  self%coords%getvar(2,itime = 1))
            jacobian(2,3,:) = matmul(ddzeta, self%coords%getvar(2,itime = 1))

            jacobian(3,1,:) = matmul(ddxi,   self%coords%getvar(3,itime = 1))
            jacobian(3,2,:) = matmul(ddeta,  self%coords%getvar(3,itime = 1))
            jacobian(3,3,:) = matmul(ddzeta, self%coords%getvar(3,itime = 1))


            !
            ! Add coordinate system scaling to jacobian matrix
            !
            select case (self%coordinate_system)
                case (CARTESIAN)
                    scaling_row2 = ONE
                case (CYLINDRICAL)
                    ! TODO: probably wrong number of nodes here for edges!!!!!
                    scaling_row2 = self%interp_coords(:,1)
                case default
                    user_msg = "element%interpolate_metrics_edge: Invalid coordinate system."
                    call chidg_signal(FATAL,user_msg)
            end select



            !
            ! Apply coorindate system scaling
            !
            jacobian(2,1,:) = jacobian(2,1,:)*scaling_row2
            jacobian(2,2,:) = jacobian(2,2,:)*scaling_row2
            jacobian(2,3,:) = jacobian(2,3,:)*scaling_row2


            !
            ! Invert jacobian matrix at each interpolation node
            !
            do inode = 1,nnodes_edge
                self%edge_metric(:,:,inode,iedge) = inv_3x3(jacobian(:,:,inode))
            end do

        end do !iedge


    end subroutine interpolate_metrics_edge
    !*************************************************************************************







    !>  Compute matrices containing gradients of basis/test function
    !!  at each quadrature node.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine interpolate_gradients(self)
        class(element_t),   intent(inout)   :: self

        character(:),   allocatable :: user_msg
        integer(ik)                 :: iterm,inode

        integer(ik)                                 :: nnodes
        real(rk),   allocatable,    dimension(:,:)  :: ddxi, ddeta, ddzeta

        nnodes = self%basis_s%nnodes_elem()
        ddxi   = self%basis_s%interpolator_element('ddxi')
        ddeta  = self%basis_s%interpolator_element('ddeta')
        ddzeta = self%basis_s%interpolator_element('ddzeta')

        do iterm = 1,self%nterms_s
            do inode = 1,nnodes
                self%grad1(inode,iterm) = self%metric(1,1,inode) * ddxi(inode,iterm)  + &
                                          self%metric(2,1,inode) * ddeta(inode,iterm) + &
                                          self%metric(3,1,inode) * ddzeta(inode,iterm)

                self%grad2(inode,iterm) = self%metric(1,2,inode) * ddxi(inode,iterm)  + &
                                          self%metric(2,2,inode) * ddeta(inode,iterm) + &
                                          self%metric(3,2,inode) * ddzeta(inode,iterm)

                self%grad3(inode,iterm) = self%metric(1,3,inode) * ddxi(inode,iterm)  + &
                                          self%metric(2,3,inode) * ddeta(inode,iterm) + &
                                          self%metric(3,3,inode) * ddzeta(inode,iterm)
            end do
        end do

        !
        ! Check for acceptable element
        !
        if (any(ieee_is_nan(self%grad1)) .or. &
            any(ieee_is_nan(self%grad2)) .or. &
            any(ieee_is_nan(self%grad3)) ) then
            user_msg = "element%interpolate_gradients: Element failed to produce valid &
                        gradient information. Element quality is likely not reasonable."
            call chidg_signal(FATAL,user_msg)
        end if

        self%grad1_trans = transpose(self%grad1)
        self%grad2_trans = transpose(self%grad2)
        self%grad3_trans = transpose(self%grad3)

    end subroutine interpolate_gradients
    !******************************************************************************************



    !>  Compute matrices containing gradients of basis/test function
    !!  at each quadrature node.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine interpolate_gradients_edge(self)
        class(element_t),   intent(inout)   :: self

        character(:),   allocatable :: user_msg
        integer(ik)                 :: iterm, inode, iedge, nnodes_edge
        real(rk), allocatable, dimension(:,:)   :: ddxi, ddeta, ddzeta

        nnodes_edge = self%basis_s%nnodes_edge()

        do iedge = 1,NEDGES

            ddxi   = self%basis_s%interpolator_edge('ddxi',  iedge)
            ddeta  = self%basis_s%interpolator_edge('ddeta', iedge)
            ddzeta = self%basis_s%interpolator_edge('ddzeta',iedge)

            do iterm = 1,self%nterms_s
                do inode = 1,nnodes_edge
                    self%edge_grad1(inode,iterm,iedge) = self%edge_metric(1,1,inode,iedge) * ddxi(inode,iterm)  + &
                                                         self%edge_metric(2,1,inode,iedge) * ddeta(inode,iterm) + &
                                                         self%edge_metric(3,1,inode,iedge) * ddzeta(inode,iterm)

                    self%edge_grad2(inode,iterm,iedge) = self%edge_metric(1,2,inode,iedge) * ddxi(inode,iterm)  + &
                                                         self%edge_metric(2,2,inode,iedge) * ddeta(inode,iterm) + &
                                                         self%edge_metric(3,2,inode,iedge) * ddzeta(inode,iterm)

                    self%edge_grad3(inode,iterm,iedge) = self%edge_metric(1,3,inode,iedge) * ddxi(inode,iterm)  + &
                                                         self%edge_metric(2,3,inode,iedge) * ddeta(inode,iterm) + &
                                                         self%edge_metric(3,3,inode,iedge) * ddzeta(inode,iterm)
                end do
            end do

        end do !iedge

        !
        ! Check for acceptable element
        !
        if (any(ieee_is_nan(self%edge_grad1)) .or. &
            any(ieee_is_nan(self%edge_grad2)) .or. &
            any(ieee_is_nan(self%edge_grad3)) ) then
            user_msg = "element%interpolate_gradients_edge: Element failed to produce valid &
                        gradient information. Element quality is likely not reasonable."
            call chidg_signal(FATAL,user_msg)
        end if


    end subroutine interpolate_gradients_edge
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
    subroutine interpolate_coords(self)
        class(element_t),   intent(inout)   :: self

        integer(ik)                             :: nnodes
        real(rk),   dimension(:),   allocatable :: coord1, coord2, coord3
        integer(ik)                             :: inode

        nnodes = self%basis_c%nnodes_elem()

        !
        ! compute coordinates associated with quadrature points
        !
        coord1 = matmul(self%basis_c%interpolator_element('Value'),self%coords%getvar(1,itime = 1))
        coord2 = matmul(self%basis_c%interpolator_element('Value'),self%coords%getvar(2,itime = 1))
        coord3 = matmul(self%basis_c%interpolator_element('Value'),self%coords%getvar(3,itime = 1))


        !
        ! Initialize each point with coordinates
        !
        do inode = 1,nnodes
            self%interp_coords(inode,1:3) = [coord1(inode), coord2(inode), coord3(inode)]
        end do

    end subroutine interpolate_coords
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

        integer(ik)             :: iterm
        real(rk),   allocatable :: temp(:,:)

        self%invmass = ZERO
        self%mass    = ZERO

        temp = transpose(self%basis_s%interpolator_element('Value'))


        !
        ! Multiply rows by quadrature weights and cell jacobians
        !
        do iterm = 1,self%nterms_s
            temp(iterm,:) = temp(iterm,:)*(self%basis_s%weights_element())*(self%jinv)
        end do


        !
        ! Perform the matrix multiplication of the transpose val matrix by
        ! the standard matrix. This produces the mass matrix. I think...
        !
        self%mass = matmul(temp,self%basis_s%interpolator_element('Value'))


        !
        ! Compute and store the inverted mass matrix
        !
        self%invmass = inv(self%mass)



    end subroutine compute_mass_matrix
    !******************************************************************************************






    !>  Initialize ALE data from nodal displacements.
    !!
    !!  With domain-global representations of node displacements and velocities passed
    !!  in, update the ALE data by:
    !!      1: locate displacements and velocities in the global array that correspond to the local element.
    !!      2: (re)allocate ale storage and initialize data containers
    !!      3: compute ALE element nodes as: (ref_nodes + displacements). Compute modal representation.
    !!      4: set ALE element node velocities. Comput modal representation.
    !!
    !!  @author Eric Wolf (AFRL)
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2017
    !!
    !-------------------------------------------------------------------------------------
    subroutine set_displacements_velocities(self,dnodes,vnodes)
        class(element_t),   intent(inout)   :: self
        real(rk),           intent(in)      :: dnodes(:,:)
        real(rk),           intent(in)      :: vnodes(:,:)

        integer(ik)             :: ipt, npts, inode, ierr
        real(rk),   allocatable :: modes1(:), modes2(:), modes3(:), dnodes_l(:,:)


        ! Check if reference geometry has been initialized
        if (.not. self%geom_initialized) call chidg_signal(FATAL,"element%set_displacements_velocities: reference geometry has not been initialized. Make sure to call element%init_geom.")


        !
        ! Accumulate local node displacements and velocities
        !
        npts = size(self%node_coords,1)
        if (allocated(self%node_coords_vel)) deallocate(self%node_coords_vel)
        allocate(dnodes_l(npts,3), self%node_coords_vel(npts,3), stat=ierr)
        if (ierr /= 0) call AllocationError

        do ipt = 1,npts
            ! Get node index
            inode = self%connectivity(ipt)

            ! Assemble local node disp/vel from global
            dnodes_l(ipt,:)  = dnodes(inode,:)
            self%node_coords_vel(ipt,:)  = vnodes(inode,:)
        end do !ipt


        !
        ! Compute deformed coordinates at element nodes: (undeformed + perturbation)
        !
        self%node_coords_def = (self%node_coords + dnodes_l)


        !
        ! ALE (re)initialization
        !
        call self%coords_def%init(self%nterms_c,self%spacedim,self%ntime,self%idomain_g,self%idomain_l,self%ielement_g,self%ielement_l)
        call self%coords_vel%init(self%nterms_c,self%spacedim,self%ntime,self%idomain_g,self%idomain_l,self%ielement_g,self%ielement_l)


        !
        ! Compute modal representation of coordinates
        !
        modes1 = matmul(self%nodes_to_modes,self%node_coords_def(:,1))
        modes2 = matmul(self%nodes_to_modes,self%node_coords_def(:,2))
        modes3 = matmul(self%nodes_to_modes,self%node_coords_def(:,3))

        call self%coords_def%setvar(1,itime = 1,vals = modes1)
        call self%coords_def%setvar(2,itime = 1,vals = modes2)
        call self%coords_def%setvar(3,itime = 1,vals = modes3)



        !
        ! Compute modal representation of grid velocity
        !
        modes1 = matmul(self%nodes_to_modes,self%node_coords_vel(:,1))
        modes2 = matmul(self%nodes_to_modes,self%node_coords_vel(:,2))
        modes3 = matmul(self%nodes_to_modes,self%node_coords_vel(:,3))

        call self%coords_vel%setvar(1,itime = 1,vals = modes1)
        call self%coords_vel%setvar(2,itime = 1,vals = modes2)
        call self%coords_vel%setvar(3,itime = 1,vals = modes3)


    end subroutine set_displacements_velocities
    !**************************************************************************************









    !>
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine update_interpolations_ale(self)
        class(element_t),       intent(inout)      :: self

        call self%interpolate_coords_ale()
        call self%interpolate_metrics_ale()

    end subroutine update_interpolations_ale
    !*******************************************************************************************






    !>  Compute ALE quantities at interpolation nodes from modes.
    !!      coords_def     -> interp_coords_def
    !!      coords_vel -> interp_coords_vel
    !!
    !!  @author Eric Wolf (AFRL)
    !!  @date   7/14/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine interpolate_coords_ale(self)
        class(element_t),   intent(inout)   :: self

        real(rk),   allocatable :: val(:,:)


        !
        ! Retrieve interpolator
        !
        val = self%basis_c%interpolator_element('Value')

        !
        ! compute cartesian coordinates associated with quadrature points
        !
        self%interp_coords_def(:,1) = matmul(val,self%coords_def%getvar(1,itime = 1))
        self%interp_coords_def(:,2) = matmul(val,self%coords_def%getvar(2,itime = 1))
        self%interp_coords_def(:,3) = matmul(val,self%coords_def%getvar(3,itime = 1))


        !
        ! Grid velocity
        ! compute cartesian coordinates associated with quadrature points
        !
        self%interp_coords_vel(:,1) = matmul(val,self%coords_vel%getvar(1,itime = 1))
        self%interp_coords_vel(:,2) = matmul(val,self%coords_vel%getvar(2,itime = 1))
        self%interp_coords_vel(:,3) = matmul(val,self%coords_vel%getvar(3,itime = 1))

    end subroutine interpolate_coords_ale
    !****************************************************************************************






    !>  Compute ALE metrics at interpolation nodes.
    !!
    !!  @author Eric Wolf (AFRL)
    !!
    !----------------------------------------------------------------------------------------
    subroutine interpolate_metrics_ale(self)
        class(element_t),    intent(inout)   :: self

        integer(ik)                             :: inode, nnodes, ierr
        character(:),               allocatable :: coordinate_system, user_msg

        real(rk),   dimension(:),   allocatable ::                  &
            dd1_dxidxi,   dd1_detadeta,   dd1_dzetadzeta,           &
            dd2_dxidxi,   dd2_detadeta,   dd2_dzetadzeta,           &
            dd3_dxidxi,   dd3_detadeta,   dd3_dzetadzeta,           &
            dd1_dxideta,  dd1_dxidzeta,   dd1_detadzeta,            &
            dd2_dxideta,  dd2_dxidzeta,   dd2_detadzeta,            &
            dd3_dxideta,  dd3_dxidzeta,   dd3_detadzeta,            &
            ale_g_ddxi,   ale_g_ddeta,    ale_g_ddzeta,             &
            jinv_grad1, jinv_grad2, jinv_grad3,   &
            jinv_def_grad1,   jinv_def_grad2,   jinv_def_grad3,     &
            fvals, temp, scaling_row2, weights

        real(rk),   dimension(:,:), allocatable ::  &
            val,                                    &
            ddxi,    ddeta,    ddzeta,              &
            dxidxi,  detadeta, dzetadzeta,          &
            dxideta, dxidzeta, detadzeta,           &
            D_matrix

        real(rk), dimension(:,:,:), allocatable ::  &
            jacobian_ale, jacobian(:,:,:)



        !
        ! Retrieve interpolators
        !
        nnodes  = self%basis_c%nnodes_elem()
        weights = self%basis_c%weights_element()
        ddxi    = self%basis_c%interpolator_element('ddxi')
        ddeta   = self%basis_c%interpolator_element('ddeta')
        ddzeta  = self%basis_c%interpolator_element('ddzeta')

        dxidxi     = self%basis_c%interpolator_element('dxidxi')
        detadeta   = self%basis_c%interpolator_element('detadeta')
        dzetadzeta = self%basis_c%interpolator_element('dzetadzeta')

        dxideta    = self%basis_c%interpolator_element('dxideta')
        dxidzeta   = self%basis_c%interpolator_element('dxidzeta')
        detadzeta  = self%basis_c%interpolator_element('detadzeta')



        ! First derivatives
        allocate(jacobian_ale(3,3,nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError
        jacobian_ale(1,1,:) = matmul(ddxi,   self%coords_def%getvar(1,itime = 1))
        jacobian_ale(1,2,:) = matmul(ddeta,  self%coords_def%getvar(1,itime = 1))
        jacobian_ale(1,3,:) = matmul(ddzeta, self%coords_def%getvar(1,itime = 1))

        jacobian_ale(2,1,:) = matmul(ddxi,   self%coords_def%getvar(2,itime = 1))
        jacobian_ale(2,2,:) = matmul(ddeta,  self%coords_def%getvar(2,itime = 1))
        jacobian_ale(2,3,:) = matmul(ddzeta, self%coords_def%getvar(2,itime = 1))

        jacobian_ale(3,1,:) = matmul(ddxi,   self%coords_def%getvar(3,itime = 1))
        jacobian_ale(3,2,:) = matmul(ddeta,  self%coords_def%getvar(3,itime = 1))
        jacobian_ale(3,3,:) = matmul(ddzeta, self%coords_def%getvar(3,itime = 1))


        !
        ! Add coordinate system scaling to jacobian matrix
        !
        allocate(scaling_row2(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError

        select case (self%coordinate_system)
            case (CARTESIAN)
                scaling_row2 = ONE
            case (CYLINDRICAL)
                scaling_row2 = self%interp_coords_def(:,1)
            case default
                user_msg = "element%interpolate_metrics_ale: Invalid coordinate system."
                call chidg_signal(FATAL,user_msg)
        end select



        !
        ! Apply coorindate system scaling
        !
        jacobian_ale(2,1,:) = jacobian_ale(2,1,:)*scaling_row2
        jacobian_ale(2,2,:) = jacobian_ale(2,2,:)*scaling_row2
        jacobian_ale(2,3,:) = jacobian_ale(2,3,:)*scaling_row2




        !
        ! Compute inverse cell mapping jacobian
        !
        do inode = 1,nnodes
            self%jinv_def(inode) = det_3x3(jacobian_ale(:,:,inode))
        end do



        !
        ! Check for negative jacobians
        !
        user_msg = "element%interpolate_metrics_ale: Negative element jacobians. &
                    Check element quality and origntation."
        if (any(self%jinv_def < ZERO)) call chidg_signal(FATAL,user_msg)


        !
        ! Compute element volume
        !
        self%vol_ale = abs(sum(self%jinv_def * weights))


        !
        ! Compute element deformation gradient: dX/dx
        !   dX/dx = [dxi/dx][dX/dxi]
        !
        do inode = 1,nnodes
            D_matrix = matmul(jacobian_ale(:,:,inode),self%metric(:,:,inode))
            self%ale_Dinv(:,:,inode) = inv_3x3(D_matrix)

            ! Invert jacobian_ale for use later in routine
            jacobian_ale(:,:,inode) = inv_3x3(jacobian_ale(:,:,inode))
        end do


        !
        ! Compute volume scaling: deformed/undeformed
        !
        self%ale_g = self%jinv_def/self%jinv


        !
        ! Project ale_g to solution basis
        !
        val  = self%basis_s%interpolator_element('Value')
        fvals = self%ale_g * weights * self%jinv
        temp = matmul(transpose(val),fvals)
        self%ale_g_modes = matmul(self%invmass,temp)



        ! Second/mixed derivatives
        dd1_dxidxi     = matmul(dxidxi,     self%coords%getvar(1,itime = 1))
        dd1_detadeta   = matmul(detadeta,   self%coords%getvar(1,itime = 1))
        dd1_dzetadzeta = matmul(dzetadzeta, self%coords%getvar(1,itime = 1))
        dd1_dxideta    = matmul(dxideta,    self%coords%getvar(1,itime = 1))
        dd1_dxidzeta   = matmul(dxidzeta,   self%coords%getvar(1,itime = 1))
        dd1_detadzeta  = matmul(detadzeta,  self%coords%getvar(1,itime = 1))

        dd2_dxidxi     = matmul(dxidxi,     self%coords%getvar(2,itime = 1))
        dd2_detadeta   = matmul(detadeta,   self%coords%getvar(2,itime = 1))
        dd2_dzetadzeta = matmul(dzetadzeta, self%coords%getvar(2,itime = 1))
        dd2_dxideta    = matmul(dxideta,    self%coords%getvar(2,itime = 1))
        dd2_dxidzeta   = matmul(dxidzeta,   self%coords%getvar(2,itime = 1))
        dd2_detadzeta  = matmul(detadzeta,  self%coords%getvar(2,itime = 1))

        dd3_dxidxi     = matmul(dxidxi,     self%coords%getvar(3,itime = 1))
        dd3_detadeta   = matmul(detadeta,   self%coords%getvar(3,itime = 1))
        dd3_dzetadzeta = matmul(dzetadzeta, self%coords%getvar(3,itime = 1))
        dd3_dxideta    = matmul(dxideta,    self%coords%getvar(3,itime = 1))
        dd3_dxidzeta   = matmul(dxidzeta,   self%coords%getvar(3,itime = 1))
        dd3_detadzeta  = matmul(detadzeta,  self%coords%getvar(3,itime = 1))

        jinv_grad1 = dd1_dxidxi*self%metric(1,1,:)    +  dd1_dxideta*self%metric(2,1,:)    +  dd1_dxidzeta*self%metric(3,1,:)   +  &
                     dd2_dxidxi*self%metric(1,2,:)    +  dd2_dxideta*self%metric(2,2,:)    +  dd2_dxidzeta*self%metric(3,2,:)   +  &
                     dd3_dxidxi*self%metric(1,3,:)    +  dd3_dxideta*self%metric(2,3,:)    +  dd3_dxidzeta*self%metric(3,3,:)

        jinv_grad2 = dd1_dxideta*self%metric(1,1,:)   +  dd1_detadeta*self%metric(2,1,:)   +  dd1_detadzeta*self%metric(3,1,:)  +  &
                     dd2_dxideta*self%metric(1,2,:)   +  dd2_detadeta*self%metric(2,2,:)   +  dd2_detadzeta*self%metric(3,2,:)  +  &
                     dd3_dxideta*self%metric(1,3,:)   +  dd3_detadeta*self%metric(2,3,:)   +  dd3_detadzeta*self%metric(3,3,:)

        jinv_grad3 = dd1_dxidzeta*self%metric(1,1,:)  +  dd1_detadzeta*self%metric(2,1,:)  +  dd1_dzetadzeta*self%metric(3,1,:) +  &
                     dd2_dxidzeta*self%metric(1,2,:)  +  dd2_detadzeta*self%metric(2,2,:)  +  dd2_dzetadzeta*self%metric(3,2,:) +  &
                     dd3_dxidzeta*self%metric(1,3,:)  +  dd3_detadzeta*self%metric(2,3,:)  +  dd3_dzetadzeta*self%metric(3,3,:)


        ! Second/mixed derivatives
        dd1_dxidxi     = matmul(dxidxi,     self%coords_def%getvar(1,itime = 1))
        dd1_detadeta   = matmul(detadeta,   self%coords_def%getvar(1,itime = 1))
        dd1_dzetadzeta = matmul(dzetadzeta, self%coords_def%getvar(1,itime = 1))
        dd1_dxideta    = matmul(dxideta,    self%coords_def%getvar(1,itime = 1))
        dd1_dxidzeta   = matmul(dxidzeta,   self%coords_def%getvar(1,itime = 1))
        dd1_detadzeta  = matmul(detadzeta,  self%coords_def%getvar(1,itime = 1))

        dd2_dxidxi     = matmul(dxidxi,     self%coords_def%getvar(2,itime = 1))
        dd2_detadeta   = matmul(detadeta,   self%coords_def%getvar(2,itime = 1))
        dd2_dzetadzeta = matmul(dzetadzeta, self%coords_def%getvar(2,itime = 1))
        dd2_dxideta    = matmul(dxideta,    self%coords_def%getvar(2,itime = 1))
        dd2_dxidzeta   = matmul(dxidzeta,   self%coords_def%getvar(2,itime = 1))
        dd2_detadzeta  = matmul(detadzeta,  self%coords_def%getvar(2,itime = 1))

        dd3_dxidxi     = matmul(dxidxi,     self%coords_def%getvar(3,itime = 1))
        dd3_detadeta   = matmul(detadeta,   self%coords_def%getvar(3,itime = 1))
        dd3_dzetadzeta = matmul(dzetadzeta, self%coords_def%getvar(3,itime = 1))
        dd3_dxideta    = matmul(dxideta,    self%coords_def%getvar(3,itime = 1))
        dd3_dxidzeta   = matmul(dxidzeta,   self%coords_def%getvar(3,itime = 1))
        dd3_detadzeta  = matmul(detadzeta,  self%coords_def%getvar(3,itime = 1))

        jinv_def_grad1 = dd1_dxidxi*jacobian_ale(1,1,:)    +  dd1_dxideta*jacobian_ale(2,1,:)    +  dd1_dxidzeta*jacobian_ale(3,1,:)   +  &
                         dd2_dxidxi*jacobian_ale(1,2,:)    +  dd2_dxideta*jacobian_ale(2,2,:)    +  dd2_dxidzeta*jacobian_ale(3,2,:)   +  &
                         dd3_dxidxi*jacobian_ale(1,3,:)    +  dd3_dxideta*jacobian_ale(2,3,:)    +  dd3_dxidzeta*jacobian_ale(3,3,:)

        jinv_def_grad2 = dd1_dxideta*jacobian_ale(1,1,:)   +  dd1_detadeta*jacobian_ale(2,1,:)   +  dd1_detadzeta*jacobian_ale(3,1,:)  +  &
                         dd2_dxideta*jacobian_ale(1,2,:)   +  dd2_detadeta*jacobian_ale(2,2,:)   +  dd2_detadzeta*jacobian_ale(3,2,:)  +  &
                         dd3_dxideta*jacobian_ale(1,3,:)   +  dd3_detadeta*jacobian_ale(2,3,:)   +  dd3_detadzeta*jacobian_ale(3,3,:)

        jinv_def_grad3 = dd1_dxidzeta*jacobian_ale(1,1,:)  +  dd1_detadzeta*jacobian_ale(2,1,:)  +  dd1_dzetadzeta*jacobian_ale(3,1,:) +  &
                         dd2_dxidzeta*jacobian_ale(1,2,:)  +  dd2_detadzeta*jacobian_ale(2,2,:)  +  dd2_dzetadzeta*jacobian_ale(3,2,:) +  &
                         dd3_dxidzeta*jacobian_ale(1,3,:)  +  dd3_detadzeta*jacobian_ale(2,3,:)  +  dd3_dzetadzeta*jacobian_ale(3,3,:)


        !
        ! Apply Quotiend Rule for computing gradient of det_jacobian_grid
        !
        !   det_jacobian_grid = jinv_def/jinv
        !
        !   grad(det_jacobian_grid) = [grad(jinv_def)*jinv - jinv_def*grad(jinv)] / [jinv*jinv]
        !
        ale_g_ddxi   = (jinv_def_grad1*self%jinv  -  self%jinv_def*jinv_grad1)/(self%jinv**TWO)
        ale_g_ddeta  = (jinv_def_grad2*self%jinv  -  self%jinv_def*jinv_grad2)/(self%jinv**TWO)
        ale_g_ddzeta = (jinv_def_grad3*self%jinv  -  self%jinv_def*jinv_grad3)/(self%jinv**TWO)


        ! Transform into gradient in physical space(undeformed geometry)
        do inode = 1,size(ale_g_ddxi)
            self%ale_g_grad1(inode) = self%metric(1,1,inode) * ale_g_ddxi(inode)  + &
                                      self%metric(2,1,inode) * ale_g_ddeta(inode) + &
                                      self%metric(3,1,inode) * ale_g_ddzeta(inode)

            self%ale_g_grad2(inode) = self%metric(1,2,inode) * ale_g_ddxi(inode)  + &
                                      self%metric(2,2,inode) * ale_g_ddeta(inode) + &
                                      self%metric(3,2,inode) * ale_g_ddzeta(inode)

            self%ale_g_grad3(inode) = self%metric(1,3,inode) * ale_g_ddxi(inode)  + &
                                      self%metric(2,3,inode) * ale_g_ddeta(inode) + &
                                      self%metric(3,3,inode) * ale_g_ddzeta(inode)
        end do



    end subroutine interpolate_metrics_ale
    !********************************************************************************************************











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

        real(rk)                   :: polyvals(self%nterms_c)
        integer(ik)                :: iterm, spacedim

        ! Evaluate polynomial modes at node location
        spacedim = self%spacedim
        do iterm = 1,self%nterms_c
            polyvals(iterm)  = polynomial_val(spacedim,self%nterms_c,iterm,[xi,eta,zeta])
        end do
        
        ! Evaluate x from dot product of modes and polynomial values
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

        real(rk)                   :: polyvals(self%nterms_c)
        integer(ik)                :: iterm, spacedim


        ! Evaluate polynomial modes at node location
        spacedim = self%spacedim
        do iterm = 1,self%nterms_c
            polyvals(iterm)  = polynomial_val(spacedim,self%nterms_c,iterm,[xi,eta,zeta])
        end do


        ! Evaluate x from dot product of modes and polynomial values
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

        real(rk)                   :: polyvals(self%nterms_c)
        integer(ik)                :: iterm, spacedim


        ! Evaluate polynomial modes at node location
        spacedim = self%spacedim
        do iterm = 1,self%nterms_c
            polyvals(iterm)  = polynomial_val(spacedim,self%nterms_c,iterm,[xi,eta,zeta])
        end do


        ! Evaluate x from dot product of modes and polynomial values
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

        real(rk)    :: val1, val2, val3
        real(rk)    :: phys_point(3)
        real(rk)    :: polyvals(self%nterms_c)
        integer(ik) :: iterm

        !
        ! Evaluate polynomial modes at node location
        !
        do iterm = 1,self%nterms_c
            polyvals(iterm) = polynomial_val(self%spacedim,self%nterms_c,iterm,[xi,eta,zeta])
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
        phys_point = [val1,val2,val3]


    end function physical_point
    !******************************************************************************************









    !>  Compute a coordinate value, based on the location in reference space (xi, eta, zeta)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  elem    Element containing coordinate expansion
    !!  @param[in]  type    'Reference' or 'ALE'
    !!  @param[in]  icoord  Integer corresponding to coordinate index 1(x), 2(y), 3(z)
    !!  @param[in]  xi      Real value for xi-coordinate
    !!  @param[in]  eta     Real value for eta-coordinate
    !!  @param[in]  zeta    Real value for zeta-coordinate
    !!
    !!  @author Mayank Sharma + MAtteo Ugolotti
    !!  @date   11/5/2016
    !!
    !!  TODO: TEST type 'Reference' and 'ALE'
    !!
    !-----------------------------------------------------------------------------------------
    function grid_point(self,coord_index,location,coordinate_frame) result(val)
        class(element_t),   intent(in)  :: self
        integer(ik),        intent(in)  :: coord_index
        real(rk),           intent(in)  :: location(3)
        character(*),       intent(in)  :: coordinate_frame

        real(rk)    :: val, polyvals(self%nterms_c)
        integer(ik) :: iterm

        if (coord_index > 3)             call chidg_signal(FATAL,"element%grid_point -- coord_index exceeded 3 physical coordinates")
        if (.not. self%geom_initialized) call chidg_signal(FATAL,"element%grid_point: geometry not initialized")


        ! Evaluate polynomial modes at node location
        do iterm = 1,self%nterms_c
            polyvals(iterm) = polynomial_val(self%spacedim,self%nterms_c,iterm,location)
        end do


        ! Evaluate mesh point from dot product of modes and polynomial values
        select case(trim(coordinate_frame))
            case('Undeformed')
                val = dot_product(self%coords%getvar(coord_index,itime = 1), polyvals)
            case('Deformed')
                val = dot_product(self%coords_def%getvar(coord_index,itime = 1), polyvals)
            case default
                call chidg_signal(FATAL,"element%grid_point: Invalid input for coordinate_frame.")
        end select


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
    function metric_point(self,coord,coordinate_frame,coordinate_scaling) result(res)
        class(element_t),   intent(in)              :: self
        real(rk),           intent(in)              :: coord(3)
        character(*),       intent(in), optional    :: coordinate_frame
        logical,            intent(in), optional    :: coordinate_scaling
        
        real(rk)                            :: r, res(3,3), dXdxi(3,3)
        real(rk), dimension(self%nterms_c)  :: ddxi, ddeta, ddzeta
        integer(ik)                         :: iterm
        logical                             :: scale_metric
        type(densevector_t)                 :: modal_coords
        character(:),   allocatable         :: frame_selector


        !
        ! Set default values for optional inputs
        !
        frame_selector = 'Undeformed'
        scale_metric   = .true.


        !
        ! Overwrite optional values if inputs are present
        !
        if (present(coordinate_scaling)) scale_metric   = coordinate_scaling
        if (present(coordinate_frame)  ) frame_selector = coordinate_frame


        !
        ! Retrieve modal coordinates
        !
        select case(trim(frame_selector))
            case('Undeformed')
                modal_coords = self%coords
            case('Deformed')
                modal_coords = self%coords_def
        end select


        !
        ! Evaluate basis modes at node location
        !
        do iterm = 1,self%nterms_c
            ddxi(iterm)   = dpolynomial_val(self%spacedim,self%nterms_c,iterm,coord,XI_DIR  )
            ddeta(iterm)  = dpolynomial_val(self%spacedim,self%nterms_c,iterm,coord,ETA_DIR )
            ddzeta(iterm) = dpolynomial_val(self%spacedim,self%nterms_c,iterm,coord,ZETA_DIR)
        end do


        !
        ! Evaluate mesh point from dot product of modes and polynomial values
        !
        dXdxi(1,1) = dot_product(ddxi,   modal_coords%getvar(ivar=1, itime=1))
        dXdxi(1,2) = dot_product(ddeta,  modal_coords%getvar(ivar=1, itime=1))
        dXdxi(1,3) = dot_product(ddzeta, modal_coords%getvar(ivar=1, itime=1))
                                         
        dXdxi(2,1) = dot_product(ddxi,   modal_coords%getvar(ivar=2, itime=1))
        dXdxi(2,2) = dot_product(ddeta,  modal_coords%getvar(ivar=2, itime=1))
        dXdxi(2,3) = dot_product(ddzeta, modal_coords%getvar(ivar=2, itime=1))
                                         
        dXdxi(3,1) = dot_product(ddxi,   modal_coords%getvar(ivar=3, itime=1))
        dXdxi(3,2) = dot_product(ddeta,  modal_coords%getvar(ivar=3, itime=1))
        dXdxi(3,3) = dot_product(ddzeta, modal_coords%getvar(ivar=3, itime=1))


        !
        ! Apply coordinate system scaling
        !
        if (scale_metric) then
            select case(self%coordinate_system)
                case(CYLINDRICAL) 
                    select case(frame_selector)
                        case('Undeformed')
                            r = self%grid_point(DIR_R,coord,frame_selector)
                        case('Deformed')
                            r = self%grid_point(DIR_R,coord,frame_selector)
                    end select
                    dXdxi(2,:) = dXdxi(2,:) * r
            end select
        end if


        !
        ! Invert transformation
        !
        res = inv_3x3(dXdxi)

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

        real(rk)    :: temp,val
        real(rk)    :: ale_g,ale_g_grad(3),ale_dinv(3,3),interp_coords_vel(3)
        real(rk)    :: polyvals(q%nterms())
        integer(ik) :: iterm, spacedim


        ! evaluate polynomial modes at node location
        spacedim = self%spacedim
        do iterm = 1,q%nterms()
            polyvals(iterm)  = polynomial_val(spacedim,q%nterms(),iterm,[xi,eta,zeta])
        end do


        ! evaluate x from dot product of modes and polynomial values
        temp = dot_product(q%getvar(ivar,itime),polyvals)

        !call self%ale_point(xi,eta,zeta,ale_g,ale_g_grad,ale_Dinv,grid_vel)
        !val = temp/det_jacobian_grid
        val = temp/dot_product(self%ale_g_modes,polyvals)


    end function solution_point
    !******************************************************************************************




!    !>  Compute a variable value, based on the location in reference space (xi, eta, zeta)
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   2/1/2016
!    !!
!    !!  @param[in]  elem    Element that the solution expansion is associated with.
!    !!  @param[in]  q       Solution expansion for a given element.
!    !!  @param[in]  ivar    Integer corresponding to variable index.
!    !!  @param[in]  xi      Real value for xi-coordinate.
!    !!  @param[in]  eta     Real value for eta-coordinate.
!    !!  @param[in]  zeta    Real value for zeta-coordinate.
!    !!
!    !!  @author Mayank Sharma + Matteo Ugolotti
!    !!  @date   11/5/2016
!    !!
!    !----------------------------------------------------------------------------------------
!    function derivative_point(self,q,ivar,itime,xi,eta,zeta,dir) result(val)
!        class(element_t),       intent(in)      :: self
!        class(densevector_t),   intent(in)      :: q
!        integer(ik),            intent(in)      :: ivar
!        integer(ik),            intent(in)      :: itime
!        real(rk),               intent(in)      :: xi,eta,zeta
!        integer(ik),            intent(in)      :: dir
!
!        real(rk)        :: val
!        real(rk)        :: ddxi(q%nterms()), ddeta(q%nterms()), ddzeta(q%nterms()), &
!                           deriv(q%nterms())
!        real(rk)        :: metric(3,3), jinv, dxi_dx, dxi_dy, dxi_dz, &
!                           deta_dx, deta_dy, deta_dz, dzeta_dx, dzeta_dy, dzeta_dz
!        integer(ik)     :: iterm, spacedim
!
!
!        !
!        ! Evaluate polynomial mode derivatives at node location
!        !
!        spacedim = self%spacedim
!        do iterm = 1,q%nterms()
!            ddxi(iterm)   = dpolynomial_val(spacedim,q%nterms(),iterm,[xi,eta,zeta],XI_DIR)
!            ddeta(iterm)  = dpolynomial_val(spacedim,q%nterms(),iterm,[xi,eta,zeta],ETA_DIR)
!            ddzeta(iterm) = dpolynomial_val(spacedim,q%nterms(),iterm,[xi,eta,zeta],ZETA_DIR)
!        end do
!
!
!        !
!        ! Compute metrics at node
!        !
!        metric(1,1) = self%metric_point(DIR_1,XI_DIR,  xi,eta,zeta)
!        metric(2,1) = self%metric_point(DIR_2,XI_DIR,  xi,eta,zeta)
!        metric(3,1) = self%metric_point(DIR_3,XI_DIR,  xi,eta,zeta)
!        metric(1,2) = self%metric_point(DIR_1,ETA_DIR, xi,eta,zeta)
!        metric(2,2) = self%metric_point(DIR_2,ETA_DIR, xi,eta,zeta)
!        metric(3,2) = self%metric_point(DIR_3,ETA_DIR, xi,eta,zeta)
!        metric(1,3) = self%metric_point(DIR_1,ZETA_DIR,xi,eta,zeta)
!        metric(2,3) = self%metric_point(DIR_2,ZETA_DIR,xi,eta,zeta)
!        metric(3,3) = self%metric_point(DIR_3,ZETA_DIR,xi,eta,zeta)
!
!
!        !
!        ! Compute inverse cell mapping jacobian
!        !
!        jinv = metric(1,1)*metric(2,2)*metric(3,3) - metric(1,2)*metric(2,1)*metric(3,3) - &
!                     metric(1,1)*metric(2,3)*metric(3,2) + metric(1,3)*metric(2,1)*metric(3,2) + &
!                     metric(1,2)*metric(2,3)*metric(3,1) - metric(1,3)*metric(2,2)*metric(3,1)
!
!
!
!
!
!        do iterm = 1,self%nterms_s
!            if (dir == DIR_1) then
!                dxi_dx   = metric(2,2)*metric(3,3) - metric(2,3)*metric(3,2)
!                deta_dx  = metric(2,3)*metric(3,1) - metric(2,1)*metric(3,3)
!                dzeta_dx = metric(2,1)*metric(3,2) - metric(2,2)*metric(3,1)
!                deriv(iterm) = dxi_dx   * ddxi(iterm)   * (ONE/jinv) + &
!                               deta_dx  * ddeta(iterm)  * (ONE/jinv) + &
!                               dzeta_dx * ddzeta(iterm) * (ONE/jinv)
!
!            else if (dir == DIR_2) then
!                dxi_dy   = metric(1,3)*metric(3,2) - metric(1,2)*metric(3,3)
!                deta_dy  = metric(1,1)*metric(3,3) - metric(1,3)*metric(3,1)
!                dzeta_dy = metric(1,2)*metric(3,1) - metric(1,1)*metric(3,2)
!                deriv(iterm) = dxi_dy   * ddxi(iterm)   * (ONE/jinv) + &
!                               deta_dy  * ddeta(iterm)  * (ONE/jinv) + &
!                               dzeta_dy * ddzeta(iterm) * (ONE/jinv)
!
!            else if (dir == DIR_3) then
!                dxi_dz   = metric(1,2)*metric(2,3) - metric(1,3)*metric(2,2)
!                deta_dz  = metric(1,3)*metric(2,1) - metric(1,1)*metric(2,3)
!                dzeta_dz = metric(1,1)*metric(2,2) - metric(1,2)*metric(2,1)
!                deriv(iterm) = dxi_dz   * ddxi(iterm)   * (ONE/jinv) + &
!                               deta_dz  * ddeta(iterm)  * (ONE/jinv) + &
!                               dzeta_dz * ddzeta(iterm) * (ONE/jinv)
!
!            else
!                call chidg_signal(FATAL,"element%derivative_point: Invalid value for 'dir' parameter. (1,2,3).")
!            end if
!        end do
!
!
!
!        !
!        ! Evaluate x from dot product of modes and polynomial values
!        !
!        val = dot_product(q%getvar(ivar,itime),deriv)
!
!    end function derivative_point
!    !*****************************************************************************************







    
    !>  Compute a computational location(xi,eta,zeta), based on the location in 
    !!  physical space (x,y,z), (r,theta,z)
    !!
    !!  NOTE: Will return a location, even if the newton solve did not converge. If the
    !!        newton solve did not converge, the point values will be NaN, so the result
    !!        should be checked for NaN's before being used. The recommended approach is:
    !!       
    !!      use ieee_arithmetic,    only: ieee_is_nan
    !!      valid_point = (any(ieee_is_nan(result)))
    !!        
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/23/2016
    !!
    !!  @param[in]  coord   Coordinate in physical space (x, y, z)
    !!  @result     coord_comp  Coordinate in the element-local coordinate system[xi,eta,zeta]
    !!
    !-----------------------------------------------------------------------------------------
    function computational_point(self,coord) result(coord_comp)
        class(element_t),   intent(in)  :: self
        real(rk),           intent(in)  :: coord(3)

        integer(ik) :: inewton
        real(rk)    :: coord_comp(3), coord_phys(3), minv(3,3), R(3), dcoord(3), res, tol

        !
        ! Newton iteration to find the donor local coordinates
        !
        tol = 100000._rk*RKTOL
        coord_comp = ZERO
        do inewton = 1,20


            ! Compute local physical coordinates as a function of xi,eta,zeta
            ! Return inverted jacobian matrix at that location.
            coord_phys = self%physical_point(coord_comp(1),coord_comp(2),coord_comp(3))
            minv       = self%metric_point(coord_comp, coordinate_frame='Undeformed', coordinate_scaling=.false.)

            ! Assemble residual vector
            R = -(coord_phys - coord)

            ! Solve linear system for Newton update
            dcoord = matmul(minv,R)

            ! Apply Newton update
            coord_comp = coord_comp + dcoord

            ! Exit if converged
            res = norm2(R)
            if ( res < tol ) then
                exit
            end if


            !
            ! Limit computational coordinates, in case they go out of bounds.
            !
            if ( coord_comp(1) >  ONE ) coord_comp(1) =  ONE
            if ( coord_comp(1) < -ONE ) coord_comp(1) = -ONE
            if ( coord_comp(2) >  ONE ) coord_comp(2) =  ONE
            if ( coord_comp(2) < -ONE ) coord_comp(2) = -ONE
            if ( coord_comp(3) >  ONE ) coord_comp(3) =  ONE
            if ( coord_comp(3) < -ONE ) coord_comp(3) = -ONE

            if ( inewton == 20 ) then
                coord_comp = ieee_value(1._rk,ieee_quiet_nan)
            end if

        end do ! inewton


    end function computational_point
    !*****************************************************************************************






    !>  Compute ALE grid quantities, based on the location in reference space (xi, eta, zeta)
    !!
    !!  @author Eric Wolf
    !!  @date   7/21/2017
    !!
    !!  @param[in]  elem            Element that the solution expansion is associated with.
    !!  @param[in]  xi              Real value for xi-coordinate.
    !!  @param[in]  eta             Real value for eta-coordinate.
    !!  @param[in]  zeta            Real value for zeta-coordinate.
    !!  @param[inout]  ale_g        Determinant of grid Jacobian
    !!  @param[inout]  ale_Dinv     Inverse of grid Jacobian
    !!  @param[inout]  interp_coords_vel Grid velocities
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine ale_point(self,xi,eta,zeta,ale_g,ale_g_grad,ale_Dinv,interp_coords_vel) 
        class(element_t),       intent(in)      :: self
        real(rk),               intent(in)      :: xi,eta,zeta
        real(rk),               intent(inout)   :: ale_g
        real(rk),               intent(inout)   :: ale_g_grad(3)
        real(rk),               intent(inout)   :: ale_Dinv(3,3)
        real(rk),               intent(inout)   :: interp_coords_vel(3)

        integer(ik)                             :: iterm, itime
        real(rk)                                :: metric(3,3), jinv, metric_ale(3,3), jinv_def
        real(rk),   dimension(self%nterms_s)    :: val, ddxi, ddeta, ddzeta, grad1, grad2, grad3

        
        !
        ! Compute metrics at node
        !
        metric     = self%metric_point([xi,eta,zeta], coordinate_frame='Undeformed', coordinate_scaling=.true.)
        metric_ale = self%metric_point([xi,eta,zeta], coordinate_frame='Deformed',   coordinate_scaling=.true.)

        !
        ! Compute inverse cell mapping jacobian
        !
        jinv = ONE/det_3x3(metric)
        jinv_def   = ONE/det_3x3(metric_ale)


        !
        ! Compute volume scaling and deformation gradient
        !
        ale_g = jinv_def/jinv
        ale_Dinv = matmul(metric_ale,inv_3x3(metric))


        ! evaluate polynomial modes at node location
        do iterm = 1,self%nterms_s
            val(iterm)    = polynomial_val( self%spacedim,self%nterms_s,iterm,[xi,eta,zeta])
            ddxi(iterm)   = dpolynomial_val(self%spacedim,self%nterms_s,iterm,[xi,eta,zeta],XI_DIR)
            ddeta(iterm)  = dpolynomial_val(self%spacedim,self%nterms_s,iterm,[xi,eta,zeta],ETA_DIR)
            ddzeta(iterm) = dpolynomial_val(self%spacedim,self%nterms_s,iterm,[xi,eta,zeta],ZETA_DIR)
        end do


        do iterm = 1,self%nterms_s
            grad1(iterm) = metric(1,1) * ddxi(iterm)  + &
                           metric(2,1) * ddeta(iterm) + &
                           metric(3,1) * ddzeta(iterm)

            grad2(iterm) = metric(1,2) * ddxi(iterm)  + &
                           metric(2,2) * ddeta(iterm) + &
                           metric(3,2) * ddzeta(iterm)

            grad3(iterm) = metric(1,3) * ddxi(iterm)  + &
                           metric(2,3) * ddeta(iterm) + &
                           metric(3,3) * ddzeta(iterm)
        end do

        ! Evaluate volume scaling and gradient
        ale_g         = dot_product(val,   self%ale_g_modes)
        ale_g_grad(1) = dot_product(grad1, self%ale_g_modes)
        ale_g_grad(2) = dot_product(grad2, self%ale_g_modes)
        ale_g_grad(3) = dot_product(grad3, self%ale_g_modes)

        ! Evaluate grid velocities
        interp_coords_vel(1) = dot_product(val,self%coords_vel%getvar(1,itime = 1))
        interp_coords_vel(2) = dot_product(val,self%coords_vel%getvar(2,itime = 1))
        interp_coords_vel(3) = dot_product(val,self%coords_vel%getvar(3,itime = 1))

    end subroutine ale_point
    !*****************************************************************************************







    !>  Project a function to the solution basis. Return modal coefficients.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/25/2016
    !!
    !-----------------------------------------------------------------------------------------
    function project(self,fcn) result(fmodes)
        class(element_t),       intent(in)      :: self
        class(function_t),      intent(inout)   :: fcn

        real(rk),       allocatable     :: fmodes(:)

        type(point_t),  allocatable     :: pts(:)
        real(rk),       allocatable     :: fvals(:), temp(:)
        real(rk)                        :: time


        ! Call function for evaluation at quadrature nodes and multiply by quadrature weights
        time  = 0._rk
        fvals = fcn%compute(time,point_t(self%interp_coords)) * self%basis_s%weights_element() * self%jinv


        ! Project
        temp = matmul(transpose(self%basis_s%interpolator_element('Value')),fvals)
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
        integer(ik),    allocatable :: face_indices(:)
        integer(ik)                 :: face_index, cindex, eindex, iface_test
        logical                     :: node_matches, face_match, &
                                       corner_one_in_face, corner_two_in_face, &
                                       corner_three_in_face, corner_four_in_face

        !
        ! Get nodes from connectivity
        !
        !   Given global indices of corner nodes such as:
        !       corner_indices = [23, 24, 37, 38]
        !
        !   We want to find their element-local connectivity indices:
        !       corner_position = [1, 2, 5, 6]
        !   
        !
        do cindex = 1,size(corner_indices)
            do eindex = 1,size(self%connectivity)


                node_matches = (corner_indices(cindex) == self%connectivity(eindex))

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

            face_indices = face_corners(iface_test,:,self%element_type)
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


        user_msg = "element%get_face_from_corners: Couldn't find a face index that matched &
                    the provided corner indices."
        if (.not. face_match) call chidg_signal(FATAL,user_msg)


    end function get_face_from_corners
    !******************************************************************************************




    subroutine destructor(self)
        type(element_t), intent(inout) :: self


    end subroutine



end module type_element
