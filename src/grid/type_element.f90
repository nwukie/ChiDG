module type_element
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,XI_MIN,XI_MAX,ETA_MIN, &
                                      ETA_MAX,ZETA_MIN,ZETA_MAX,ONE,ZERO,THIRD, &
                                      DIR_1, DIR_2, DIR_3, DIR_THETA, XI_DIR, ETA_DIR, ZETA_DIR, &
                                      TWO_DIM, THREE_DIM, RKTOL, VALID_POINT, INVALID_POINT, NO_PMM_ASSIGNED, &
                                      ZERO, CARTESIAN, CYLINDRICAL
    use mod_quadrature,         only: GQ, get_quadrature
    use mod_grid,               only: get_element_mapping, face_corners
    use mod_reference_elements, only: get_reference_element, ref_elems
    use mod_polynomial,         only: polynomial_val, dpolynomial_val
    use mod_inv,                only: inv
    use mod_io,                 only: gq_rule


    use type_point
    use type_densevector,           only: densevector_t
    use type_quadrature,            only: quadrature_t
    use type_function,              only: function_t
    use type_element_connectivity,  only: element_connectivity_t
    use type_reference_element,     only: reference_element_t
    use DNAD_D
    use ieee_arithmetic,            only: ieee_value, ieee_quiet_nan
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

        ! Element location
        integer(ik)     :: element_location(4)              ! [idomain_g, idomain_l, ielement_g, ielement_l], useful for nonblocking send
        integer(ik)     :: idomain_g                        ! Global index of parent domain
        integer(ik)     :: idomain_l                        ! Proc-local index of parent domain
        integer(ik)     :: ielement_g                       ! Domain-global index of element
        integer(ik)     :: ielement_l                       ! Proc-local index of the element

        ! Element data
        integer(ik)     :: element_data(8)                  ! [element_type, spacedim, coordinate_system, neqns, nterms_s, nterms_c, ntime, interpolation_level]
        integer(ik)     :: element_type                     ! 1=linear, 2=quadratic, 3=cubic, 4=quartic, etc.
        integer(ik)     :: spacedim                         ! Number of spatial dimensions for the element
        integer(ik)     :: coordinate_system                ! CARTESIAN, CYLINDRICAL. parameters from mod_constants
        integer(ik)     :: neqns                            ! Number of equations being solved
        integer(ik)     :: nterms_s                         ! Number of terms in solution expansion.  
        integer(ik)     :: nterms_c                         ! Number of terms in coordinate expansion. 
        integer(ik)     :: ntime                            ! Number of time levels in solution
        integer(ik)     :: interpolation_level              ! 1=lowest, 2-> are higher

        ! Element quadrature points, mesh points and modes
        integer(ik),    allocatable     :: connectivity(:)      ! Integer indices of the associated nodes in block node list
        real(rk),       allocatable     :: quad_pts(:,:)        ! Coordinates of discrete quadrature points
        real(rk),       allocatable     :: elem_pts(:,:)        ! Coordinates of discrete points defining element
        real(rk),       allocatable     :: dnodes_l(:,:)        ! Node displacements, local element ordering
        real(rk),       allocatable     :: vnodes_l(:,:)        ! Node velocities,    local element ordering
        real(rk),       allocatable     :: nodes_to_modes(:,:)  ! Transformation matrix for converting nodal values to modal coefficients
        type(densevector_t)             :: coords               ! Modal expansion of coordinates (nterms_var,(x,y,z))

        ! Element metric terms
        real(rk), allocatable           :: metric(:,:,:)                ! metric matrix for each quadrature node    (mat_i,mat_j,quad_pt)
        real(rk), allocatable           :: jinv(:)                      ! volume jacobian at quadrature nodes

        !ALE mesh motion terms - to be computed and updated according to a TBD mesh motion algorithm 
        ! Grid displacements/new Cartesian coordinates analogous to elem_pts
        ! Grid Velocities
        ! Grid motion inverse Jacobian
        ! Grid motion Jacobian determinant
        integer(ik)                     :: pmm_ID = NO_PMM_ASSIGNED


        real(rk), allocatable           :: ale_quad_pts(:,:)
        real(rk), allocatable           :: ale_elem_pts(:,:)
        real(rk), allocatable           :: ale_vel_elem_pts(:,:)
        type(densevector_t)             :: ale_coords               ! Modal representation of cartesian coordinates (nterms_var,(x,y,z))
        type(densevector_t)             :: ale_vel_coords           ! Modal representation of cartesian coordinates (nterms_var,(x,y,z))
        real(rk), allocatable           :: grid_vel(:,:)
        real(rk), allocatable           :: jacobian_grid(:,:,:)
        real(rk), allocatable           :: inv_jacobian_grid(:,:,:)
        real(rk), allocatable           :: det_jacobian_grid(:)
        real(rk), allocatable           :: det_jacobian_grid_grad1(:)
        real(rk), allocatable           :: det_jacobian_grid_grad2(:)
        real(rk), allocatable           :: det_jacobian_grid_grad3(:)
        real(rk), allocatable           :: det_jacobian_grid_modes(:)

        real(rk), allocatable           :: jinv_ale(:)                      ! jacobian terms at quadrature nodes

        ! Matrices of physical gradients of basis/test functions
        real(rk), allocatable           :: grad1(:,:)           ! Grad of basis functions in at quadrature nodes
        real(rk), allocatable           :: grad2(:,:)           ! Grad of basis functions in at quadrature nodes
        real(rk), allocatable           :: grad3(:,:)           ! Grad of basis functions in at quadrature nodes
        real(rk), allocatable           :: grad1_trans(:,:)     ! transpose grad1
        real(rk), allocatable           :: grad2_trans(:,:)     ! transpose grad2
        real(rk), allocatable           :: grad3_trans(:,:)     ! transpose grad3

        ! Reference element and interpolators
        type(reference_element_t),  pointer :: basis_s => null()  ! Pointer to solution basis and interpolator
        type(reference_element_t),  pointer :: basis_c => null()  ! Pointer to coordinate basis and interpolator

        ! Element-local mass, inverse mass matrices
        real(rk), allocatable           :: mass(:,:)        
        real(rk), allocatable           :: invmass(:,:)
        real(rk), allocatable           :: mass_c(:,:)        
        real(rk), allocatable           :: invmass_c(:,:)

        ! Element volume, approx. size of bounding box
        real(rk)                        :: vol
        real(rk)                        :: vol_ale
        real(rk)                        :: h(3)     


        ! A psudo-timestep for each equation in the element. Used in the nonlinear solver. 
        ! Quasi-Newton, for example.
        real(rk),   allocatable         :: dtau(:)

        ! Logical tests
        logical :: geom_initialized = .false.
        logical :: numInitialized  = .false.


    contains

        ! Initialization procedures
        procedure, public   :: init_geom
        procedure, public   :: init_ale
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

        procedure, public   :: metric_point_ale
        procedure, public   :: ale_point

        ! Compute a projection of a function onto the solution basis
        procedure, public   :: project


        ! Get connected face
        procedure, public   :: get_face_from_corners

        ! Private utility procedure
        procedure           :: compute_element_matrices
        procedure           :: compute_mass_matrix
        procedure           :: compute_mass_matrix_c
        procedure           :: compute_quadrature_gradients
        procedure           :: compute_quadrature_metrics
        procedure           :: compute_quadrature_coords


        ! ALE procedures
        procedure, public   :: update_element_ale
        procedure           :: update_geom_ale
        procedure           :: compute_quadrature_coords_ale
        procedure           :: compute_quadrature_metrics_ale

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
    subroutine init_geom(self,nodes,connectivity,etype,location,coord_system)
        class(element_t),               intent(inout)   :: self
        real(rk),                       intent(in)      :: nodes(:,:)
        integer(ik),                    intent(in)      :: connectivity(:)
        integer(ik),                    intent(in)      :: etype
        integer(ik),                    intent(in)      :: location(4)
        character(*),                   intent(in)      :: coord_system

        character(:),   allocatable :: user_msg
        real(rk),       allocatable :: nodes_l(:,:), dnodes(:,:), vnodes(:,:)
        real(rk),       allocatable :: modes1(:), modes2(:), modes3(:)
        real(rk)                    :: xmin, xmax, xwidth,  &
                                       ymin, ymax, ywidth,  &
                                       zmin, zmax, zwidth
        integer(ik)                 :: ierr, nterms_c, ipt, npts_1d, npts, &
                                       mapping, inode, spacedim, ref_ID_c
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
        spacedim=3
        self%nodes_to_modes = self%basis_c%nodes_to_modes
        nterms_c = size(self%nodes_to_modes,1)
        self%nterms_c = nterms_c


        user_msg = "element%init_geom: mapping and points do not match."
        if (nterms_c /= size(nodes_l,1)) call chidg_signal(FATAL,user_msg)


        !
        ! Allocate storage
        !
        allocate(self%elem_pts(nterms_c,3),stat=ierr)
        call self%coords%init(nterms_c,3,ntime,self%idomain_g,self%idomain_l,self%ielement_g,self%ielement_l)
        self%spacedim = spacedim
        self%elem_pts = nodes_l

        
        !
        ! Compute modal expansion of element coordinates
        !
        modes1 = matmul(self%nodes_to_modes,self%elem_pts(:,1))
        modes2 = matmul(self%nodes_to_modes,self%elem_pts(:,2))
        modes3 = matmul(self%nodes_to_modes,self%elem_pts(:,3))

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
        ! ALE initialization
        !   Default: zero displacements/velocities
        !
        dnodes = nodes
        vnodes = nodes
        dnodes = ZERO
        vnodes = ZERO
        call self%init_ale(dnodes,vnodes)


        !
        ! Store element_data(1-2)
        !
        self%element_data(1) = self%element_type
        self%element_data(2) = self%spacedim
        self%element_data(3) = self%coordinate_system



    end subroutine init_geom
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
    subroutine init_ale(self,dnodes,vnodes)
        class(element_t),   intent(inout)   :: self
        real(rk),           intent(in)      :: dnodes(:,:)
        real(rk),           intent(in)      :: vnodes(:,:)

        integer(ik)             :: ipt, npts, inode, ierr
        integer(ik), parameter  :: ntime = 1
        real(rk),   allocatable :: modes1(:), modes2(:), modes3(:)


        ! Check if reference geometry has been initialized
        if (.not. self%geom_initialized) call chidg_signal(FATAL,"element%init_ale: reference geometry has not been initialized. Make sure to call element%init_geom.")



        !
        ! Accumulate local node displacements and velocities
        !
        npts = size(self%elem_pts,1)
        if (allocated(self%dnodes_l)) deallocate(self%dnodes_l, self%vnodes_l)
        allocate(self%dnodes_l(npts,3), self%vnodes_l(npts,3), stat=ierr)
        if (ierr /= 0) call AllocationError

        do ipt = 1,npts
            ! Get node index
            !inode = self%connectivity%get_element_node(ipt)
            inode = self%connectivity(ipt)

            ! Assemble local node disp/vel from global
            self%dnodes_l(ipt,:)  = dnodes(inode,:)
            self%vnodes_l(ipt,:)  = vnodes(inode,:)
        end do !ipt




        !
        ! ALE (re)initialization
        !
        if (allocated(self%ale_elem_pts)) deallocate(self%ale_elem_pts,self%ale_vel_elem_pts)
        allocate(self%ale_elem_pts(self%nterms_c,3), self%ale_vel_elem_pts(self%nterms_c,3),stat=ierr)
        if (ierr /= 0) call AllocationError
        call self%ale_coords%init(    self%nterms_c,3,ntime,self%idomain_g,self%idomain_l,self%ielement_g,self%ielement_l)
        call self%ale_vel_coords%init(self%nterms_c,3,ntime,self%idomain_g,self%idomain_l,self%ielement_g,self%ielement_l)



        ! Compute ALE grid nodes: (reference + perturbation)
        self%ale_elem_pts = (self%elem_pts + self%dnodes_l)

        ! Compute ALE grid modes
        modes1 = matmul(self%nodes_to_modes,self%ale_elem_pts(:,1))
        modes2 = matmul(self%nodes_to_modes,self%ale_elem_pts(:,2))
        modes3 = matmul(self%nodes_to_modes,self%ale_elem_pts(:,3))
        call self%ale_coords%setvar(1,itime = 1,vals = modes1)
        call self%ale_coords%setvar(2,itime = 1,vals = modes2)
        call self%ale_coords%setvar(3,itime = 1,vals = modes3)
        



        ! Set ALE velocities
        self%ale_vel_elem_pts = self%vnodes_l

        ! Compute ALE velocity modes
        modes1 = matmul(self%nodes_to_modes,self%ale_vel_elem_pts(:,1))
        modes2 = matmul(self%nodes_to_modes,self%ale_vel_elem_pts(:,2))
        modes3 = matmul(self%nodes_to_modes,self%ale_vel_elem_pts(:,3))
        call self%ale_vel_coords%setvar(1,itime = 1,vals = modes1)
        call self%ale_vel_coords%setvar(2,itime = 1,vals = modes2)
        call self%ale_vel_coords%setvar(3,itime = 1,vals = modes3)



    end subroutine init_ale
    !**************************************************************************************










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
    subroutine init_sol(self,interpolation,level,nterms_s,nfields,ntime)
        class(element_t),   intent(inout) :: self
        character(*),       intent(in)    :: interpolation
        integer(ik),        intent(in)    :: level
        integer(ik),        intent(in)    :: nterms_s
        integer(ik),        intent(in)    :: nfields
        integer(ik),        intent(in)    :: ntime

        integer(ik) :: ierr
        integer(ik) :: nnodes
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
                                         !nterms_rule  = max(self%nterms_c, nterms_s))
                                         nterms_rule  = nterms_s)
        ref_ID_c = get_reference_element(element_type = self%element_type,  &
                                         polynomial   = 'Legendre',         &
                                         nterms       = self%nterms_c,      &
                                         node_set     = interpolation,      &
                                         level        = level,              &
                                         !nterms_rule  = max(self%nterms_c, nterms_s))
                                         nterms_rule  = nterms_s)
        self%basis_s => ref_elems(ref_ID_s)
        self%basis_c => ref_elems(ref_ID_c)



        !
        ! (Re)Allocate storage for element data structures
        !
        if (allocated(self%jinv)) &
            deallocate( self%jinv,                      &
                        self%metric,                    &
                        self%quad_pts,                  &
                        self%ale_quad_pts,              &
                        self%grad1,                     &
                        self%grad2,                     &
                        self%grad3,                     &
                        self%grad1_trans,               &
                        self%grad2_trans,               &
                        self%grad3_trans,               &
                        self%mass,                      &
                        self%invmass,                   &
                        self%mass_c,                    &
                        self%invmass_c,                 &
                        self%jinv_ale,                  &
                        self%grid_vel,                  &
                        self%jacobian_grid,             &
                        self%inv_jacobian_grid,         &
                        self%det_jacobian_grid,         &
                        self%det_jacobian_grid_grad1,   &
                        self%det_jacobian_grid_grad2,   &
                        self%det_jacobian_grid_grad3,   &
                        self%det_jacobian_grid_modes,   &
                        self%dtau                       &
                        )
            

        nnodes = ref_elems(ref_ID_s)%nnodes_ie()
        allocate(self%jinv(nnodes),                             &
                 self%metric(3,3,nnodes),                       &
                 self%quad_pts(nnodes,3),                       &
                 self%ale_quad_pts(nnodes,3),                   &
                 self%grad1(nnodes,nterms_s),                   &
                 self%grad2(nnodes,nterms_s),                   &
                 self%grad3(nnodes,nterms_s),                   &
                 self%grad1_trans(nterms_s,nnodes),             &
                 self%grad2_trans(nterms_s,nnodes),             &
                 self%grad3_trans(nterms_s,nnodes),             &
                 self%mass(nterms_s,nterms_s),                  &
                 self%invmass(nterms_s,nterms_s),               &
                 self%mass_c(self%nterms_c,self%nterms_c),      &
                 self%invmass_c(self%nterms_c,self%nterms_c),   &
                 self%jinv_ale(nnodes),                         &
                 self%grid_vel(nnodes,3),                       &
                 self%jacobian_grid(nnodes,3,3),                &
                 self%inv_jacobian_grid(nnodes,3,3),            &
                 self%det_jacobian_grid(nnodes),                &
                 self%det_jacobian_grid_grad1(nnodes),          &
                 self%det_jacobian_grid_grad2(nnodes),          &
                 self%det_jacobian_grid_grad3(nnodes),          &
                 self%det_jacobian_grid_modes(self%nterms_s),   &
                 self%dtau(nfields), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Call element metric and matrix calculation routines
        !
        call self%compute_element_matrices()
        call self%update_element_ale()

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
    !*****************************************************************************************







!    !>  Assign quadrature instances for solution modes (GQ) and mesh modes (GQMESH)
!    !!      self%gq
!    !!      self%gqmesh
!    !!
!    !!  TODO: would be good to eliminate pointers in the element data type and just 
!    !!        use integer indices to a global array of quadrature instances.
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   2/1/2016
!    !!
!    !!
!    !-----------------------------------------------------------------------------------------
!    subroutine assign_quadrature(self)
!        use mod_quadrature,     only: compute_nnodes_gq
!        class(element_t),   intent(inout)   :: self
!
!        character(:), allocatable   :: user_msg
!        integer(ik)                 :: nnodes_face, nnodes_vol, igq_s, igq_f
!
!
!        associate (spacedim => self%spacedim, nterms_s => self%nterms_s, nterms_c => self%nterms_c)
!
!            user_msg = "element%assign_quadrature: coordinate expansion not defined."
!            if (nterms_c == 0) call chidg_signal(FATAL,user_msg)
!
!            !
!            ! Get number of quadrature nodes
!            !
!            call compute_nnodes_gq(spacedim,nterms_s,nterms_c,nnodes_face,nnodes_vol)
!
!
!            !
!            ! Get solution quadrature instance
!            !
!            call get_quadrature(spacedim,nterms_s,nnodes_vol,nnodes_face,igq_s)
!            self%gq => GQ(igq_s)
!
!
!            !
!            ! Get coordinate quadrature instance
!            !
!            call get_quadrature(spacedim,nterms_c,nnodes_vol,nnodes_face,igq_f)
!            self%gqmesh => GQ(igq_f)
!
!
!        end associate
!
!    end subroutine assign_quadrature
!    !*****************************************************************************************









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
        !call self%compute_mass_matrix_c()

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

        integer(ik)                 :: inode, nnodes, ierr
        character(:),   allocatable :: coordinate_system

        real(rk),   dimension(:),   allocatable             ::  &
            d1dxi, d1deta, d1dzeta,                             &
            d2dxi, d2deta, d2dzeta,                             &
            d3dxi, d3deta, d3dzeta,                             &
            scaling_12, scaling_13, scaling_23, scaling_123, weights

        real(rk),   dimension(:,:), allocatable :: val, ddxi, ddeta, ddzeta

        nnodes  = self%basis_c%nnodes_ie()
        weights = self%basis_c%weights()
        val     = self%basis_c%interpolator('Value')
        ddxi    = self%basis_c%interpolator('ddxi')
        ddeta   = self%basis_c%interpolator('ddeta')
        ddzeta  = self%basis_c%interpolator('ddzeta')

        !
        ! Compute element metric terms
        !
        d1dxi   = matmul(ddxi,   self%coords%getvar(1,itime = 1))
        d1deta  = matmul(ddeta,  self%coords%getvar(1,itime = 1))
        d1dzeta = matmul(ddzeta, self%coords%getvar(1,itime = 1))

        d2dxi   = matmul(ddxi,   self%coords%getvar(2,itime = 1))
        d2deta  = matmul(ddeta,  self%coords%getvar(2,itime = 1))
        d2dzeta = matmul(ddzeta, self%coords%getvar(2,itime = 1))

        d3dxi   = matmul(ddxi,   self%coords%getvar(3,itime = 1))
        d3deta  = matmul(ddeta,  self%coords%getvar(3,itime = 1))
        d3dzeta = matmul(ddzeta, self%coords%getvar(3,itime = 1))



        !
        ! Define area/volume scaling for coordinate system
        !   Cartesian:
        !       12 = x-y  ;  13 = x-z  ;  23 = y-z
        !
        !   Cylindrical
        !       12 = r-theta  ;  13 = r-z      ;  23 = theta-z
        !
        allocate(scaling_12(nnodes), scaling_13(nnodes), scaling_23(nnodes), scaling_123(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError

        !select case (trim(self%coordinate_system))
        select case (self%coordinate_system)
            case (CARTESIAN)
                scaling_12  = ONE
                scaling_13  = ONE
                scaling_23  = ONE
                scaling_123 = ONE
            case (CYLINDRICAL)
                scaling_12  = self%quad_pts(:,1)
                scaling_13  = ONE
                scaling_23  = self%quad_pts(:,1)
                scaling_123 = self%quad_pts(:,1)
            case default
                call chidg_signal_one(FATAL,"element%compute_quadrature_metrics: Invalid coordinate system. Choose 'Cartesian' or 'Cylindrical'.",self%coordinate_system)
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
        self%vol = abs(sum(self%jinv * weights))


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
        integer(ik)                         :: iterm,inode

        integer(ik)                                 :: nnodes
        real(rk),   allocatable,    dimension(:,:)  :: ddxi, ddeta, ddzeta

        nnodes = self%basis_s%nnodes_ie()
        ddxi   = self%basis_s%interpolator('ddxi')
        ddeta  = self%basis_s%interpolator('ddeta')
        ddzeta = self%basis_s%interpolator('ddzeta')

        do iterm = 1,self%nterms_s
            do inode = 1,nnodes
                self%grad1(inode,iterm) = self%metric(1,1,inode) * ddxi(inode,iterm) + &
                                          self%metric(2,1,inode) * ddeta(inode,iterm) + &
                                          self%metric(3,1,inode) * ddzeta(inode,iterm)

                self%grad2(inode,iterm) = self%metric(1,2,inode) * ddxi(inode,iterm) + &
                                          self%metric(2,2,inode) * ddeta(inode,iterm) + &
                                          self%metric(3,2,inode) * ddzeta(inode,iterm)

                self%grad3(inode,iterm) = self%metric(1,3,inode) * ddxi(inode,iterm) + &
                                          self%metric(2,3,inode) * ddeta(inode,iterm) + &
                                          self%metric(3,3,inode) * ddzeta(inode,iterm)
            end do
        end do


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

        integer(ik)                             :: nnodes
        real(rk),   dimension(:),   allocatable :: coord1, coord2, coord3
        integer(ik)                             :: inode

        nnodes = self%basis_c%nnodes_ie()

        !
        ! compute coordinates associated with quadrature points
        !
        coord1 = matmul(self%basis_c%interpolator('Value'),self%coords%getvar(1,itime = 1))
        coord2 = matmul(self%basis_c%interpolator('Value'),self%coords%getvar(2,itime = 1))
        coord3 = matmul(self%basis_c%interpolator('Value'),self%coords%getvar(3,itime = 1))


        !
        ! Initialize each point with coordinates
        !
        do inode = 1,nnodes
            self%quad_pts(inode,1:3) = [coord1(inode), coord2(inode), coord3(inode)]
        end do

    end subroutine compute_quadrature_coords
    !*****************************************************************************************









!    !>  Compute element volume.
!    !!
!    !!  Approach: 
!    !!
!    !!  @author Nathan A. Wukie (AFRL)
!    !!  @date   7/5/2017
!    !!
!    !!
!    !-----------------------------------------------------------------------------------------
!    subroutine compute_volume(self,ref_ID)
!        class(element_t),   intent(inout)   :: self
!        integer(ik),        intent(in)      :: ref_ID
!
!
!        integer(ik)                 :: inode, nnodes, ierr
!        character(:),   allocatable :: coordinate_system
!
!        real(rk),   dimension(:),   allocatable             ::  &
!            d1dxi, d1deta, d1dzeta,                             &
!            d2dxi, d2deta, d2dzeta,                             &
!            d3dxi, d3deta, d3dzeta,                             &
!            scaling_12, scaling_13, scaling_23, scaling_123, weights
!
!        real(rk),   dimension(:,:), allocatable :: ddxi, ddeta, ddzeta
!
!        nnodes  = self%basis_c%nnodes_ie()
!        weights = self%basis_c%weights()
!        ddxi    = self%basis_c%interpolator('ddxi')
!        ddeta   = self%basis_c%interpolator('ddeta')
!        ddzeta  = self%basis_c%interpolator('ddzeta')
!
!        !
!        ! Compute element metric terms
!        !
!        d1dxi   = matmul(ddxi,   self%coords%getvar(1,itime = 1))
!        d1deta  = matmul(ddeta,  self%coords%getvar(1,itime = 1))
!        d1dzeta = matmul(ddzeta, self%coords%getvar(1,itime = 1))
!
!        d2dxi   = matmul(ddxi,   self%coords%getvar(2,itime = 1))
!        d2deta  = matmul(ddeta,  self%coords%getvar(2,itime = 1))
!        d2dzeta = matmul(ddzeta, self%coords%getvar(2,itime = 1))
!
!        d3dxi   = matmul(ddxi,   self%coords%getvar(3,itime = 1))
!        d3deta  = matmul(ddeta,  self%coords%getvar(3,itime = 1))
!        d3dzeta = matmul(ddzeta, self%coords%getvar(3,itime = 1))
!
!
!
!        !
!        ! Define area/volume scaling for coordinate system
!        !   Cartesian:
!        !       12 = x-y  ;  13 = x-z  ;  23 = y-z
!        !
!        !   Cylindrical
!        !       12 = r-theta  ;  13 = r-z      ;  23 = theta-z
!        !
!        allocate(scaling_12(nnodes), scaling_13(nnodes), scaling_23(nnodes), scaling_123(nnodes), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!        select case (self%coordinate_system)
!            case ('Cartesian')
!                scaling_12  = ONE
!                scaling_13  = ONE
!                scaling_23  = ONE
!                scaling_123 = ONE
!            case ('Cylindrical')
!                scaling_12  = self%quad_pts(:,1)
!                scaling_13  = ONE
!                scaling_23  = self%quad_pts(:,1)
!                scaling_123 = self%quad_pts(:,1)
!            case default
!                call chidg_signal_one(FATAL,"element%compute_quadrature_metrics: Invalid coordinate system. Choose 'Cartesian' or 'Cylindrical'.",self%coordinate_system)
!        end select
!
!
!        !
!        ! Compute inverse cell mapping jacobian
!        !
!        self%jinv = scaling_123*(d1dxi*d2deta*d3dzeta  -  d1deta*d2dxi*d3dzeta - &
!                                 d1dxi*d2dzeta*d3deta  +  d1dzeta*d2dxi*d3deta + &
!                                 d1deta*d2dzeta*d3dxi  -  d1dzeta*d2deta*d3dxi)
!
!        !
!        ! Check for negative jacobians
!        !
!        if (any(self%jinv < ZERO)) call chidg_signal(FATAL,"element%compute_quadrature_metrics: Negative element jacobians. Check element quality and orientation.")
!
!
!        !
!        ! Compute element volume
!        !
!        self%vol = abs(sum(self%jinv * weights))
!
!
!
!    end subroutine compute_volume
!    !*****************************************************************************************









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

        temp = transpose(self%basis_s%interpolator('Value'))


        !
        ! Multiply rows by quadrature weights and cell jacobians
        !
        do iterm = 1,self%nterms_s
            temp(iterm,:) = temp(iterm,:)*(self%basis_s%weights())*(self%jinv)
        end do


        !
        ! Perform the matrix multiplication of the transpose val matrix by
        ! the standard matrix. This produces the mass matrix. I think...
        !
        self%mass = matmul(temp,self%basis_s%interpolator('Value'))


        !
        ! Compute and store the inverted mass matrix
        !
        self%invmass = inv(self%mass)



    end subroutine compute_mass_matrix
    !******************************************************************************************





    !>  Compute element-local mass matrix for coordinate expansion
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_mass_matrix_c(self)
        class(element_t), intent(inout) :: self

        integer(ik)             :: iterm
        real(rk),   allocatable :: temp(:,:)

        self%invmass_c = ZERO
        self%mass_c    = ZERO

        temp = transpose(self%basis_c%interpolator('Value'))


        !
        ! Multiply rows by quadrature weights and cell jacobians
        !
        do iterm = 1,self%nterms_c
            temp(iterm,:) = temp(iterm,:)*(self%basis_c%weights())*(self%jinv)
        end do


        !
        ! Perform the matrix multiplication of the transpose val matrix by
        ! the standard matrix. This produces the mass matrix. I think...
        !
        self%mass_c = matmul(temp,self%basis_c%interpolator('Value'))


        !
        ! Compute and store the inverted mass matrix
        !
        if (self%ielement_g==1) then
        do iterm = 1, self%nterms_c
        print *, self%mass_c(iterm,iterm)
        end do
        print *, self%nterms_c
        print *, size(self%jinv)
        print *, size(self%basis_c%weights())
        print *, size(temp,1)
        print *, size(temp,2)
        end if
        self%invmass_c = inv(self%mass_c)



    end subroutine compute_mass_matrix_c
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

        real(rk)                :: val1, val2, val3
        real(rk)                :: phys_point(3)
        real(rk)                :: polyvals(self%nterms_c)
        integer(ik)             :: iterm, spacedim



        !
        ! Evaluate polynomial modes at node location
        !
        spacedim = self%spacedim
        do iterm = 1,self%nterms_c
            polyvals(iterm)  = polynomial_val(spacedim,self%nterms_c,iterm,[xi,eta,zeta])
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
        !call phys_point%set(val1,val2,val3) 
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
    function grid_point(self,type,icoord,xi,eta,zeta) result(val)
        class(element_t),   intent(in)  :: self
        character(*),       intent(in)  :: type
        integer(ik),        intent(in)  :: icoord
        real(rk),           intent(in)  :: xi, eta, zeta

        real(rk)        :: val
        real(rk)        :: polyvals(self%nterms_c)
        integer(ik)     :: iterm, spacedim

        if (icoord > 3)                  call chidg_signal(FATAL,"element%grid_point -- icoord exceeded 3 physical coordinates")
        if (.not. self%geom_initialized) call chidg_signal(FATAL,"element%grid_point: geometry not initialized")


        ! Evaluate polynomial modes at node location
        spacedim = self%spacedim
        do iterm = 1,self%nterms_c
                polyvals(iterm) = polynomial_val(spacedim,self%nterms_c,iterm,[xi,eta,zeta])
        end do


        ! Evaluate mesh point from dot product of modes and polynomial values
        if (type == 'Reference') then
            val = dot_product(self%coords%getvar(icoord,itime = 1), polyvals)
        else if (type == 'ALE') then
            val = dot_product(self%ale_coords%getvar(icoord,itime = 1), polyvals)
        end if


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
        real(rk)        :: polyvals(self%nterms_c)
        integer(ik)     :: iterm, spacedim


        if (phys_dir > 3) call chidg_signal(FATAL,"element%metric_point: phys_dir exceeded 3 physical coordinates")
        if (comp_dir > 3) call chidg_signal(FATAL,"element%metric_point: comp_dir exceeded 3 physical coordinates")



        !
        ! Evaluate polynomial modes at node location
        !
        spacedim = self%spacedim
        do iterm = 1,self%nterms_c
            polyvals(iterm) = dpolynomial_val(spacedim,self%nterms_c,iterm,[xi,eta,zeta],comp_dir)
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
                if (self%coordinate_system == CARTESIAN) then

                else if (self%coordinate_system == CYLINDRICAL) then
                    if (phys_dir == DIR_THETA) then
                        r = self%grid_point('Reference',1,xi,eta,zeta)
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

        real(rk)                   :: temp,val
        real(rk)                   :: polyvals(q%nterms())
        integer(ik)                :: iterm, spacedim


        ! evaluate polynomial modes at node location
        spacedim = self%spacedim
        do iterm = 1,q%nterms()
            polyvals(iterm)  = polynomial_val(spacedim,q%nterms(),iterm,[xi,eta,zeta])
        end do


        ! evaluate x from dot product of modes and polynomial values
        temp = dot_product(q%getvar(ivar,itime),polyvals)
        val = temp/dot_product(self%det_jacobian_grid_modes,polyvals)


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
        real(rk)        :: ddxi(q%nterms()), ddeta(q%nterms()), ddzeta(q%nterms()), &
                           deriv(q%nterms())
        real(rk)        :: metric(3,3), jinv, dxi_dx, dxi_dy, dxi_dz, &
                           deta_dx, deta_dy, deta_dz, dzeta_dx, dzeta_dy, dzeta_dz
        integer(ik)     :: iterm, spacedim


        !
        ! Evaluate polynomial mode derivatives at node location
        !
        spacedim = self%spacedim
        do iterm = 1,q%nterms()
            ddxi(iterm)   = dpolynomial_val(spacedim,q%nterms(),iterm,[xi,eta,zeta],XI_DIR)
            ddeta(iterm)  = dpolynomial_val(spacedim,q%nterms(),iterm,[xi,eta,zeta],ETA_DIR)
            ddzeta(iterm) = dpolynomial_val(spacedim,q%nterms(),iterm,[xi,eta,zeta],ZETA_DIR)
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
    !!  NOTE: Will return a location, even if the newton solve did not converge. If the
    !!        newton solve did not converge, the point values will be NaN, so the result
    !!        should be checked for NaN's before being used. The recommended approach is:
    !!       
    !!      use ieee_arithmetic,    only: ieee_is_nan
    !!      valid_point = (any(ieee_is_nan(result)))
    !!        
    !!
    !!      
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

        !type(point_t)       :: loc, point_n
        real(rk)            :: loc(3), point_n(3)
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
            R(1) = -(point_n(1) - coord1)
            R(2) = -(point_n(2) - coord2)
            R(3) = -(point_n(3) - coord3)


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
                !loc%status = VALID_POINT  ! point found
                !call loc%set(xi,eta,zeta)
                loc = [xi, eta, zeta]
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
                !loc%status = INVALID_POINT  ! point not found
                loc = ieee_value(1._rk,ieee_quiet_nan)
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


        ! Call function for evaluation at quadrature nodes and multiply by quadrature weights
        time  = 0._rk
        fvals = fcn%compute(time,point_t(self%quad_pts)) * self%basis_s%weights() * self%jinv


        ! Project
        temp = matmul(transpose(self%basis_s%interpolator('Value')),fvals)
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





    !>
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------
    subroutine update_element_ale(self)
        class(element_t),       intent(inout)      :: self

        call self%update_geom_ale()
        call self%compute_quadrature_coords_ale()
        call self%compute_quadrature_metrics_ale()

    end subroutine update_element_ale
    !******************************************************************************



    !>
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------
    subroutine update_geom_ale(self)
        class(element_t),               intent(inout)   :: self

        real(rk),   allocatable :: xmodes(:), ymodes(:), zmodes(:)

        ! We are assuming that the ale_elem_pts have already been updated according to mesh motion
        ! A worker has a pointer to the mesh, which it has used to update this field


        !
        ! Compute mesh x,y,z modes
        !
        xmodes = matmul(self%nodes_to_modes,self%ale_elem_pts(:,1))
        ymodes = matmul(self%nodes_to_modes,self%ale_elem_pts(:,2))
        zmodes = matmul(self%nodes_to_modes,self%ale_elem_pts(:,3))

        call self%ale_coords%setvar(1,itime = 1,vals = xmodes)
        call self%ale_coords%setvar(2,itime = 1,vals = ymodes)
        call self%ale_coords%setvar(3,itime = 1,vals = zmodes)



        !
        ! Compute grid velocity modes
        !
        xmodes = matmul(self%nodes_to_modes,self%ale_vel_elem_pts(:,1))
        ymodes = matmul(self%nodes_to_modes,self%ale_vel_elem_pts(:,2))
        zmodes = matmul(self%nodes_to_modes,self%ale_vel_elem_pts(:,3))

        call self%ale_vel_coords%setvar(1,itime = 1,vals = xmodes)
        call self%ale_vel_coords%setvar(2,itime = 1,vals = ymodes)
        call self%ale_vel_coords%setvar(3,itime = 1,vals = zmodes)


    end subroutine update_geom_ale
    !***********************************************************************************************************
    




    !>
    !!
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine compute_quadrature_coords_ale(self)
        class(element_t),   intent(inout)   :: self

        real(rk),   allocatable :: val(:,:)


        !
        ! Retrieve interpolator
        !
        val = self%basis_c%interpolator('Value')

        !
        ! compute cartesian coordinates associated with quadrature points
        !
        self%ale_quad_pts(:,1) = matmul(val,self%ale_coords%getvar(1,itime = 1))
        self%ale_quad_pts(:,2) = matmul(val,self%ale_coords%getvar(2,itime = 1))
        self%ale_quad_pts(:,3) = matmul(val,self%ale_coords%getvar(3,itime = 1))


        !
        ! Grid velocity
        ! compute cartesian coordinates associated with quadrature points
        !
        self%grid_vel(:,1) = matmul(val,self%ale_vel_coords%getvar(1,itime = 1))
        self%grid_vel(:,2) = matmul(val,self%ale_vel_coords%getvar(2,itime = 1))
        self%grid_vel(:,3) = matmul(val,self%ale_vel_coords%getvar(3,itime = 1))

    end subroutine compute_quadrature_coords_ale
    !****************************************************************************************



    !>
    !!
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine compute_quadrature_metrics_ale(self)
        class(element_t),    intent(inout)   :: self

        integer(ik)             :: inode
        integer(ik)             :: nnodes
        character(:),   allocatable :: coordinate_system

        integer(ik)                             :: ierr
        real(rk),   dimension(:),   allocatable ::  &
            d1dxi, d1deta, d1dzeta,                 &
            d2dxi, d2deta, d2dzeta,                 &
            d3dxi, d3deta, d3dzeta,                 &
            d1dxi_ale, d1deta_ale, d1dzeta_ale,                 &
            d2dxi_ale, d2deta_ale, d2dzeta_ale,                 &
            d3dxi_ale, d3deta_ale, d3dzeta_ale,                 &
            scaling_12, scaling_13, scaling_23, scaling_123,    &
            fvals, temp, weights

        real(rk),   dimension(:,:), allocatable :: val, ddxi, ddeta, ddzeta, tempmat
        real(rk), dimension(:,:,:), allocatable :: jacobian_matrix, jacobian_matrix_ale

        !
        ! Retrieve interpolators
        !
        nnodes  = self%basis_c%nnodes_ie()
        weights = self%basis_c%weights()
        ddxi    = self%basis_c%interpolator('ddxi')
        ddeta   = self%basis_c%interpolator('ddeta')
        ddzeta  = self%basis_c%interpolator('ddzeta')

        d1dxi   = matmul(ddxi,  self%coords%getvar(1,itime = 1))
        d1deta  = matmul(ddeta, self%coords%getvar(1,itime = 1))
        d1dzeta = matmul(ddzeta,self%coords%getvar(1,itime = 1))

        d2dxi   = matmul(ddxi,  self%coords%getvar(2,itime = 1))
        d2deta  = matmul(ddeta, self%coords%getvar(2,itime = 1))
        d2dzeta = matmul(ddzeta,self%coords%getvar(2,itime = 1))

        d3dxi   = matmul(ddxi,  self%coords%getvar(3,itime = 1))
        d3deta  = matmul(ddeta, self%coords%getvar(3,itime = 1))
        d3dzeta = matmul(ddzeta,self%coords%getvar(3,itime = 1))


        d1dxi_ale   = matmul(ddxi,  self%ale_coords%getvar(1,itime = 1))
        d1deta_ale  = matmul(ddeta, self%ale_coords%getvar(1,itime = 1))
        d1dzeta_ale = matmul(ddzeta,self%ale_coords%getvar(1,itime = 1))

        d2dxi_ale   = matmul(ddxi,  self%ale_coords%getvar(2,itime = 1))
        d2deta_ale  = matmul(ddeta, self%ale_coords%getvar(2,itime = 1))
        d2dzeta_ale = matmul(ddzeta,self%ale_coords%getvar(2,itime = 1))

        d3dxi_ale   = matmul(ddxi,  self%ale_coords%getvar(3,itime = 1))
        d3deta_ale  = matmul(ddeta, self%ale_coords%getvar(3,itime = 1))
        d3dzeta_ale = matmul(ddzeta,self%ale_coords%getvar(3,itime = 1))

        !
        ! Define area/volume scaling for coordinate system
        !   Cartesian:
        !       12 = x-y  ;  13 = x-z  ;  23 = y-z
        !
        !   Cylindrical
        !       12 = r-theta  ;  13 = r-z      ;  23 = theta-z
        !
        allocate(scaling_12(nnodes), scaling_13(nnodes), scaling_23(nnodes), scaling_123(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError

        select case (self%coordinate_system)
            case (CARTESIAN)
                scaling_12  = ONE
                scaling_13  = ONE
                scaling_23  = ONE
                scaling_123 = ONE
            case (CYLINDRICAL)
                scaling_12  = self%quad_pts(:,1)
                scaling_13  = ONE
                scaling_23  = self%quad_pts(:,1)
                scaling_123 = self%quad_pts(:,1)
            case default
                call chidg_signal(FATAL,"element%compute_quadrature_metrics_ale: Invalid coordinate system. Choose 'Cartesian' or 'Cylindrical'.")
        end select


        !
        ! Compute inverse cell mapping jacobian
        !
        self%jinv_ale = scaling_123*(d1dxi_ale*d2deta_ale*d3dzeta_ale  -  d1deta_ale*d2dxi_ale*d3dzeta_ale - &
                                     d1dxi_ale*d2dzeta_ale*d3deta_ale  +  d1dzeta_ale*d2dxi_ale*d3deta_ale + &
                                     d1deta_ale*d2dzeta_ale*d3dxi_ale  -  d1dzeta_ale*d2deta_ale*d3dxi_ale)

        !
        ! Check for negative jacobians
        !
        if (any(self%jinv_ale < ZERO)) call chidg_signal(FATAL,"element%compute_quadrature_metrics_ale: Negative element jacobians. Check element quality and orientation.")


        !
        ! Compute element volume
        !
        self%vol_ale = abs(sum(self%jinv_ale * weights))


        allocate(jacobian_matrix(nnodes,3,3), jacobian_matrix_ale(nnodes,3,3), tempmat(3,3))
        do inode = 1,nnodes
            jacobian_matrix(inode,1,1) = d1dxi(inode)
            jacobian_matrix(inode,1,2) = d1deta(inode)
            jacobian_matrix(inode,1,3) = d1dzeta(inode)
                                     
            jacobian_matrix(inode,2,1) = d2dxi(inode)
            jacobian_matrix(inode,2,2) = d2deta(inode)
            jacobian_matrix(inode,2,3) = d2dzeta(inode)
                                     
            jacobian_matrix(inode,3,1) = d3dxi(inode)
            jacobian_matrix(inode,3,2) = d3deta(inode)
            jacobian_matrix(inode,3,3) = d3dzeta(inode)


            tempmat = inv(jacobian_matrix(inode,:,:))
            jacobian_matrix(inode,:,:) = tempmat 

            jacobian_matrix_ale(inode,1,1) = d1dxi_ale(inode)
            jacobian_matrix_ale(inode,1,2) = d1deta_ale(inode)
            jacobian_matrix_ale(inode,1,3) = d1dzeta_ale(inode)
                                         
            jacobian_matrix_ale(inode,2,1) = d2dxi_ale(inode)
            jacobian_matrix_ale(inode,2,2) = d2deta_ale(inode)
            jacobian_matrix_ale(inode,2,3) = d2dzeta_ale(inode)
                                         
            jacobian_matrix_ale(inode,3,1) = d3dxi_ale(inode)
            jacobian_matrix_ale(inode,3,2) = d3deta_ale(inode)
            jacobian_matrix_ale(inode,3,3) = d3dzeta_ale(inode)


            self%jacobian_grid(inode,:,:) = matmul(jacobian_matrix_ale(inode,:,:),jacobian_matrix(inode,:,:))
            self%inv_jacobian_grid(inode,:,:) = inv(self%jacobian_grid(inode,:,:))
        end do



        self%det_jacobian_grid = self%jinv_ale/self%jinv
        fvals = self%det_jacobian_grid * weights * self%jinv


        !
        ! Project
        !
        val  = self%basis_s%interpolator('Value')
        temp = matmul(transpose(val),fvals)
        self%det_jacobian_grid_modes = matmul(self%invmass,temp)


        self%det_jacobian_grid_grad1 = matmul(self%grad1,self%det_jacobian_grid_modes)
        self%det_jacobian_grid_grad2 = matmul(self%grad2,self%det_jacobian_grid_modes)
        self%det_jacobian_grid_grad3 = matmul(self%grad3,self%det_jacobian_grid_modes)


    end subroutine compute_quadrature_metrics_ale
    !********************************************************************************************************









    !> Compute ALE coordinate metric term at a given point in computational space
    !!
    !!  @author Eric Wolf
    !!  @date   7/21/2017
    !!
    !!  @param[in]  elem        element_t containing the geometry definition and data
    !!  @param[in]  phys_dir    physical coordinate being differentiated
    !!  @param[in]  comp_dir    Computational coordinate being differentiated with respect to
    !!  @param[in]  xi          Computational coordinate - xi
    !!  @param[in]  eta         Computational coordinate - eta
    !!  @param[in]  zeta        Computational coordinate - zeta
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function metric_point_ale(self,phys_dir,comp_dir,xi,eta,zeta,scale) result(val)
        class(element_t),   intent(in)              :: self
        integer(ik),        intent(in)              :: phys_dir
        integer(ik),        intent(in)              :: comp_dir
        real(rk),           intent(in)              :: xi, eta, zeta
        logical,            intent(in), optional    :: scale
        
        real(rk)        :: val, r
        real(rk)        :: polyvals(self%nterms_c)
        integer(ik)     :: iterm, spacedim


        if (phys_dir > 3) call chidg_signal(FATAL,"element%metric_point: phys_dir exceeded 3 physical coordinates")
        if (comp_dir > 3) call chidg_signal(FATAL,"element%metric_point: comp_dir exceeded 3 physical coordinates")



        !
        ! Evaluate polynomial modes at node location
        !
        spacedim = self%spacedim
        do iterm = 1,self%nterms_c
            polyvals(iterm) = dpolynomial_val(spacedim,self%nterms_c,iterm,[xi,eta,zeta],comp_dir)
        end do


        !
        ! Evaluate mesh point from dot product of modes and polynomial values
        !
        val = dot_product(self%ale_coords%getvar(phys_dir, itime=1), polyvals)



        !
        ! Apply scaling due to coordinate system.
        !
        if (present(scale)) then
            if (scale) then
                if (self%coordinate_system == CARTESIAN) then

                else if (self%coordinate_system == CYLINDRICAL) then
                    if (phys_dir == DIR_THETA) then
                        r = self%grid_point('Reference',1,xi,eta,zeta)
                        val = val * r
                    end if
                end if
            end if
        end if



    end function metric_point_ale
    !*****************************************************************************************




    !>  Compute ALE grid quantities, based on the location in reference space (xi, eta, zeta)
    !!
    !!  @author Eric Wolf
    !!  @date   7/21/2017
    !!
    !!  @param[in]  elem    Element that the solution expansion is associated with.
    !!  @param[in]  xi      Real value for xi-coordinate.
    !!  @param[in]  eta     Real value for eta-coordinate.
    !!  @param[in]  zeta    Real value for zeta-coordinate.
    !!  @param[inout]  det_jacobian_grid       Determinant of grid Jacobian
    !!  @param[inout]  inv_jacobian_grid       Inverse of grid Jacobian
    !!  @param[inout]  grid_vel                Grid velocities
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine ale_point(self,xi,eta,zeta,det_jacobian_grid,det_jacobian_grid_grad,inv_jacobian_grid,grid_vel) 
        class(element_t),       intent(in)      :: self
        real(rk),               intent(in)      :: xi,eta,zeta
        real(rk),               intent(inout)   :: det_jacobian_grid
        real(rk),               intent(inout)   :: det_jacobian_grid_grad(3)
        real(rk),               intent(inout)   :: inv_jacobian_grid(3,3)
        real(rk),               intent(inout)   :: grid_vel(3)

        real(rk)        :: metric(3,3), jinv, metric_ale(3,3), jinv_ale, polyval(self%nterms_s)
        real(rk), dimension(self%nterms_s)  :: ddxi, ddeta, ddzeta, grad1, grad2, grad3
        integer(ik)     :: iterm, spacedim, itime


        !
        ! Evaluate polynomial mode derivatives at node location
        !
        spacedim = self%spacedim
        
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

        metric_ale(1,1) = self%metric_point_ale(DIR_1,XI_DIR,  xi,eta,zeta)
        metric_ale(2,1) = self%metric_point_ale(DIR_2,XI_DIR,  xi,eta,zeta)
        metric_ale(3,1) = self%metric_point_ale(DIR_3,XI_DIR,  xi,eta,zeta)
        metric_ale(1,2) = self%metric_point_ale(DIR_1,ETA_DIR, xi,eta,zeta)
        metric_ale(2,2) = self%metric_point_ale(DIR_2,ETA_DIR, xi,eta,zeta)
        metric_ale(3,2) = self%metric_point_ale(DIR_3,ETA_DIR, xi,eta,zeta)
        metric_ale(1,3) = self%metric_point_ale(DIR_1,ZETA_DIR,xi,eta,zeta)
        metric_ale(2,3) = self%metric_point_ale(DIR_2,ZETA_DIR,xi,eta,zeta)
        metric_ale(3,3) = self%metric_point_ale(DIR_3,ZETA_DIR,xi,eta,zeta)

        !
        ! Compute inverse cell mapping jacobian
        !
        jinv = metric(1,1)*metric(2,2)*metric(3,3) - metric(1,2)*metric(2,1)*metric(3,3) - &
               metric(1,1)*metric(2,3)*metric(3,2) + metric(1,3)*metric(2,1)*metric(3,2) + &
               metric(1,2)*metric(2,3)*metric(3,1) - metric(1,3)*metric(2,2)*metric(3,1)


        jinv_ale = metric_ale(1,1)*metric_ale(2,2)*metric_ale(3,3) - metric_ale(1,2)*metric_ale(2,1)*metric_ale(3,3) - &
                   metric_ale(1,1)*metric_ale(2,3)*metric_ale(3,2) + metric_ale(1,3)*metric_ale(2,1)*metric_ale(3,2) + &
                   metric_ale(1,2)*metric_ale(2,3)*metric_ale(3,1) - metric_ale(1,3)*metric_ale(2,2)*metric_ale(3,1)


        det_jacobian_grid = jinv_ale/jinv
        inv_jacobian_grid = matmul(inv(metric_ale),metric)


        ! evaluate polynomial modes at node location
        do iterm = 1,self%nterms_s
            polyval(iterm)  = polynomial_val(spacedim,self%nterms_s,iterm,[xi,eta,zeta])
        end do


        ! evaluate grid velocities from dot product of modes and polynomial values
        grid_vel(1) = dot_product(polyval,self%ale_vel_coords%getvar(1,itime = 1))
        grid_vel(2) = dot_product(polyval,self%ale_vel_coords%getvar(2,itime = 1))
        grid_vel(3) = dot_product(polyval,self%ale_vel_coords%getvar(3,itime = 1))

        do iterm = 1,self%nterms_s
            ddxi(iterm)   = dpolynomial_val(spacedim,self%nterms_s,iterm,[xi,eta,zeta],XI_DIR)
            ddeta(iterm)  = dpolynomial_val(spacedim,self%nterms_s,iterm,[xi,eta,zeta],ETA_DIR)
            ddzeta(iterm) = dpolynomial_val(spacedim,self%nterms_s,iterm,[xi,eta,zeta],ZETA_DIR)
        end do

        metric = inv(metric)
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

        det_jacobian_grid_grad(1) = dot_product(grad1, self%det_jacobian_grid_modes)
        det_jacobian_grid_grad(2) = dot_product(grad2, self%det_jacobian_grid_modes)
        det_jacobian_grid_grad(3) = dot_product(grad3, self%det_jacobian_grid_modes)


    end subroutine ale_point
    !*****************************************************************************************







    subroutine destructor(self)
        type(element_t), intent(inout) :: self


    end subroutine

end module type_element
