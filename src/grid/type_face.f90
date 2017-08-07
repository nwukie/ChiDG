module type_face
#include  <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX,                 &
                                      ZETA_MIN, ZETA_MAX, XI_DIR, ETA_DIR, ZETA_DIR,    &
                                      NO_INTERIOR_NEIGHBOR, NO_PROC,                    &
                                      ZERO, ONE, TWO, ORPHAN, NO_PMM_ASSIGNED, CARTESIAN, CYLINDRICAL
    use type_reference_element, only: reference_element_t

    use type_point,             only: point_t
    use type_element,           only: element_t
    use type_quadrature,        only: quadrature_t
    use type_densevector,       only: densevector_t
    use mod_inv,                only: inv
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
    !!  @author Mayank Sharma
    !!  @date   11/12/2016
    !!
    !------------------------------------------------------------------------------------------
    type, public :: face_t

        integer(ik)        :: spacedim        ! Number of spatial dimensions

        ! Self information
        integer(ik)        :: ftype           ! INTERIOR, BOUNDARY, CHIMERA, ORPHAN 
        integer(ik)        :: ChiID = 0       ! Identifier for domain-local Chimera interfaces

        integer(ik)        :: bc_ID      = 0  ! Index for bc state group data%bc_state_group(bc_ID)
        integer(ik)        :: group_ID   = 0  ! Index for bc patch group mesh%bc_patch_group(group_ID)
        integer(ik)        :: patch_ID   = 0  ! Index for bc patch 
        integer(ik)        :: face_ID    = 0  ! Index for bc patch face


        integer(ik)        :: pmm_ID = NO_PMM_ASSIGNED

        ! Owner-element information
        integer(ik)        :: face_location(5)! [idomain_g, idomain_l, iparent_g, iparent_l, iface]
        integer(ik)        :: idomain_g       ! Global index of the parent domain
        integer(ik)        :: idomain_l       ! Processor-local index of the parent domain
        integer(ik)        :: iparent_g       ! Domain-global index of the parent element
        integer(ik)        :: iparent_l       ! Processor-local index of the parent element
        integer(ik)        :: iface           ! XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, etc
        integer(ik)        :: neqns           ! Number of equations in equationset_t
        integer(ik)        :: nterms_s        ! Number of terms in solution polynomial expansion
        integer(ik)        :: ntime
        integer(ik)        :: coordinate_system    ! CARTESIAN, CYLINDRICAL. parameters from mod_constants.f90


        ! Neighbor information
        integer(ik)        :: neighbor_location(5) = 0
        integer(ik)        :: ineighbor_proc       = NO_PROC    ! MPI processor rank of the neighboring element
        integer(ik)        :: ineighbor_domain_g   = 0          ! Global index of the neighboring element's domain
        integer(ik)        :: ineighbor_domain_l   = 0          ! Processor-local index of the neighboring element's domain
        integer(ik)        :: ineighbor_element_g  = 0          ! Domain-global index of the neighboring element
        integer(ik)        :: ineighbor_element_l  = 0          ! Processor-local index of the neighboring element
        integer(ik)        :: ineighbor_face       = 0
        integer(ik)        :: ineighbor_neqns      = 0
        integer(ik)        :: ineighbor_nterms_s   = 0
        integer(ik)        :: recv_comm            = 0
        integer(ik)        :: recv_domain          = 0
        integer(ik)        :: recv_element         = 0

        ! Neighbor information: if neighbor is off-processor
        real(rk)                        :: neighbor_h(3)           ! Approximate size of neighbor bounding box
        real(rk),           allocatable :: neighbor_grad1(:,:)     ! Grad of basis functions in at quadrature nodes
        real(rk),           allocatable :: neighbor_grad2(:,:)     ! Grad of basis functions in at quadrature nodes
        real(rk),           allocatable :: neighbor_grad3(:,:)     ! Grad of basis functions in at quadrature nodes
        real(rk),           allocatable :: neighbor_br2_face(:,:)  ! Matrix for computing/obtaining br2 modes at face nodes
        real(rk),           allocatable :: neighbor_br2_vol(:,:)   ! Matrix for computing/obtaining br2 modes at volume nodes
        real(rk),           allocatable :: neighbor_invmass(:,:)    

        ! Neighbor ALE: if neighbor is off-processor
        real(rk), allocatable           :: neighbor_grid_vel(:,:)
        real(rk), allocatable           :: neighbor_inv_jacobian_grid(:,:,:)
        real(rk), allocatable           :: neighbor_det_jacobian_grid(:)
        real(rk), allocatable           :: neighbor_det_jacobian_grid_grad1(:)
        real(rk), allocatable           :: neighbor_det_jacobian_grid_grad2(:)
        real(rk), allocatable           :: neighbor_det_jacobian_grid_grad3(:)


        ! Chimera face offset. For periodic boundary condition.
        logical                         :: periodic_offset  = .false.
        real(rk)                        :: chimera_offset_1 = 0._rk
        real(rk)                        :: chimera_offset_2 = 0._rk
        real(rk)                        :: chimera_offset_3 = 0._rk


        ! Geometry
        type(densevector_t)             :: coords               ! Modal expansion of coordinates 
        real(rk),           allocatable :: quad_pts(:,:)        ! Discrete coordinates at quadrature nodes

        ! Metric terms
        real(rk),           allocatable :: jinv(:)              ! array of inverse element jacobians on the face
        real(rk),           allocatable :: metric(:,:,:)        ! Face metric terms
        real(rk),           allocatable :: norm(:,:)            ! Face normal vector - scaled by differential area
        real(rk),           allocatable :: unorm(:,:)           ! Face normal vector - unit length


        ! Matrices of cartesian gradients of basis/test functions
        real(rk),           allocatable :: grad1(:,:)           ! Deriv of basis functions in at quadrature nodes
        real(rk),           allocatable :: grad2(:,:)           ! Deriv of basis functions in at quadrature nodes
        real(rk),           allocatable :: grad3(:,:)           ! Deriv of basis functions in at quadrature nodes


        ! Quadrature matrices
        type(reference_element_t), pointer  :: basis_s => null()
        type(reference_element_t), pointer  :: basis_c => null()


        ! BR2 matrix
        real(rk),           allocatable :: br2_face(:,:)
        real(rk),           allocatable :: br2_vol(:,:)

        ! Face area
        real(rk)                        :: total_area
        real(rk),           allocatable :: differential_areas(:)

        ! ALE
        real(rk),   allocatable         :: ale_quad_pts(:,:)
        real(rk),   allocatable         :: ale_elem_pts(:,:)
        type(densevector_t)             :: ale_coords                   ! Modal representation of cartesian coordinates (nterms_var,(x,y,z))
        type(densevector_t)             :: ale_vel_coords               ! Modal representation of cartesian coordinates (nterms_var,(x,y,z))
        real(rk), allocatable           :: grid_vel(:,:)
        real(rk), allocatable           :: jacobian_grid(:,:,:)
        real(rk), allocatable           :: inv_jacobian_grid(:,:,:)
        real(rk), allocatable           :: det_jacobian_grid(:)
        real(rk), allocatable           :: det_jacobian_grid_grad1(:)
        real(rk), allocatable           :: det_jacobian_grid_grad2(:)
        real(rk), allocatable           :: det_jacobian_grid_grad3(:)
        real(rk), allocatable           :: det_jacobian_grid_modes(:)

        real(rk), allocatable           :: jinv_ale(:)                      ! jacobian terms at quadrature nodes



        ! Logical tests
        logical :: geom_initialized    = .false.
        logical :: neighborInitialized = .false.
        logical :: numInitialized      = .false.



    contains

        procedure           :: init_geom
        procedure           :: init_face_ale_coords
        procedure           :: init_sol

        procedure           :: init_neighbor

        procedure           :: compute_quadrature_metrics       ! Compute metric terms at quadrature nodes
        procedure           :: compute_quadrature_normals       ! Compute normals at quadrature nodes
        procedure           :: compute_quadrature_coords        ! Compute cartesian coordinates at quadrature nodes
        procedure           :: compute_quadrature_gradients     ! Compute gradients in cartesian coordinates

        
        ! ALE procedures
        procedure, public   :: update_face_ale
        procedure           :: update_face_ale_coords
        procedure           :: update_face_ale_quadrature
        procedure           :: compute_quadrature_coords_ale
        procedure           :: compute_quadrature_metrics_ale

        procedure           :: get_neighbor_element_g           ! Return neighbor element index
        procedure           :: get_neighbor_element_l           ! Return neighbor element index
        procedure           :: get_neighbor_face                ! Return neighbor face index

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
        class(face_t),      intent(inout), target   :: self
        integer(ik),        intent(in)              :: iface
        type(element_t),    intent(in)              :: elem



        !
        ! Set indices
        !
        self%ftype     = ORPHAN
        self%spacedim  = elem%spacedim


        !
        ! Set owner element
        !
        self%idomain_g     = elem%idomain_g
        self%idomain_l     = elem%idomain_l
        self%iparent_g     = elem%ielement_g
        self%iparent_l     = elem%ielement_l
        self%iface         = iface
        self%face_location = [elem%idomain_g, elem%idomain_l, elem%ielement_g, elem%ielement_l, iface]


        !
        ! No neighbor associated at this point
        !
        self%ineighbor_domain_g  = NO_INTERIOR_NEIGHBOR
        self%ineighbor_domain_l  = NO_INTERIOR_NEIGHBOR
        self%ineighbor_element_g = NO_INTERIOR_NEIGHBOR
        self%ineighbor_element_l = NO_INTERIOR_NEIGHBOR
        self%ineighbor_proc      = NO_PROC
        

        !
        ! Set modal representation of element coordinates, displacements, velocities:
        !   1: set reference coordinates
        !   2: set default ALE (displacements, velocities)
        !
        self%coords = elem%coords
        call self%init_face_ale_coords(elem)


        !
        ! Set coordinate system, confirm initialization.
        !
        self%coordinate_system = elem%coordinate_system
        self%geom_initialized  = .true.

    end subroutine init_geom
    !******************************************************************************************









    !>  Initialize ALE data from nodal displacements.
    !!
    !!  @author Eric Wolf (AFRL)
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2017
    !!
    !--------------------------------------------------------------------------------------
    subroutine init_face_ale_coords(self,elem)
        class(face_t),      intent(inout)   :: self
        type(element_t),    intent(in)  :: elem

        self%ale_coords     = elem%ale_coords
        self%ale_vel_coords = elem%ale_vel_coords

    end subroutine init_face_ale_coords 
    !**************************************************************************************



    !>  Initialize ALE data from nodal displacements.
    !!
    !!  @author Eric Wolf (AFRL)
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2017
    !!
    !--------------------------------------------------------------------------------------
    subroutine update_face_ale_coords(self,elem)
        class(face_t),      intent(inout)   :: self
        type(element_t),    intent(in)  :: elem

        call self%init_face_ale_coords(elem)
        self%det_jacobian_grid_modes = elem%det_jacobian_grid_modes

    end subroutine update_face_ale_coords
    !**************************************************************************************






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
        integer(ik),    intent(in)      :: ineighbor_proc
        integer(ik),    intent(in)      :: ineighbor_neqns
        integer(ik),    intent(in)      :: ineighbor_nterms_s


        self%ftype               = ftype
        self%ineighbor_domain_g  = ineighbor_domain_g
        self%ineighbor_domain_l  = ineighbor_domain_l
        self%ineighbor_element_g = ineighbor_element_g
        self%ineighbor_element_l = ineighbor_element_l
        self%ineighbor_face      = ineighbor_face
        self%ineighbor_neqns     = ineighbor_neqns
        self%ineighbor_nterms_s  = ineighbor_nterms_s
        self%ineighbor_proc      = ineighbor_proc

        self%neighbor_location = [ineighbor_domain_g, ineighbor_domain_l, ineighbor_element_g, ineighbor_element_l, ineighbor_face]

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
        real(rk), allocatable   :: tmp(:,:), val_e(:,:), val_f(:,:)

        !
        ! Set indices and associate quadrature instances.
        !
        self%neqns      = elem%neqns
        self%nterms_s   = elem%nterms_s
        self%ntime      = elem%ntime
        self%basis_s    => elem%basis_s
        self%basis_c    => elem%basis_c



        !
        ! (Re)Allocate storage for face data structures.
        !
        if (allocated(self%jinv))                               &
            deallocate(self%jinv,                               &
                       self%quad_pts,                           &
                       self%metric,                             &
                       self%norm,                               &
                       self%unorm,                              &
                       self%ale_quad_pts,                       &
                       self%jinv_ale,                           &
                       self%jacobian_grid,                      &
                       self%inv_jacobian_grid,                  &
                       self%det_jacobian_grid,                  &
                       self%det_jacobian_grid_grad1,            &
                       self%det_jacobian_grid_grad2,            &
                       self%det_jacobian_grid_grad3,            &
                       self%det_jacobian_grid_modes,            &
                       self%grid_vel,                           &
                       self%neighbor_inv_jacobian_grid,         &
                       self%neighbor_det_jacobian_grid,         &
                       self%neighbor_det_jacobian_grid_grad1,   &
                       self%neighbor_det_jacobian_grid_grad2,   &
                       self%neighbor_det_jacobian_grid_grad3,   &
                       self%neighbor_grid_vel,                  &
                       self%grad1,                              &
                       self%grad2,                              &
                       self%grad3                               &
                       ) 



        nnodes = self%basis_s%nnodes_if()
        allocate(self%jinv(nnodes),                                 &
                 self%quad_pts(nnodes,3),                           &
                 self%metric(3,3,nnodes),                           &
                 self%norm(nnodes,3),                               &
                 self%unorm(nnodes,3),                              &
                 self%ale_quad_pts(nnodes,3),                       &
                 self%jinv_ale(nnodes),                             &
                 self%jacobian_grid(nnodes,3,3),                    &
                 self%inv_jacobian_grid(nnodes,3,3),                &
                 self%det_jacobian_grid(nnodes),                    &
                 self%det_jacobian_grid_grad1(nnodes),              &
                 self%det_jacobian_grid_grad2(nnodes),              &
                 self%det_jacobian_grid_grad3(nnodes),              &
                 self%det_jacobian_grid_modes(self%nterms_s),       &
                 self%grid_vel(nnodes,3),                           &
                 self%neighbor_inv_jacobian_grid(nnodes,3,3),       &
                 self%neighbor_det_jacobian_grid(nnodes),           &
                 self%neighbor_det_jacobian_grid_grad1(nnodes),     &
                 self%neighbor_det_jacobian_grid_grad2(nnodes),     &
                 self%neighbor_det_jacobian_grid_grad3(nnodes),     &
                 self%neighbor_grid_vel(nnodes,3),                  &
                 self%grad1(nnodes,self%nterms_s),                  &
                 self%grad2(nnodes,self%nterms_s),                  &
                 self%grad3(nnodes,self%nterms_s),                  &
                 stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Compute metrics, normals, node coordinates
        !
        call self%compute_quadrature_coords()
        call self%compute_quadrature_metrics()
        call self%compute_quadrature_normals()
        call self%compute_quadrature_gradients()

        call self%update_face_ale(elem)

        !
        ! Compute BR2 matrix
        !
        val_e = self%basis_s%interpolator('Value')
        val_f = self%basis_s%interpolator('Value',self%iface)

        tmp = matmul(elem%invmass,transpose(val_f))
        self%br2_face = matmul(val_f,tmp)
        self%br2_vol  = matmul(val_e,tmp)



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

        integer(ik)                 :: inode, iface, nnodes, ierr
        character(:),   allocatable :: coordinate_system

        real(rk),   dimension(:),   allocatable ::  &
            d1dxi, d1deta, d1dzeta,                 &
            d2dxi, d2deta, d2dzeta,                 &
            d3dxi, d3deta, d3dzeta,                 &
            scaling_12, scaling_13, scaling_23, scaling_123, &
            invjac, weights

        real(rk),   dimension(:,:), allocatable :: val, ddxi, ddeta, ddzeta

        iface  = self%iface
        nnodes  = self%basis_c%nnodes_if()
        weights = self%basis_c%weights()
        val     = self%basis_c%interpolator('Value', self%iface)
        ddxi    = self%basis_c%interpolator('ddxi',  self%iface)
        ddeta   = self%basis_c%interpolator('ddeta', self%iface)
        ddzeta  = self%basis_c%interpolator('ddzeta',self%iface)

        !
        ! Evaluate directional derivatives of coordinates at quadrature nodes.
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
                call chidg_signal_one(FATAL,"face%compute_quadrature_metrics: Invalid coordinate system. Choose 'Cartesian' or 'Cylindrical'.",self%coordinate_system)
        end select

        !
        ! compute inverse cell mapping jacobian terms
        !
        invjac = scaling_123 * (d1dxi*d2deta*d3dzeta - d1deta*d2dxi*d3dzeta - &
                                d1dxi*d2dzeta*d3deta + d1dzeta*d2dxi*d3deta + &
                                d1deta*d2dzeta*d3dxi - d1dzeta*d2deta*d3dxi)
        self%jinv = invjac



        !
        ! At each quadrature node, compute metric terms.
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





        !
        ! Check for negative jacobians
        !
        if (any(self%jinv < ZERO)) call chidg_signal(FATAL,"face%compute_quadrature_metrics: Negative element jacobians detected. Check element quality and orientation.")



    end subroutine compute_quadrature_metrics
    !******************************************************************************************












    !> Compute normal vector components at face quadrature nodes
    !!
    !!  NOTE: be sure to differentiate between normals self%norm and unit-normals self%unorm
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti  
    !!  @date   11/5/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_quadrature_normals(self)
        class(face_t),  intent(inout)   :: self

        integer(ik)                 :: inode, iface, nnodes, ierr
        character(:),   allocatable :: coordinate_system

        real(rk),   dimension(:),   allocatable ::              &
            d1dxi, d1deta, d1dzeta,                             &
            d2dxi, d2deta, d2dzeta,                             &
            d3dxi, d3deta, d3dzeta,                             &
            scaling_12, scaling_13, scaling_23, scaling_123,    &
            norm_mag, weights

        real(rk),   dimension(:,:),   allocatable ::  &
            val, ddxi, ddeta, ddzeta


        iface   = self%iface
        nnodes  = self%basis_c%nnodes_if()
        weights = self%basis_c%weights(iface)
        val     = self%basis_c%interpolator('Value', iface)
        ddxi    = self%basis_c%interpolator('ddxi',  iface)
        ddeta   = self%basis_c%interpolator('ddeta', iface)
        ddzeta  = self%basis_c%interpolator('ddzeta',iface)


        
        !
        ! Evaluate directional derivatives of coordinates at quadrature nodes.
        !
        d1dxi   = matmul(ddxi,   self%coords%getvar(1,itime=1))
        d1deta  = matmul(ddeta,  self%coords%getvar(1,itime=1))
        d1dzeta = matmul(ddzeta, self%coords%getvar(1,itime=1))

        d2dxi   = matmul(ddxi,   self%coords%getvar(2,itime=1))
        d2deta  = matmul(ddeta,  self%coords%getvar(2,itime=1))
        d2dzeta = matmul(ddzeta, self%coords%getvar(2,itime=1))

        d3dxi   = matmul(ddxi,   self%coords%getvar(3,itime=1))
        d3deta  = matmul(ddeta,  self%coords%getvar(3,itime=1))
        d3dzeta = matmul(ddzeta, self%coords%getvar(3,itime=1))



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
                call chidg_signal(FATAL,"face%compute_quadrature_normals: Invalid coordinate system. Choose 'Cartesian' or 'Cylindrical'.")
        end select



        !
        ! Compute normal vectors for each face
        !
        select case (self%iface)
            case (XI_MIN, XI_MAX)

                do inode = 1,nnodes
                    self%norm(inode,XI_DIR)   = scaling_23(inode)*d2deta(inode)*d3dzeta(inode) - scaling_23(inode)*d2dzeta(inode)*d3deta(inode)
                    self%norm(inode,ETA_DIR)  = scaling_13(inode)*d1dzeta(inode)*d3deta(inode) - scaling_13(inode)*d1deta(inode)*d3dzeta(inode)
                    self%norm(inode,ZETA_DIR) = scaling_12(inode)*d1deta(inode)*d2dzeta(inode) - scaling_12(inode)*d1dzeta(inode)*d2deta(inode)
                end do

            case (ETA_MIN, ETA_MAX)

                do inode = 1,nnodes
                    self%norm(inode,XI_DIR)   = scaling_23(inode)*d2dzeta(inode)*d3dxi(inode)  - scaling_23(inode)*d2dxi(inode)*d3dzeta(inode)
                    self%norm(inode,ETA_DIR)  = scaling_13(inode)*d1dxi(inode)*d3dzeta(inode)  - scaling_13(inode)*d1dzeta(inode)*d3dxi(inode)
                    self%norm(inode,ZETA_DIR) = scaling_12(inode)*d1dzeta(inode)*d2dxi(inode)  - scaling_12(inode)*d1dxi(inode)*d2dzeta(inode)
                end do

            case (ZETA_MIN, ZETA_MAX)

                do inode = 1,nnodes
                    self%norm(inode,XI_DIR)   = scaling_23(inode)*d2dxi(inode)*d3deta(inode)   - scaling_23(inode)*d3dxi(inode)*d2deta(inode)
                    self%norm(inode,ETA_DIR)  = scaling_13(inode)*d3dxi(inode)*d1deta(inode)   - scaling_13(inode)*d1dxi(inode)*d3deta(inode)
                    self%norm(inode,ZETA_DIR) = scaling_12(inode)*d1dxi(inode)*d2deta(inode)   - scaling_12(inode)*d2dxi(inode)*d1deta(inode)
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


        !
        ! The 'norm' component is really a normal vector scaled by the FACE inverse jacobian.
        ! This is really a differential area scaling. We can compute the area 
        ! scaling(jinv for the face, different than jinv for the element),
        ! by taking the magnitude of the 'norm' vector.
        !
        self%differential_areas = norm_mag

        !
        ! Compute the total face area by integrating the differential areas over the face
        !
        self%total_area = sum(abs(self%differential_areas * weights))


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
    subroutine compute_quadrature_gradients(self)
        class(face_t),      intent(inout)   :: self

        integer(ik)                                 :: iterm,inode,iface,nnodes
        real(rk),   allocatable,    dimension(:,:)  :: ddxi, ddeta, ddzeta

        iface  = self%iface
        nnodes = self%basis_s%nnodes_if()
        ddxi   = self%basis_s%interpolator('ddxi',  iface)
        ddeta  = self%basis_s%interpolator('ddeta', iface)
        ddzeta = self%basis_s%interpolator('ddzeta',iface)



        do iterm = 1,self%nterms_s
            do inode = 1,nnodes

                self%grad1(inode,iterm) = self%metric(1,1,inode) * ddxi(inode,iterm)   + &
                                          self%metric(2,1,inode) * ddeta(inode,iterm)  + &
                                          self%metric(3,1,inode) * ddzeta(inode,iterm)

                self%grad2(inode,iterm) = self%metric(1,2,inode) * ddxi(inode,iterm)   + &
                                          self%metric(2,2,inode) * ddeta(inode,iterm)  + &
                                          self%metric(3,2,inode) * ddzeta(inode,iterm)

                self%grad3(inode,iterm) = self%metric(1,3,inode) * ddxi(inode,iterm)   + &
                                          self%metric(2,3,inode) * ddeta(inode,iterm)  + &
                                          self%metric(3,3,inode) * ddzeta(inode,iterm)

            end do
        end do

    end subroutine compute_quadrature_gradients
    !*******************************************************************************************







    !> Compute cartesian coordinates at face quadrature nodes
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti  
    !!  @date   11/5/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_quadrature_coords(self)
        class(face_t),  intent(inout)   :: self

        integer(ik) :: iface, inode
        real(rk),   allocatable, dimension(:)   :: c1, c2, c3
        real(rk),   allocatable, dimension(:,:) :: val

        iface = self%iface

        !
        ! compute real coordinates associated with quadrature points
        !
        val = self%basis_c%interpolator('Value',iface)
        c1 = matmul(val,self%coords%getvar(1,itime = 1))
        c2 = matmul(val,self%coords%getvar(2,itime = 1))
        c3 = matmul(val,self%coords%getvar(3,itime = 1))

        
        !
        ! For each quadrature node, store real coordinates
        !
        do inode = 1,self%basis_s%nnodes_if()
            self%quad_pts(inode,1:3) = [c1(inode), c2(inode), c3(inode)]
        end do !inode

    end subroutine compute_quadrature_coords
    !******************************************************************************************







    !>  Return neighbor element index
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
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
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
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

    !>  Initialize ALE data from nodal displacements.
    !!
    !!  @author Eric Wolf (AFRL)
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2017
    !!
    !--------------------------------------------------------------------------------------
    subroutine update_face_ale(self,elem)
        class(face_t),      intent(inout)   :: self
        type(element_t),    intent(in)  :: elem

        call self%update_face_ale_coords(elem)
        call self%update_face_ale_quadrature()

    end subroutine update_face_ale
    !**************************************************************************************


    !>
    !!
    !!  @author Eric Wolf (AFRL)
    !!  @date   7/5/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine update_face_ale_quadrature(self)
        class(face_t),       intent(inout)      :: self

        call self%compute_quadrature_coords_ale()
        call self%compute_quadrature_metrics_ale()

    end subroutine update_face_ale_quadrature
    !******************************************************************************************



    !>
    !!
    !!
    !!  @author Eric Wolf (AFRL)
    !!  @date   7/5/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_quadrature_coords_ale(self)
        class(face_t),   intent(inout)   :: self


        integer(ik)                             :: nnodes, inode
        real(rk),   allocatable, dimension(:)   :: x, y, z, vg1, vg2, vg3
        real(rk),   allocatable, dimension(:,:) :: val

        nnodes = self%basis_s%nnodes_if()
        val    = self%basis_c%interpolator('Value',self%iface)

        !
        ! compute cartesian coordinates associated with quadrature points
        !
        x = matmul(val,self%ale_coords%getvar(1,itime = 1))
        y = matmul(val,self%ale_coords%getvar(2,itime = 1))
        z = matmul(val,self%ale_coords%getvar(3,itime = 1))

        !
        ! Initialize each point with cartesian coordinates
        !
        do inode = 1,nnodes
            self%ale_quad_pts(inode,1:3) = [x(inode), y(inode), z(inode)]
        end do


        ! Grid velocity

        ! compute cartesian coordinates associated with quadrature points
        !
        vg1 = matmul(val,self%ale_vel_coords%getvar(1,itime = 1))
        vg2 = matmul(val,self%ale_vel_coords%getvar(2,itime = 1))
        vg3 = matmul(val,self%ale_vel_coords%getvar(3,itime = 1))


        !
        ! Initialize each point with cartesian coordinates
        !
        do inode = 1,nnodes
            self%grid_vel(inode,1) = vg1(inode)
            self%grid_vel(inode,2) = vg2(inode)
            self%grid_vel(inode,3) = vg3(inode)
        end do 


    end subroutine compute_quadrature_coords_ale
    !****************************************************************************************


    !> Compute metric terms and cell jacobians at face quadrature nodes
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  TODO: Generalize 2D physical coordinates. Currently assumes x-y.
    !!
    !----------------------------------------------------------------------------------------
    subroutine compute_quadrature_metrics_ale(self)
        class(face_t),  intent(inout)   :: self

        integer(ik)                 :: inode, iface, ierr
        integer(ik)                 :: nnodes

        real(rk),   dimension(:),   allocatable ::  &
            d1dxi, d1deta, d1dzeta,                 &
            d2dxi, d2deta, d2dzeta,                 &
            d3dxi, d3deta, d3dzeta,                 &
            d1dxi_ale, d1deta_ale, d1dzeta_ale,                 &
            d2dxi_ale, d2deta_ale, d2dzeta_ale,                 &
            d3dxi_ale, d3deta_ale, d3dzeta_ale,                 &
            invjac_ale,                                         &
            scaling_12, scaling_13, scaling_23, scaling_123,    &
            fvals, temp, weights

        real(rk),   dimension(:,:), allocatable :: val, ddxi, ddeta, ddzeta, tempmat
        real(rk), dimension(:,:,:), allocatable :: jacobian_matrix, jacobian_matrix_ale

        !
        ! Retrieve interpolators
        !
        iface   = self%iface
        nnodes  = self%basis_c%nnodes_if()
        weights = self%basis_c%weights()
        ddxi    = self%basis_c%interpolator('ddxi',iface)
        ddeta   = self%basis_c%interpolator('ddeta',iface)
        ddzeta  = self%basis_c%interpolator('ddzeta',iface)

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

        !self%det_jacobian_grid       = matmul(self%basis_s%interpolator('Value',iface),self%det_jacobian_grid_modes)
        self%det_jacobian_grid       = self%jinv_ale/self%jinv

        self%det_jacobian_grid_grad1 = matmul(self%grad1,self%det_jacobian_grid_modes)
        self%det_jacobian_grid_grad2 = matmul(self%grad2,self%det_jacobian_grid_modes)
        self%det_jacobian_grid_grad3 = matmul(self%grad3,self%det_jacobian_grid_modes)

    end subroutine compute_quadrature_metrics_ale
    !*****************************************************************************************************






    subroutine destructor(self)
        type(face_t), intent(inout) :: self


    end subroutine




end module type_face
