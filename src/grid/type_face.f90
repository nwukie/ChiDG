module type_face
#include  <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX,                 &
                                      ZETA_MIN, ZETA_MAX, XI_DIR, ETA_DIR, ZETA_DIR,    &
                                      NO_INTERIOR_NEIGHBOR, NO_PROC,                    &
                                      ZERO, ONE, TWO, ORPHAN, NO_MM_ASSIGNED, CARTESIAN, CYLINDRICAL, NO_ID
    use type_reference_element, only: reference_element_t
    use type_element,           only: element_t
    use type_densevector,       only: densevector_t
    use mod_inv,                only: inv, inv_3x3
    use mod_determinant,        only: det_3x3
    use ieee_arithmetic,        only: ieee_is_nan
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

        ! Self information
        integer(ik)             :: spacedim         ! Number of spatial dimensions
        integer(ik)             :: ftype            ! INTERIOR, BOUNDARY, CHIMERA, ORPHAN 
        integer(ik)             :: ChiID    = NO_ID ! Identifier for domain-local Chimera interfaces
        integer(ik)             :: bc_ID    = NO_ID ! Index for bc state group data%bc_state_group(bc_ID)
        integer(ik)             :: group_ID = NO_ID ! Index for bc patch group mesh%bc_patch_group(group_ID)
        integer(ik)             :: patch_ID = NO_ID ! Index for bc patch 
        integer(ik)             :: face_ID  = NO_ID ! Index for bc patch face
        integer(ik)             :: mm_ID    = NO_MM_ASSIGNED

        ! Owner-element information
        integer(ik)             :: face_location(5)! [idomain_g, idomain_l, iparent_g, iparent_l, iface]
        integer(ik)             :: idomain_g       ! Global index of the parent domain
        integer(ik)             :: idomain_l       ! Processor-local index of the parent domain
        integer(ik)             :: iparent_g       ! Domain-global index of the parent element
        integer(ik)             :: iparent_l       ! Processor-local index of the parent element
        integer(ik)             :: iface           ! XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, etc
        integer(ik)             :: nfields         ! Number of equations in equationset_t
        integer(ik)             :: nterms_s        ! Number of terms in solution polynomial expansion
        integer(ik)             :: dof_start       ! Starting DOF index in ChiDG-global index 
        integer(ik)             :: dof_local_start ! Starting DOF index in ChiDG-local index 
        integer(ik)             :: ntime
        integer(ik)             :: coordinate_system    ! CARTESIAN, CYLINDRICAL. parameters from mod_constants.f90


        ! Neighbor information
        integer(ik)             :: neighbor_location(7)      = 0         ! [idomain_g, idomain_l, ielement_g, ielement_l, iface, dof_start, dof_local_start]
        integer(ik)             :: ineighbor_proc            = NO_PROC   ! MPI processor rank of the neighboring element
        integer(ik)             :: ineighbor_domain_g        = 0         ! Global index of the neighboring element's domain
        integer(ik)             :: ineighbor_domain_l        = 0         ! Processor-local index of the neighboring element's domain
        integer(ik)             :: ineighbor_element_g       = 0         ! Domain-global index of the neighboring element
        integer(ik)             :: ineighbor_element_l       = 0         ! Processor-local index of the neighboring element
        integer(ik)             :: ineighbor_face            = 0
        integer(ik)             :: ineighbor_nfields         = 0
        integer(ik)             :: ineighbor_nterms_s        = 0
        integer(ik)             :: ineighbor_nnodes_r        = 0
        integer(ik)             :: ineighbor_ntime           = 0
        integer(ik)             :: ineighbor_dof_start       = NO_ID
        integer(ik)             :: ineighbor_dof_local_start = NO_ID
        integer(ik)             :: ineighbor_pelem_ID        = NO_ID
        integer(ik)             :: recv_comm                 = NO_ID
        integer(ik)             :: recv_domain               = NO_ID
        integer(ik)             :: recv_element              = NO_ID
        integer(ik)             :: recv_dof                  = NO_ID


        ! Neighbor information: if neighbor is off-processor
        real(rk)                :: neighbor_h(3)           ! Approximate size of neighbor bounding box
        real(rk),   allocatable :: neighbor_grad1(:,:)     ! Grad of basis functions in at quadrature nodes
        real(rk),   allocatable :: neighbor_grad2(:,:)     ! Grad of basis functions in at quadrature nodes
        real(rk),   allocatable :: neighbor_grad3(:,:)     ! Grad of basis functions in at quadrature nodes
        real(rk),   allocatable :: neighbor_br2_face(:,:)  ! Matrix for computing/obtaining br2 modes at face nodes
        real(rk),   allocatable :: neighbor_br2_vol(:,:)   ! Matrix for computing/obtaining br2 modes at volume nodes
        real(rk),   allocatable :: neighbor_invmass(:,:)    
        real(rk),   allocatable :: neighbor_coords(:,:)         ! Modal representation of neighbor's coordinates
        real(rk),   allocatable :: neighbor_dgrad1_dx(:,:,:,:)  ! Derivatives of grad of basis functions in at quadrature nodes wrt to neighbor nodes
        real(rk),   allocatable :: neighbor_dgrad2_dx(:,:,:,:)  ! Derivatives of grad of basis functions in at quadrature nodes wrt to neighbor nodes
        real(rk),   allocatable :: neighbor_dgrad3_dx(:,:,:,:)  ! Derivatives of grad of basis functions in at quadrature nodes wrt to neighbor nodes
        real(rk),   allocatable :: neighbor_dbr2_f_dx(:,:,:,:)  ! Derivatives of br2_face matrix wrt to neigbor nodes
        real(rk),   allocatable :: neighbor_dnorm_dx(:,:,:,:)   ! Derivatives of neighbor face wrt neighbor support nodes


        ! Neighbor ALE: if neighbor is off-processor
        real(rk),   allocatable :: neighbor_interp_coords_vel(:,:)   
        real(rk),   allocatable :: neighbor_ale_Dinv(:,:,:)
        real(rk),   allocatable :: neighbor_ale_g(:)
        real(rk),   allocatable :: neighbor_ale_g_grad1(:)
        real(rk),   allocatable :: neighbor_ale_g_grad2(:)
        real(rk),   allocatable :: neighbor_ale_g_grad3(:)


        ! Chimera face offset. For periodic boundary condition.
        logical                 :: periodic_offset  = .false.
        real(rk)                :: chimera_offset_1 = 0._rk
        real(rk)                :: chimera_offset_2 = 0._rk
        real(rk)                :: chimera_offset_3 = 0._rk


        ! Modal representations of element coordinates/velocity
        type(densevector_t)     :: coords                   ! Modal expansion of coordinates 
        type(densevector_t)     :: ale_coords               ! Modal representation of cartesian coordinates (nterms_var,(x,y,z))
        type(densevector_t)     :: ale_vel_coords           ! Modal representation of cartesian coordinates (nterms_var,(x,y,z))


        ! Element data at interpolation nodes
        real(rk),   allocatable :: interp_coords(:,:)       ! Undeformed coordinates at face interpolation nodes
        real(rk),   allocatable :: interp_coords_def(:,:)   ! Deformed coordinates at face interpolation nodes
        real(rk),   allocatable :: interp_coords_vel(:,:)   ! Coordinate velocities at face interpolation nodes
        real(rk),   allocatable :: jinv(:)                  ! Volume scaling: Undeformed/Reference
        real(rk),   allocatable :: jinv_def(:)              ! Volume scaling: Deformed/Reference
        real(rk),   allocatable :: metric(:,:,:)            ! Face metric terms  : undeformed face
        real(rk),   allocatable :: norm(:,:)                ! Face normal vector : scaled by differential area : undeformed face
        real(rk),   allocatable :: unorm(:,:)               ! Face normal vector : unit length : undeformed face
        real(rk),   allocatable :: unorm_def(:,:)           ! Face normal vector : unit length : deformed face


        ! Matrices of cartesian gradients of basis/test functions
        real(rk),   allocatable :: grad1(:,:)           ! Deriv of basis functions in at interpolation nodes
        real(rk),   allocatable :: grad2(:,:)           ! Deriv of basis functions in at interpolation nodes
        real(rk),   allocatable :: grad3(:,:)           ! Deriv of basis functions in at interpolation nodes


        ! BR2 matrix
        real(rk),   allocatable :: br2_face(:,:)
        real(rk),   allocatable :: br2_vol(:,:)


        ! Face area
        real(rk)                :: total_area
        real(rk)                :: centroid(3)
        real(rk),   allocatable :: differential_areas(:)
        real(rk),   allocatable :: ale_area_ratio(:)

        ! Smoothed h-field
        real(rk),   allocatable :: h_smooth(:,:)        ! (ngq, 3)
        real(rk),   allocatable :: size_smooth(:)       ! (ngq)


        ! Arbitrary Lagrangian Eulerian data
        !   : This defines a mapping from some deformed element back to the original
        !   : undeformed element with the idea that the governing equations are transformed
        !   : and solved on the undeformed element.
        real(rk),   allocatable :: ale_Dinv(:,:,:)
        real(rk),   allocatable :: ale_g(:)
        real(rk),   allocatable :: ale_g_grad1(:)
        real(rk),   allocatable :: ale_g_grad2(:)
        real(rk),   allocatable :: ale_g_grad3(:)
        real(rk),   allocatable :: ale_g_modes(:)


        ! Grid geometry sensitivities, adjoint-based
        !   : Computes the derivatives of the metrics and jinv wrt reference grid nodes
        real(rk),       allocatable :: dmetric_dx(:,:,:,:,:)    ! Derivatives of inverted jacobian matrix for each quadrature node wrt to each element node (mat_i,mat_j,quad_pt,diff_node,ncoords)
        real(rk),       allocatable :: djinv_dx(:,:,:)          ! Derivative of differential volume ratio wrt element's node coordinates. (quad_pt,diff_nodes,ncoords) 
        real(rk),       allocatable :: dnorm_dx(:,:,:,:)        ! Derivatives of face normal vector : scaled by differential area : undeformed face (quad_pt,dir,diff_node,ncoords)
        real(rk),       allocatable :: dgrad1_dx(:,:,:,:)       ! Derivative of grad of basis functions in at quadrature nodes wrt grid nodes 
        real(rk),       allocatable :: dgrad2_dx(:,:,:,:)       ! Derivative of grad of basis functions in at quadrature nodes wrt grid nodes
        real(rk),       allocatable :: dgrad3_dx(:,:,:,:)       ! Derivative of grad of basis functions in at quadrature nodes wrt grid nodes
        real(rk),       allocatable :: dbr2_v_dx(:,:,:,:)       ! Derivative of br2_vol matrix
        real(rk),       allocatable :: dbr2_f_dx(:,:,:,:)       ! Derivative of br2_face matrix
        


        ! Solution/Coordinate basis objects
        type(reference_element_t), pointer  :: basis_s => null()
        type(reference_element_t), pointer  :: basis_c => null()


        ! Logical tests
        logical :: geom_initialized    = .false.
        logical :: sol_initialized     = .false.
        logical :: neighborInitialized = .false.


    contains

        ! Initialization procedures
        procedure, public   :: init_geom
        procedure, public   :: init_sol

        ! Undeformed element procedures
        procedure, public   :: update_interpolations
        procedure, private  :: interpolate_coords        
        procedure, private  :: interpolate_metrics       
        procedure, private  :: interpolate_normals       
        procedure, private  :: interpolate_gradients     

        procedure           :: compute_projected_areas

        
        ! Deformed element/ALE procedures
        procedure, public   :: set_displacements_velocities
        procedure, public   :: update_interpolations_ale
        procedure, private  :: interpolate_normals_ale  
        procedure, private  :: interpolate_coords_ale
        procedure, private  :: interpolate_metrics_ale

        ! Neighbor data procedures
        procedure           :: set_neighbor             ! Set neighbor location data
        procedure           :: get_neighbor_element_g   ! Return neighbor element index
        procedure           :: get_neighbor_element_l   ! Return neighbor element index
        procedure           :: get_neighbor_face        ! Return neighbor face index

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
        self%ineighbor_domain_g        = NO_INTERIOR_NEIGHBOR
        self%ineighbor_domain_l        = NO_INTERIOR_NEIGHBOR
        self%ineighbor_element_g       = NO_INTERIOR_NEIGHBOR
        self%ineighbor_element_l       = NO_INTERIOR_NEIGHBOR
        self%ineighbor_face            = NO_INTERIOR_NEIGHBOR
        self%ineighbor_proc            = NO_PROC
        self%ineighbor_dof_start       = NO_ID
        self%ineighbor_dof_local_start = NO_ID
        self%neighbor_location = [self%ineighbor_domain_g, self%ineighbor_domain_l,   &
                                  self%ineighbor_element_g, self%ineighbor_element_l, &
                                  self%ineighbor_face, self%ineighbor_dof_start,      &
                                  self%ineighbor_dof_local_start]
        

        !
        ! Set modal representation of element coordinates, displacements, velocities:
        !   1: set reference coordinates
        !   2: set default ALE (displacements, velocities)
        !
        self%coords = elem%coords
        call self%set_displacements_velocities(elem)


        !
        ! Set coordinate system, confirm initialization.
        !
        self%coordinate_system = elem%coordinate_system
        self%geom_initialized  = .true.

    end subroutine init_geom
    !******************************************************************************************








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
        self%nfields    = elem%nfields
        self%nterms_s   = elem%nterms_s
        self%dof_start  = elem%dof_start
        self%ntime      = elem%ntime
        self%basis_s    => elem%basis_s
        self%basis_c    => elem%basis_c



        !
        ! (Re)Allocate storage for face data structures.
        !
        if (allocated(self%jinv))                   &
            deallocate(self%jinv,                   &
                       self%jinv_def,               &
                       self%interp_coords,          &
                       self%metric,                 &
                       self%norm,                   &
                       self%unorm,                  &
                       self%unorm_def,              &
                       self%interp_coords_def,      &
                       self%ale_Dinv,               &
                       self%ale_g,                  &
                       self%ale_g_grad1,            &
                       self%ale_g_grad2,            &
                       self%ale_g_grad3,            &
                       self%ale_g_modes,            &
                       self%interp_coords_vel,      &
                       self%neighbor_ale_Dinv,      &
                       self%neighbor_ale_g,         &
                       self%neighbor_ale_g_grad1,   &
                       self%neighbor_ale_g_grad2,   &
                       self%neighbor_ale_g_grad3,   &
                       self%neighbor_interp_coords_vel,  &
                       self%grad1,                  &
                       self%grad2,                  &
                       self%grad3,                  &
                       self%h_smooth,               &
                       self%size_smooth             &
                       ) 



        nnodes = self%basis_s%nnodes_face()
        allocate(self%jinv(nnodes),                     &
                 self%jinv_def(nnodes),                 &
                 self%interp_coords(nnodes,3),          &
                 self%metric(3,3,nnodes),               &
                 self%norm(nnodes,3),                   &
                 self%unorm(nnodes,3),                  &
                 self%unorm_def(nnodes,3),              &
                 self%interp_coords_def(nnodes,3),      &
                 self%ale_Dinv(3,3,nnodes),             &
                 self%ale_g(nnodes),                    &
                 self%ale_g_grad1(nnodes),              &
                 self%ale_g_grad2(nnodes),              &
                 self%ale_g_grad3(nnodes),              &
                 self%ale_g_modes(self%nterms_s),       &
                 self%interp_coords_vel(nnodes,3),      &
                 self%neighbor_ale_Dinv(3,3,nnodes),    &
                 self%neighbor_ale_g(nnodes),           &
                 self%neighbor_ale_g_grad1(nnodes),     &
                 self%neighbor_ale_g_grad2(nnodes),     &
                 self%neighbor_ale_g_grad3(nnodes),     &
                 self%neighbor_interp_coords_vel(nnodes,3),  &
                 self%grad1(nnodes,self%nterms_s),      &
                 self%grad2(nnodes,self%nterms_s),      &
                 self%grad3(nnodes,self%nterms_s),      &
                 self%h_smooth(nnodes,3),               &
                 self%size_smooth(nnodes),              &
                 stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Compute metrics, normals, node coordinates
        !
        call self%update_interpolations()
        call self%update_interpolations_ale(elem)

        !
        ! Compute BR2 matrix
        !   val * invmass * val_trans
        !
        val_e = self%basis_s%interpolator_element('Value')
        val_f = self%basis_s%interpolator_face('Value',self%iface)

        tmp = matmul(elem%invmass,transpose(val_f))
        self%br2_face = matmul(val_f,tmp)
        self%br2_vol  = matmul(val_e,tmp)



        !
        ! Confirm face numerics were initialized
        !
        self%sol_initialized  = .true.

    end subroutine init_sol
    !******************************************************************************************






    !>  Update interpolations of data related to the undeformed element/face to the 
    !!  interpolation node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/15/2017
    !!
    !-----------------------------------------------------------------------------------------
    subroutine update_interpolations(self)
        class(face_t),  intent(inout)   :: self

        call self%interpolate_coords()
        call self%interpolate_metrics()
        call self%interpolate_normals()
        call self%interpolate_gradients()

    end subroutine update_interpolations
    !*****************************************************************************************






    !> Compute metric terms and cell jacobians at face quadrature nodes
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  TODO: Generalize 2D physical coordinates. Currently assumes x-y.
    !!
    !-----------------------------------------------------------------------------------------
    subroutine interpolate_metrics(self)
        class(face_t),  intent(inout)   :: self

        integer(ik)                 :: inode, nnodes, ierr
        character(:),   allocatable :: coordinate_system, user_msg

        real(rk),   dimension(:),       allocatable :: scaling_row2
        real(rk),   dimension(:,:),     allocatable :: val, ddxi, ddeta, ddzeta
        real(rk),   dimension(:,:,:),   allocatable :: jacobian


        nnodes  = self%basis_c%nnodes_face()
        val     = self%basis_c%interpolator_face('Value', self%iface)
        ddxi    = self%basis_c%interpolator_face('ddxi',  self%iface)
        ddeta   = self%basis_c%interpolator_face('ddeta', self%iface)
        ddzeta  = self%basis_c%interpolator_face('ddzeta',self%iface)


        !
        ! Compute element jacobian matrix at interpolation nodes
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
                user_msg = "face%interpolate_metrics: Invalid coordinate system."
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
        user_msg = "face%interpolate_metrics: Negative element &
                    volume detected. Check element quality and orientation."
        if (any(self%jinv < ZERO)) call chidg_signal_three(FATAL,user_msg,self%idomain_g,self%iparent_g,self%iface)



        !
        ! Invert jacobian matrix at each interpolation node
        !
        do inode = 1,nnodes
            self%metric(:,:,inode) = inv_3x3(jacobian(:,:,inode))
        end do



    end subroutine interpolate_metrics
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
    subroutine interpolate_normals(self)
        class(face_t),  intent(inout)   :: self

        integer(ik)                                 :: inode, nnodes, ierr
        character(:),   allocatable                 :: coordinate_system, user_msg
        real(rk),       allocatable, dimension(:)   :: norm_mag, weights


        nnodes  = self%basis_c%nnodes_face()
        weights = self%basis_c%weights_face(self%iface)


        !
        ! Compute normal vectors for each face
        !
        select case (self%iface)
            case (XI_MIN, XI_MAX)
                !self%norm(:,XI_DIR)   = self%jinv(:)*self%metric(1,1,:)
                !self%norm(:,ETA_DIR)  = self%jinv(:)*self%metric(1,2,:)
                !self%norm(:,ZETA_DIR) = self%jinv(:)*self%metric(1,3,:)
                do inode = 1,size(self%jinv)
                    self%norm(inode,:) = matmul(transpose(self%metric(:,:,inode)), [self%jinv(inode), ZERO, ZERO])
                end do

            case (ETA_MIN, ETA_MAX)
                !self%norm(:,XI_DIR)   = self%jinv(:)*self%metric(2,1,:)
                !self%norm(:,ETA_DIR)  = self%jinv(:)*self%metric(2,2,:)
                !self%norm(:,ZETA_DIR) = self%jinv(:)*self%metric(2,3,:)
                do inode = 1,size(self%jinv)
                    self%norm(inode,:) = matmul(transpose(self%metric(:,:,inode)), [ZERO, self%jinv(inode), ZERO])
                end do

            case (ZETA_MIN, ZETA_MAX)
                !self%norm(:,XI_DIR)   = self%jinv(:)*self%metric(3,1,:)
                !self%norm(:,ETA_DIR)  = self%jinv(:)*self%metric(3,2,:)
                !self%norm(:,ZETA_DIR) = self%jinv(:)*self%metric(3,3,:)
                do inode = 1,size(self%jinv)
                    self%norm(inode,:) = matmul(transpose(self%metric(:,:,inode)), [ZERO, ZERO, self%jinv(inode)])
                end do

            case default
                user_msg = "face%interpolate_normals: Invalid face index in face initialization."
                call chidg_signal(FATAL,user_msg)
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
        norm_mag = self%norm(:,1) ! Allocate to avoid DEBUG error
        norm_mag = sqrt(self%norm(:,XI_DIR)**TWO + self%norm(:,ETA_DIR)**TWO + self%norm(:,ZETA_DIR)**TWO)
        self%unorm(:,XI_DIR)   = self%norm(:,XI_DIR  )/norm_mag
        self%unorm(:,ETA_DIR)  = self%norm(:,ETA_DIR )/norm_mag
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


    end subroutine interpolate_normals
    !******************************************************************************************






    !>  Initialize ALE data from nodal displacements.
    !!
    !!  @author Eric Wolf (AFRL)
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2017
    !!
    !--------------------------------------------------------------------------------------
    subroutine set_displacements_velocities(self,elem)
        class(face_t),      intent(inout)   :: self
        type(element_t),    intent(in)  :: elem

        self%ale_coords     = elem%coords_def
        self%ale_vel_coords = elem%coords_vel

    end subroutine set_displacements_velocities
    !**************************************************************************************






    !> Compute normal vector components at face quadrature nodes: ALE
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
    subroutine interpolate_normals_ale(self)
        class(face_t),  intent(inout)   :: self

        integer(ik)                 :: inode, nnodes, ierr
        character(:),   allocatable :: coordinate_system, user_msg

        real(rk)                                :: metric_ale(3,3)
        real(rk),   dimension(:),   allocatable ::  norm_mag 
        real(rk),   dimension(:,:), allocatable ::  norm


        nnodes  = self%basis_c%nnodes_face()

        !
        ! Compute normal vectors for each face
        !
        allocate(norm(nnodes,3), stat=ierr)
        if (ierr /= 0) call AllocationError

        do inode = 1,nnodes

            ! Compute metric_ale: 
            !   dxi/dx = [Dinv]*[dxi/dX]
            !   dxi/dx = [dX/dx]*[dxi/dX]
            !metric_ale = matmul(self%ale_Dinv(:,:,inode),self%metric(:,:,inode))
            metric_ale = matmul(self%metric(:,:,inode),self%ale_Dinv(:,:,inode))

            select case (self%iface)
                case (XI_MIN, XI_MAX)
                    norm(inode,XI_DIR)   = self%jinv_def(inode)*metric_ale(1,1)
                    norm(inode,ETA_DIR)  = self%jinv_def(inode)*metric_ale(1,2)
                    norm(inode,ZETA_DIR) = self%jinv_def(inode)*metric_ale(1,3)

                case (ETA_MIN, ETA_MAX)
                    norm(inode,XI_DIR)   = self%jinv_def(inode)*metric_ale(2,1)
                    norm(inode,ETA_DIR)  = self%jinv_def(inode)*metric_ale(2,2)
                    norm(inode,ZETA_DIR) = self%jinv_def(inode)*metric_ale(2,3)

                case (ZETA_MIN, ZETA_MAX)
                    norm(inode,XI_DIR)   = self%jinv_def(inode)*metric_ale(3,1)
                    norm(inode,ETA_DIR)  = self%jinv_def(inode)*metric_ale(3,2)
                    norm(inode,ZETA_DIR) = self%jinv_def(inode)*metric_ale(3,3)

                case default
                    user_msg = "face%interpolate_normals_ale: Invalid face index in face initialization."
                    call chidg_signal(FATAL,user_msg)
            end select

        end do

        
        !
        ! Reverse normal vectors for faces XI_MIN,ETA_MIN,ZETA_MIN
        !
        if (self%iface == XI_MIN   .or. &
            self%iface == ETA_MIN  .or. &
            self%iface == ZETA_MIN) then
            norm(:,XI_DIR)   = -norm(:,XI_DIR)
            norm(:,ETA_DIR)  = -norm(:,ETA_DIR)
            norm(:,ZETA_DIR) = -norm(:,ZETA_DIR)
        end if



        !
        ! Compute vector magnitude, which is the differential area
        !
        norm_mag = self%norm(:,1) ! Allocate to avoid DEBUG error
        norm_mag = sqrt(norm(:,XI_DIR)**TWO + norm(:,ETA_DIR)**TWO + norm(:,ZETA_DIR)**TWO)
        self%unorm_def(:,XI_DIR)   = norm(:,XI_DIR  )/norm_mag
        self%unorm_def(:,ETA_DIR)  = norm(:,ETA_DIR )/norm_mag
        self%unorm_def(:,ZETA_DIR) = norm(:,ZETA_DIR)/norm_mag



        !
        ! Compute da/dA
        !
        self%ale_area_ratio = norm_mag/self%differential_areas


    end subroutine interpolate_normals_ale
    !******************************************************************************************








    !>  Compute matrices containing cartesian gradients of basis/test function
    !!  at each quadrature node.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine interpolate_gradients(self)
        class(face_t),      intent(inout)   :: self

        integer(ik)                                 :: iterm,inode,iface,nnodes
        real(rk),   allocatable,    dimension(:,:)  :: ddxi, ddeta, ddzeta

        iface  = self%iface
        nnodes = self%basis_s%nnodes_face()
        ddxi   = self%basis_s%interpolator_face('ddxi',  iface)
        ddeta  = self%basis_s%interpolator_face('ddeta', iface)
        ddzeta = self%basis_s%interpolator_face('ddzeta',iface)



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

    end subroutine interpolate_gradients
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
    subroutine interpolate_coords(self)
        class(face_t),  intent(inout)   :: self

        integer(ik) :: iface, inode
        real(rk),   allocatable, dimension(:)   :: c1, c2, c3
        real(rk),   allocatable, dimension(:,:) :: val

        iface = self%iface

        !
        ! compute real coordinates associated with quadrature points
        !
        val = self%basis_c%interpolator_face('Value',iface)
        c1 = matmul(val,self%coords%getvar(1,itime = 1))
        c2 = matmul(val,self%coords%getvar(2,itime = 1))
        c3 = matmul(val,self%coords%getvar(3,itime = 1))

        
        !
        ! For each quadrature node, store real coordinates
        !
        do inode = 1,self%basis_s%nnodes_face()
            self%interp_coords(inode,1:3) = [c1(inode), c2(inode), c3(inode)]
        end do !inode


        !
        ! Update face centroid, here we just take as an arithmetic average.
        !
        self%centroid(1) = sum(self%interp_coords(:,1))/size(self%interp_coords(:,1))
        self%centroid(2) = sum(self%interp_coords(:,2))/size(self%interp_coords(:,2))
        self%centroid(3) = sum(self%interp_coords(:,3))/size(self%interp_coords(:,3))

    end subroutine interpolate_coords
    !******************************************************************************************









    !>  Initialize ALE data from nodal displacements.
    !!
    !!  @author Eric Wolf (AFRL)
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2017
    !!
    !--------------------------------------------------------------------------------------
    subroutine update_interpolations_ale(self,elem)
        class(face_t),      intent(inout)   :: self
        type(element_t),    intent(in)      :: elem

        self%ale_g_modes = elem%ale_g_modes
        call self%interpolate_coords_ale()
        call self%interpolate_metrics_ale()
        call self%interpolate_normals_ale()


    end subroutine update_interpolations_ale
    !**************************************************************************************



    !>
    !!
    !!
    !!  @author Eric Wolf (AFRL)
    !!  @date   7/5/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine interpolate_coords_ale(self)
        class(face_t),   intent(inout)   :: self


        integer(ik)                             :: nnodes, inode
        real(rk),   allocatable, dimension(:)   :: x, y, z, vg1, vg2, vg3
        real(rk),   allocatable, dimension(:,:) :: val

        nnodes = self%basis_s%nnodes_face()
        val    = self%basis_c%interpolator_face('Value',self%iface)

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
            self%interp_coords_def(inode,1:3) = [x(inode), y(inode), z(inode)]
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
            self%interp_coords_vel(inode,1) = vg1(inode)
            self%interp_coords_vel(inode,2) = vg2(inode)
            self%interp_coords_vel(inode,3) = vg3(inode)
        end do 


    end subroutine interpolate_coords_ale
    !****************************************************************************************


    !> Compute metric terms and cell jacobians at face quadrature nodes
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  TODO: Generalize 2D physical coordinates. Currently assumes x-y.
    !!
    !----------------------------------------------------------------------------------------
    subroutine interpolate_metrics_ale(self)
        class(face_t),  intent(inout)   :: self

        integer(ik)                 :: inode, nnodes, ierr
        character(:),   allocatable :: coordinate_system, user_msg

        real(rk),   dimension(:),   allocatable ::                  &
!            dd1_dxidxi,   dd1_detadeta,   dd1_dzetadzeta,           &
!            dd2_dxidxi,   dd2_detadeta,   dd2_dzetadzeta,           &
!            dd3_dxidxi,   dd3_detadeta,   dd3_dzetadzeta,           &
!            dd1_dxideta,  dd1_dxidzeta,   dd1_detadzeta,            &
!            dd2_dxideta,  dd2_dxidzeta,   dd2_detadzeta,            &
!            dd3_dxideta,  dd3_dxidzeta,   dd3_detadzeta,            &
!            jinv_grad1, jinv_grad2, jinv_grad3,   &
!            jinv_def_grad1,   jinv_def_grad2,   jinv_def_grad3,     &
            ale_g_ddxi, ale_g_ddeta, ale_g_ddzeta, scaling_row2

        real(rk),   dimension(:,:), allocatable ::  &
            val, ddxi,    ddeta,    ddzeta,         &
            dxidxi,  detadeta, dzetadzeta,          &
            dxideta, dxidzeta, detadzeta, D_matrix

        real(rk), dimension(:,:,:), allocatable :: jacobian_ale


        !
        ! Retrieve interpolators
        !
        nnodes     = self%basis_c%nnodes_face()

        val        = self%basis_c%interpolator_face('Value',     self%iface)
        ddxi       = self%basis_c%interpolator_face('ddxi',      self%iface)
        ddeta      = self%basis_c%interpolator_face('ddeta',     self%iface)
        ddzeta     = self%basis_c%interpolator_face('ddzeta',    self%iface)

        dxidxi     = self%basis_c%interpolator_face('dxidxi',    self%iface)
        detadeta   = self%basis_c%interpolator_face('detadeta',  self%iface)
        dzetadzeta = self%basis_c%interpolator_face('dzetadzeta',self%iface)

        dxideta    = self%basis_c%interpolator_face('dxideta',   self%iface)
        dxidzeta   = self%basis_c%interpolator_face('dxidzeta',  self%iface)
        detadzeta  = self%basis_c%interpolator_face('detadzeta', self%iface)


        !
        ! Compute coordinate jacobian matrix at interpolation nodes
        !
        allocate(jacobian_ale(3,3,nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError
        jacobian_ale(1,1,:) = matmul(ddxi,   self%ale_coords%getvar(1,itime = 1))
        jacobian_ale(1,2,:) = matmul(ddeta,  self%ale_coords%getvar(1,itime = 1))
        jacobian_ale(1,3,:) = matmul(ddzeta, self%ale_coords%getvar(1,itime = 1))

        jacobian_ale(2,1,:) = matmul(ddxi,   self%ale_coords%getvar(2,itime = 1))
        jacobian_ale(2,2,:) = matmul(ddeta,  self%ale_coords%getvar(2,itime = 1))
        jacobian_ale(2,3,:) = matmul(ddzeta, self%ale_coords%getvar(2,itime = 1))

        jacobian_ale(3,1,:) = matmul(ddxi,   self%ale_coords%getvar(3,itime = 1))
        jacobian_ale(3,2,:) = matmul(ddeta,  self%ale_coords%getvar(3,itime = 1))
        jacobian_ale(3,3,:) = matmul(ddzeta, self%ale_coords%getvar(3,itime = 1))



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
        user_msg = "face%interpolate_metrics_ale: Negative element jacobians. &
                    Check element quality and origntation."
        if (any(self%jinv_def < ZERO)) call chidg_signal(FATAL,user_msg)


        !
        ! Compute element deformation gradient: dX/dx
        !   dX/dx = [dxi/dx][dX/dxi]
        !
        do inode = 1,nnodes
            D_matrix = matmul(jacobian_ale(:,:,inode),self%metric(:,:,inode))
            self%ale_Dinv(:,:,inode) = inv(D_matrix)

            ! Invert jacobian_matrix_ale for use in computing grad(jinv_ale)
            jacobian_ale(:,:,inode) = inv(jacobian_ale(:,:,inode))
        end do !inode


        ! Note, interpolator is from solution basis, not coordinate basis
        ! since ale_g_modes is in solution basis expansion.
        val    = self%basis_s%interpolator_face('Value',  self%iface)
        ddxi   = self%basis_s%interpolator_face('ddxi',   self%iface)
        ddeta  = self%basis_s%interpolator_face('ddeta',  self%iface)
        ddzeta = self%basis_s%interpolator_face('ddzeta', self%iface)

        self%ale_g    = matmul(val,    self%ale_g_modes)
        ale_g_ddxi    = matmul(ddxi,   self%ale_g_modes)
        ale_g_ddeta   = matmul(ddeta,  self%ale_g_modes)
        ale_g_ddzeta  = matmul(ddzeta, self%ale_g_modes)


        self%ale_g_grad1 = self%metric(1,1,:) * ale_g_ddxi  + &
                           self%metric(2,1,:) * ale_g_ddeta + &
                           self%metric(3,1,:) * ale_g_ddzeta

        self%ale_g_grad2 = self%metric(1,2,:) * ale_g_ddxi  + &
                           self%metric(2,2,:) * ale_g_ddeta + &
                           self%metric(3,2,:) * ale_g_ddzeta

        self%ale_g_grad3 = self%metric(1,3,:) * ale_g_ddxi  + &
                           self%metric(2,3,:) * ale_g_ddeta + &
                           self%metric(3,3,:) * ale_g_ddzeta






!        ! Second/mixed derivatives
!        dd1_dxidxi     = matmul(dxidxi,     self%coords%getvar(1,itime = 1))
!        dd1_detadeta   = matmul(detadeta,   self%coords%getvar(1,itime = 1))
!        dd1_dzetadzeta = matmul(dzetadzeta, self%coords%getvar(1,itime = 1))
!        dd1_dxideta    = matmul(dxideta,    self%coords%getvar(1,itime = 1))
!        dd1_dxidzeta   = matmul(dxidzeta,   self%coords%getvar(1,itime = 1))
!        dd1_detadzeta  = matmul(detadzeta,  self%coords%getvar(1,itime = 1))
!
!        dd2_dxidxi     = matmul(dxidxi,     self%coords%getvar(2,itime = 1))
!        dd2_detadeta   = matmul(detadeta,   self%coords%getvar(2,itime = 1))
!        dd2_dzetadzeta = matmul(dzetadzeta, self%coords%getvar(2,itime = 1))
!        dd2_dxideta    = matmul(dxideta,    self%coords%getvar(2,itime = 1))
!        dd2_dxidzeta   = matmul(dxidzeta,   self%coords%getvar(2,itime = 1))
!        dd2_detadzeta  = matmul(detadzeta,  self%coords%getvar(2,itime = 1))
!
!        dd3_dxidxi     = matmul(dxidxi,     self%coords%getvar(3,itime = 1))
!        dd3_detadeta   = matmul(detadeta,   self%coords%getvar(3,itime = 1))
!        dd3_dzetadzeta = matmul(dzetadzeta, self%coords%getvar(3,itime = 1))
!        dd3_dxideta    = matmul(dxideta,    self%coords%getvar(3,itime = 1))
!        dd3_dxidzeta   = matmul(dxidzeta,   self%coords%getvar(3,itime = 1))
!        dd3_detadzeta  = matmul(detadzeta,  self%coords%getvar(3,itime = 1))
!
!        jinv_grad1 = dd1_dxidxi*self%metric(1,1,:)     +  dd1_dxideta*self%metric(2,1,:)    +  dd1_dxidzeta*self%metric(3,1,:)   +  &
!                     dd2_dxidxi*self%metric(1,2,:)     +  dd2_dxideta*self%metric(2,2,:)    +  dd2_dxidzeta*self%metric(3,2,:)   +  &
!                     dd3_dxidxi*self%metric(1,3,:)     +  dd3_dxideta*self%metric(2,3,:)    +  dd3_dxidzeta*self%metric(3,3,:)
!
!        jinv_grad2 = dd1_dxideta*self%metric(1,1,:)    +  dd1_detadeta*self%metric(2,1,:)   +  dd1_detadzeta*self%metric(3,1,:)  +  &
!                     dd2_dxideta*self%metric(1,2,:)    +  dd2_detadeta*self%metric(2,2,:)   +  dd2_detadzeta*self%metric(3,2,:)  +  &
!                     dd3_dxideta*self%metric(1,3,:)    +  dd3_detadeta*self%metric(2,3,:)   +  dd3_detadzeta*self%metric(3,3,:)
!
!        jinv_grad3 = dd1_dxidzeta*self%metric(1,1,:)   +  dd1_detadzeta*self%metric(2,1,:)  +  dd1_dzetadzeta*self%metric(3,1,:) +  &
!                     dd2_dxidzeta*self%metric(1,2,:)   +  dd2_detadzeta*self%metric(2,2,:)  +  dd2_dzetadzeta*self%metric(3,2,:) +  &
!                     dd3_dxidzeta*self%metric(1,3,:)   +  dd3_detadzeta*self%metric(2,3,:)  +  dd3_dzetadzeta*self%metric(3,3,:)
!
!
!        ! Second/mixed derivatives
!        dd1_dxidxi     = matmul(dxidxi,     self%ale_coords%getvar(1,itime = 1))
!        dd1_detadeta   = matmul(detadeta,   self%ale_coords%getvar(1,itime = 1))
!        dd1_dzetadzeta = matmul(dzetadzeta, self%ale_coords%getvar(1,itime = 1))
!        dd1_dxideta    = matmul(dxideta,    self%ale_coords%getvar(1,itime = 1))
!        dd1_dxidzeta   = matmul(dxidzeta,   self%ale_coords%getvar(1,itime = 1))
!        dd1_detadzeta  = matmul(detadzeta,  self%ale_coords%getvar(1,itime = 1))
!
!        dd2_dxidxi     = matmul(dxidxi,     self%ale_coords%getvar(2,itime = 1))
!        dd2_detadeta   = matmul(detadeta,   self%ale_coords%getvar(2,itime = 1))
!        dd2_dzetadzeta = matmul(dzetadzeta, self%ale_coords%getvar(2,itime = 1))
!        dd2_dxideta    = matmul(dxideta,    self%ale_coords%getvar(2,itime = 1))
!        dd2_dxidzeta   = matmul(dxidzeta,   self%ale_coords%getvar(2,itime = 1))
!        dd2_detadzeta  = matmul(detadzeta,  self%ale_coords%getvar(2,itime = 1))
!
!        dd3_dxidxi     = matmul(dxidxi,     self%ale_coords%getvar(3,itime = 1))
!        dd3_detadeta   = matmul(detadeta,   self%ale_coords%getvar(3,itime = 1))
!        dd3_dzetadzeta = matmul(dzetadzeta, self%ale_coords%getvar(3,itime = 1))
!        dd3_dxideta    = matmul(dxideta,    self%ale_coords%getvar(3,itime = 1))
!        dd3_dxidzeta   = matmul(dxidzeta,   self%ale_coords%getvar(3,itime = 1))
!        dd3_detadzeta  = matmul(detadzeta,  self%ale_coords%getvar(3,itime = 1))
!
!
!
!        jinv_def_grad1 = dd1_dxidxi*jacobian_ale(1,1,:)    +  dd1_dxideta*jacobian_ale(2,1,:)    +  dd1_dxidzeta*jacobian_ale(3,1,:)   +  &
!                         dd2_dxidxi*jacobian_ale(1,2,:)    +  dd2_dxideta*jacobian_ale(2,2,:)    +  dd2_dxidzeta*jacobian_ale(3,2,:)   +  &
!                         dd3_dxidxi*jacobian_ale(1,3,:)    +  dd3_dxideta*jacobian_ale(2,3,:)    +  dd3_dxidzeta*jacobian_ale(3,3,:)
!
!        jinv_def_grad2 = dd1_dxideta*jacobian_ale(1,1,:)   +  dd1_detadeta*jacobian_ale(2,1,:)   +  dd1_detadzeta*jacobian_ale(3,1,:)  +  &
!                         dd2_dxideta*jacobian_ale(1,2,:)   +  dd2_detadeta*jacobian_ale(2,2,:)   +  dd2_detadzeta*jacobian_ale(3,2,:)  +  &
!                         dd3_dxideta*jacobian_ale(1,3,:)   +  dd3_detadeta*jacobian_ale(2,3,:)   +  dd3_detadzeta*jacobian_ale(3,3,:)
!
!        jinv_def_grad3 = dd1_dxidzeta*jacobian_ale(1,1,:)  +  dd1_detadzeta*jacobian_ale(2,1,:)  +  dd1_dzetadzeta*jacobian_ale(3,1,:) +  &
!                         dd2_dxidzeta*jacobian_ale(1,2,:)  +  dd2_detadzeta*jacobian_ale(2,2,:)  +  dd2_dzetadzeta*jacobian_ale(3,2,:) +  &
!                         dd3_dxidzeta*jacobian_ale(1,3,:)  +  dd3_detadzeta*jacobian_ale(2,3,:)  +  dd3_dzetadzeta*jacobian_ale(3,3,:)
!
!
!        !
!        ! Apply Quotiend Rule for computing gradient of det_jacobian_grid
!        !
!        !   det_jacobian_grid = jinv_ale/jinv
!        !
!        !   grad(det_jacobian_grid) = [grad(jinv_ale)*jinv - jinv_ale*grad(jinv)] / [jinv*jinv]
!        !
!        self%ale_g   = self%jinv_def/self%jinv
!        ale_g_ddxi   = (jinv_def_grad1*self%jinv  -  self%jinv_def*jinv_grad1)/(self%jinv**TWO)
!        ale_g_ddeta  = (jinv_def_grad2*self%jinv  -  self%jinv_def*jinv_grad2)/(self%jinv**TWO)
!        ale_g_ddzeta = (jinv_def_grad3*self%jinv  -  self%jinv_def*jinv_grad3)/(self%jinv**TWO)
!
!
!        ! Transform into gradient in physical space(undeformed geometry)
!        do inode = 1,size(ale_g_ddxi)
!            self%ale_g_grad1(inode) = self%metric(1,1,inode) * ale_g_ddxi(inode)  + &
!                                      self%metric(2,1,inode) * ale_g_ddeta(inode) + &
!                                      self%metric(3,1,inode) * ale_g_ddzeta(inode)
!
!            self%ale_g_grad2(inode) = self%metric(1,2,inode) * ale_g_ddxi(inode)  + &
!                                      self%metric(2,2,inode) * ale_g_ddeta(inode) + &
!                                      self%metric(3,2,inode) * ale_g_ddzeta(inode)
!
!            self%ale_g_grad3(inode) = self%metric(1,3,inode) * ale_g_ddxi(inode)  + &
!                                      self%metric(2,3,inode) * ale_g_ddeta(inode) + &
!                                      self%metric(3,3,inode) * ale_g_ddzeta(inode)
!        end do


    end subroutine interpolate_metrics_ale
    !*****************************************************************************************************





    !>
    !!
    !! @author  Eric M. Wolf
    !! @date    03/05/2019 
    !!
    !--------------------------------------------------------------------------------
    function compute_projected_areas(self) result(integral)
        class(face_t), intent(in) :: self

        real(rk)    :: integral(3)
        integer(ik) :: idir

        do idir = 1,3
            integral(idir) = sum(abs(self%norm(:,idir))*self%basis_s%weights_face(self%iface))
        end do

    end function compute_projected_areas
    !********************************************************************************



















    !>  Initialize neighbor location.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine set_neighbor(self,ftype,ineighbor_domain_g,ineighbor_domain_l,   &
                                       ineighbor_element_g,ineighbor_element_l, &
                                       ineighbor_face,ineighbor_nfields,        &
                                       ineighbor_ntime,ineighbor_nterms_s,      &
                                       ineighbor_nnodes_r,ineighbor_proc,       &
                                       ineighbor_dof_start, ineighbor_dof_local_start)
        class(face_t),  intent(inout)   :: self
        integer(ik),    intent(in)      :: ftype
        integer(ik),    intent(in)      :: ineighbor_domain_g
        integer(ik),    intent(in)      :: ineighbor_domain_l
        integer(ik),    intent(in)      :: ineighbor_element_g
        integer(ik),    intent(in)      :: ineighbor_element_l
        integer(ik),    intent(in)      :: ineighbor_face
        integer(ik),    intent(in)      :: ineighbor_nfields
        integer(ik),    intent(in)      :: ineighbor_ntime
        integer(ik),    intent(in)      :: ineighbor_nterms_s
        integer(ik),    intent(in)      :: ineighbor_nnodes_r
        integer(ik),    intent(in)      :: ineighbor_proc
        integer(ik),    intent(in)      :: ineighbor_dof_start
        integer(ik),    intent(in)      :: ineighbor_dof_local_start


        self%ftype                     = ftype
        self%ineighbor_domain_g        = ineighbor_domain_g
        self%ineighbor_domain_l        = ineighbor_domain_l
        self%ineighbor_element_g       = ineighbor_element_g
        self%ineighbor_element_l       = ineighbor_element_l
        self%ineighbor_face            = ineighbor_face
        self%ineighbor_nfields         = ineighbor_nfields
        self%ineighbor_ntime           = ineighbor_ntime
        self%ineighbor_nterms_s        = ineighbor_nterms_s
        self%ineighbor_nnodes_r        = ineighbor_nnodes_r
        self%ineighbor_proc            = ineighbor_proc
        self%ineighbor_dof_start       = ineighbor_dof_start
        self%ineighbor_dof_local_start = ineighbor_dof_local_start

        self%neighbor_location = [ineighbor_domain_g,  ineighbor_domain_l,  &
                                  ineighbor_element_g, ineighbor_element_l, &
                                  ineighbor_face, ineighbor_dof_start, ineighbor_dof_local_start]

        self%neighborInitialized = .true.

    end subroutine set_neighbor
    !*******************************************************************************************






    !>  Return neighbor element index
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
    !------------------------------------------------------------------------------------------
    function get_neighbor_face(self) result(neighbor_f)
        class(face_t),  intent(in)   ::  self

        integer(ik) :: neighbor_e
        integer(ik) :: neighbor_f

        neighbor_e = self%get_neighbor_element_l()

        if ( neighbor_e == NO_INTERIOR_NEIGHBOR ) then
            neighbor_f = NO_INTERIOR_NEIGHBOR
        else
            neighbor_f = self%ineighbor_face
        end if

    end function get_neighbor_face
    !******************************************************************************************






    subroutine destructor(self)
        type(face_t), intent(inout) :: self


    end subroutine




end module type_face
