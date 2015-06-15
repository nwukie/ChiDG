module type_element
    use mod_types,              only: rk,ik
    use mod_constants,          only: SPACEDIM,NFACES,XI_MIN,XI_MAX,ETA_MIN, &
                                      ETA_MAX,ZETA_MIN,ZETA_MAX,ONE,ZERO
    use mod_quadrature,         only: GQ,GQ_OVER,GQMESH,GQMESH_OVER
    use mod_element_mapping,    only: reference_mapping
    use mod_io,                 only: nterms_mesh3d
    use mod_polynomial,         only: polynomialVal
    use mod_inv,                only: inv

    use type_point,             only: point_t
    use type_variablesVector,   only: variablesVector_t
    use type_variables,         only: variables_t
    use type_sparserow,         only: sparserow_t
    use type_elementQuadrature, only: elementQuadrature_t


    implicit none

    ! Declare BLAS routines
    EXTERNAL    dgemv


    !=========================================================
    !=========================================================
    type, public :: element_t
        integer(kind=ik)                :: neqns
        integer(kind=ik)                :: nterms_sol
        integer(kind=ik)                :: nterms_mesh



        !=========================================================
        ! Element mesh points and modes
        !=========================================================
        type(point_t),      allocatable :: mesh_pts(:)
        type(variables_t)               :: mesh_modes       ! mesh_modes(nterms_var,(x,y,z))

        !=========================================================
        ! cartesian coordinates at quadrature values,
        ! mostly for solution initialization
        !=========================================================
        type(point_t),      allocatable :: quad_pts(:)
        type(point_t),      allocatable :: quad_pts_coll(:)


        !=========================================================
        ! jacobian terms at quadrature nodes
        !=========================================================
        real(kind=rk),      allocatable :: jinv(:)
        real(kind=rk),      allocatable :: jinv_over(:)

        !=========================================================
        ! metric matrix for each quadrature node    (mat_i,mat_j,quad_pt)
        !=========================================================
        real(kind=rk),      allocatable :: metric(:,:,:)
        real(kind=rk),      allocatable :: metric_over(:,:,:)


        !=========================================================
        !   Test function cartesian gradients
        !=========================================================
        real(kind=rk),      allocatable :: dtdx(:,:)
        real(kind=rk),      allocatable :: dtdy(:,:)
        real(kind=rk),      allocatable :: dtdz(:,:)

        !=========================================================
        ! Quadrature matrices
        !=========================================================
        type(elementQuadrature_t), pointer     :: gq_coll      => null()
        type(elementQuadrature_t), pointer     :: gqmesh_coll  => null()
        type(elementQuadrature_t), pointer     :: gq           => null()
        type(elementQuadrature_t), pointer     :: gq_over      => null()
        type(elementQuadrature_t), pointer     :: gqmesh       => null()
        type(elementQuadrature_t), pointer     :: gqmesh_over  => null()



        !=========================================================
        ! Element matrices
        !=========================================================
        real(kind=rk),      allocatable :: mass(:,:)
        real(kind=rk),      allocatable :: diagmass(:)
        real(kind=rk),      allocatable :: invmass(:,:)
        real(kind=rk),      allocatable :: Sxi(:,:)
        real(kind=rk),      allocatable :: Seta(:,:)
        real(kind=rk),      allocatable :: Szeta(:,:)

        !=========================================================
        ! Solution storage
        !=========================================================
        ! vector storage
        ! q_e
        ! var1_1    var2_1    ...   varN_1
        ! var1_2    var2_2    ...   varN_2
        ! var1_3    var2_3    ...   varN_3
        ! var1_...  var2_...  ...   varN_...

        ! q and rhs are always what the spatial discretization is working with
        ! _ref variables are used as reference states to reset values during numerical jacobian
        ! computation
        type(variables_t)               :: q
        type(variables_t)               :: q_ref    ! unperturbed solution

        type(variables_t)               :: rhs
        type(variables_t)               :: rhs_ref

        type(variables_t)               :: q_next

        type(variables_t)               :: dq
        type(variables_t)               :: dq_next

        ! Total lift mode storage for volume integrals
        type(variablesVector_t)         :: lift_g

        ! Local lift mode storage for each face
        type(variablesVector_t)         :: lift_l(NFACES)

        !=========================================================
        ! Storage of linearization for implicit solver
        !=========================================================
        type(sparserow_t)               :: linrow



    contains
        ! Initialization procedures
        procedure :: init
        procedure, public   :: compute_metrics
        procedure, public   :: compute_mass_matrix
        procedure, public   :: compute_stiffness_matrices
        procedure, public   :: initialize_solution

        ! Interpolation procedures
        procedure, public   :: compute_var


        ! Integration procedures
        procedure, public   :: integrate_volume_flux
        procedure, public   :: integrate_volume_source

        ! Utility procedures
        procedure, public   :: x
        procedure, public   :: y
        procedure, public   :: z
        procedure, public   :: solution_point


        final     :: destructor
    end type element_t
    !-----------------------
    private
contains
    
    subroutine init(self,neqns,nterms_sol,nterms_mesh,points)
        use mod_quadrature,     only: compute_nnodes_integration
        use mod_io,             only: integration_rule

        class(element_t),  intent(inout)      :: self
        type(point_t),     intent(in)         :: points(:)
        integer(kind=ik),  intent(in)         :: neqns
        integer(kind=ik),  intent(in)         :: nterms_sol
        integer(kind=ik),  intent(in)         :: nterms_mesh

        integer(kind=ik)                  :: AllocateStatus
        integer(kind=ik)                  :: nquad_nodes,nquad_nodes_over,nnodes_face,nnodes_vol, &
                                             nnodes_face_coll, nnodes_vol_coll,nquad_nodes_coll

        ! Set integer parameters
        self%neqns       = neqns
        self%nterms_sol  = nterms_sol
        self%nterms_mesh = nterms_mesh

        ! Set number of quadrature nodes
        call compute_nnodes_integration(1,nterms_sol,nnodes_face_coll,nnodes_vol_coll)
        call compute_nnodes_integration(integration_rule,nterms_sol,nnodes_face,nnodes_vol)
!        nquad_nodes      = nterms_sol
        nquad_nodes = nnodes_vol
        nquad_nodes_coll = nnodes_vol_coll
!        nquad_nodes = compute_overintegration(nterms_sol)
!        nquad_nodes_over = compute_overintegration(nterms_sol)

        ! Allocate element storage
        allocate(self%mesh_pts(nterms_mesh),stat=AllocateStatus)
        allocate(self%quad_pts(nquad_nodes),stat=AllocateStatus)
        allocate(self%quad_pts_coll(nquad_nodes_coll),stat=AllocateStatus)
        allocate(self%jinv(nquad_nodes),stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "Memory allocation error: element%init"

!        allocate(self%jinv_over(nquad_nodes_over),stat=AllocateStatus)
        allocate(self%metric(SPACEDIM,SPACEDIM,nquad_nodes),stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "Memory allocation error: element%init"

!        allocate(self%metric_over(SPACEDIM,SPACEDIM,nquad_nodes_over),stat=AllocateStatus)
        allocate(self%dtdx(nquad_nodes,nterms_sol))
        allocate(self%dtdy(nquad_nodes,nterms_sol))
        allocate(self%dtdz(nquad_nodes,nterms_sol))
        if (AllocateStatus /= 0) stop "Memory allocation error: element%init"


        allocate(self%mass(nterms_sol,nterms_sol),stat=AllocateStatus)
        allocate(self%diagmass(nterms_sol),stat=AllocateStatus)
        allocate(self%invmass(nterms_sol,nterms_sol),stat=AllocateStatus)
        allocate(self%Sxi(  nterms_sol,nterms_sol),stat=AllocateStatus)
        allocate(self%Seta( nterms_sol,nterms_sol),stat=AllocateStatus)
        allocate(self%Szeta(nterms_sol,nterms_sol),stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "Memory allocation error: element%init"


        ! Call initialization routines for derived type allocation
        call self%mesh_modes%init(nterms_mesh,SPACEDIM)  ! one set of modes for each coordinate

        ! Solution and right-hand side storage
        call self%q%init(      nterms_sol,neqns)
        call self%q_ref%init(  nterms_sol,neqns)
        call self%q_next%init( nterms_sol,neqns)

        call self%rhs%init(    nterms_sol,neqns)
        call self%rhs_ref%init(nterms_sol,neqns)

        ! 1D update vectors so they can be multiplied by the linearization matrices
        call self%dq%init(      neqns*nterms_sol,1)
        call self%dq_next%init( neqns*nterms_sol,1)

        ! global and face-local lifting mode storage
        call self%lift_g%init(nterms_sol,neqns)

        call self%lift_l(XI_MIN)%init(nterms_sol,neqns)
        call self%lift_l(XI_MAX)%init(nterms_sol,neqns)
        call self%lift_l(ETA_MIN)%init(nterms_sol,neqns)
        call self%lift_l(ETA_MAX)%init(nterms_sol,neqns)
        call self%lift_l(ZETA_MIN)%init(nterms_sol,neqns)
        call self%lift_l(ZETA_MAX)%init(nterms_sol,neqns)
        if (AllocateStatus /= 0) stop "Memory allocation error: element%init"



        ! The subblock sizes for the linearization storage will be
        ! square of size (nterms_sol*neqns  x  nterms_sol*neqns)
        call self%linrow%init(nterms_sol*neqns)

        ! Associate quadrature pointers to static objects so we can
        ! access them if we need to
        self%gq_coll     => GQ
        self%gqmesh_coll => GQMESH
        self%gq_over     => GQ_OVER
!        self%gqmesh      => GQMESH
!        self%gqmesh_over => GQMESH_OVER

        ! Associate gq with selected quadrature rule
        if (integration_rule == 1) then
            ! COLLOCATION QUADRATURE
            self%gq     =>  GQ
            self%gqmesh =>  GQMESH
        elseif (integration_rule == 2) then
            ! OVER-INTEGRATION RULE #1
            self%gq     =>  GQ_OVER
            self%gqmesh =>  GQMESH_OVER
        else
            print*, "Error: Element%init - Invalid integration rule"
            stop
        endif


        if (AllocateStatus /= 0) stop "Memory allocation error: element%init"

        ! Set element mesh points
        self%mesh_pts = points

        ! Compute polynomial representations of cartesian coordinates
        ! by multiplying the array of mesh coordinates by the inverted
        ! reference element mapping matrix
        self%mesh_modes%vals(:,1) = matmul(reference_mapping,self%mesh_pts(:)%x)
        self%mesh_modes%vals(:,2) = matmul(reference_mapping,self%mesh_pts(:)%y)
        self%mesh_modes%vals(:,3) = matmul(reference_mapping,self%mesh_pts(:)%z)


!        print*, "COMPUTING METRICS"
        call self%compute_metrics()

!        print*, "COMPUTING MASS MATRIX"
        call self%compute_mass_matrix()

!        print*, "COMPUTING STIFFNESS MATRIX"
        call self%compute_stiffness_matrices()
        call self%initialize_solution()

    end subroutine


    !============================================================
    !
    !  Convert local(xi,eta,zeta) coordinates to global coordinates(x,y,z)
    !
    !============================================================
    function x(self,xi,eta,zeta) result(xval)
        class(element_t),   intent(in)  :: self
        real(kind=rk),      intent(in)  :: xi,eta,zeta
        real(kind=rk)                   :: xval

        type(point_t)                   :: node
        real(kind=rk)                   :: polyvals(nterms_mesh3d)
        integer(kind=ik)                :: iterm

        node%x = xi
        node%y = eta
        node%z = zeta

        ! Evaluate polynomial modes at node location
        do iterm = 1,nterms_mesh3d
            polyvals(iterm)  = polynomialVal(3,nterms_mesh3d,iterm,node)
        end do

        ! Evaluate x from dot product of modes and polynomial values
        xval = dot_product(self%mesh_modes%vals(:,1),polyvals)

    end function

    function y(self,xi,eta,zeta) result(yval)
        class(element_t),   intent(in)  :: self
        real(kind=rk),      intent(in)  :: xi,eta,zeta
        real(kind=rk)                   :: yval

        type(point_t)                   :: node
        real(kind=rk)                   :: polyvals(nterms_mesh3d)
        integer(kind=ik)                :: iterm

        node%x = xi
        node%y = eta
        node%z = zeta

        ! Evaluate polynomial modes at node location
        do iterm = 1,nterms_mesh3d
            polyvals(iterm)  = polynomialVal(3,nterms_mesh3d,iterm,node)
        end do

        ! Evaluate x from dot product of modes and polynomial values
        yval = dot_product(self%mesh_modes%vals(:,2),polyvals)

    end function

    function z(self,xi,eta,zeta) result(zval)
        class(element_t),   intent(in)  :: self
        real(kind=rk),      intent(in)  :: xi,eta,zeta
        real(kind=rk)                   :: zval

        type(point_t)                   :: node
        real(kind=rk)                   :: polyvals(nterms_mesh3d)
        integer(kind=ik)                :: iterm

        node%x = xi
        node%y = eta
        node%z = zeta

        ! Evaluate polynomial modes at node location
        do iterm = 1,nterms_mesh3d
            polyvals(iterm)  = polynomialVal(3,nterms_mesh3d,iterm,node)
        end do

        ! Evaluate x from dot product of modes and polynomial values
        zval = dot_product(self%mesh_modes%vals(:,3),polyvals)

    end function

    function solution_point(self,varindex,xi,eta,zeta) result(val)
        class(element_t),   intent(in)  :: self
        integer(kind=ik),   intent(in)  :: varindex
        real(kind=rk),      intent(in)  :: xi,eta,zeta
        real(kind=rk)                   :: val

        type(point_t)                   :: node
        real(kind=rk)                   :: polyvals(self%nterms_sol)
        integer(kind=ik)                :: iterm

        node%x = xi
        node%y = eta
        node%z = zeta

        ! Evaluate polynomial modes at node location
        do iterm = 1,self%nterms_sol
            polyvals(iterm)  = polynomialVal(3,self%nterms_sol,iterm,node)
        end do

        ! Evaluate x from dot product of modes and polynomial values
        val = dot_product(self%q%vals(:,varindex),polyvals)

    end function

    !============================================================
    !
    !  Compute element metric and jacobian terms
    !
    !============================================================
    subroutine compute_metrics(self)
        class(element_t),    intent(inout)   :: self

        integer(kind=ik)    :: ixi,ieta,izeta
        integer(kind=ik)    :: inode,iterm

        real(kind=rk)       :: dxdxi(self%gq%vol%nnodes), dxdeta(self%gq%vol%nnodes), dxdzeta(self%gq%vol%nnodes)
        real(kind=rk)       :: dydxi(self%gq%vol%nnodes), dydeta(self%gq%vol%nnodes), dydzeta(self%gq%vol%nnodes)
        real(kind=rk)       :: dzdxi(self%gq%vol%nnodes), dzdeta(self%gq%vol%nnodes), dzdzeta(self%gq%vol%nnodes)
        real(kind=rk)       :: invjac(self%gq%vol%nnodes)
        real(kind=rk)       :: x(self%gq%vol%nnodes),y(self%gq%vol%nnodes),z(self%gq%vol%nnodes)

        ! compute element metric terms
        dxdxi   = matmul(self%gqmesh%vol%ddxi,  self%mesh_modes%vals(:,1))
        dxdeta  = matmul(self%gqmesh%vol%ddeta, self%mesh_modes%vals(:,1))
        dxdzeta = matmul(self%gqmesh%vol%ddzeta,self%mesh_modes%vals(:,1))


        dydxi   = matmul(self%gqmesh%vol%ddxi,  self%mesh_modes%vals(:,2))
        dydeta  = matmul(self%gqmesh%vol%ddeta, self%mesh_modes%vals(:,2))
        dydzeta = matmul(self%gqmesh%vol%ddzeta,self%mesh_modes%vals(:,2))


        dzdxi   = matmul(self%gqmesh%vol%ddxi,  self%mesh_modes%vals(:,3))
        dzdeta  = matmul(self%gqmesh%vol%ddeta, self%mesh_modes%vals(:,3))
        dzdzeta = matmul(self%gqmesh%vol%ddzeta,self%mesh_modes%vals(:,3))


        do inode = 1,self%gq%vol%nnodes
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


        ! compute cartesian coordinates associated with quadrature points
        x = matmul(self%gqmesh%vol%val,self%mesh_modes%vals(:,1))
        y = matmul(self%gqmesh%vol%val,self%mesh_modes%vals(:,2))
        z = matmul(self%gqmesh%vol%val,self%mesh_modes%vals(:,3))

        do inode = 1,self%gq%vol%nnodes
            call self%quad_pts(inode)%set_coord(x(inode),y(inode),z(inode))
        end do


        ! compute cartesian coordinates associated with collocation quadrature points
        x = matmul(self%gqmesh_coll%vol%val,self%mesh_modes%vals(:,1))
        y = matmul(self%gqmesh_coll%vol%val,self%mesh_modes%vals(:,2))
        z = matmul(self%gqmesh_coll%vol%val,self%mesh_modes%vals(:,3))

        do inode = 1,self%gq_coll%vol%nnodes
            call self%quad_pts_coll(inode)%set_coord(x(inode),y(inode),z(inode))
        end do


        ! Compute gradient of each test function term at quadrature nodes in cartesian coordinates
        ! Columns are quadrature points
        do iterm = 1,self%nterms_sol
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


    !============================================================
    !
    !  Initialize element solution
    !
    !============================================================
    subroutine initialize_solution(self)
        class(element_t),   intent(inout)   :: self

        real(kind=rk)   :: x,y,z
        real(kind=rk)   :: xval,yval,zval
        real(kind=rk), dimension(self%gq_coll%vol%nnodes)   :: &
                    avals, bvals, cvals, dvals, evals, init_vals

        integer(kind=ik)    :: inode

        ! Get quadrature x,y,z points and compute initialization function
        do inode = 1,self%gq_coll%vol%nnodes
            x = self%quad_pts_coll(inode)%x
            y = self%quad_pts_coll(inode)%y
            z = self%quad_pts_coll(inode)%z

!            xval = exp(-((x - 0._rk)**2._rk)/0.5_rk)
!            yval = exp(-((y - 2._rk)**2._rk)/0.5_rk)
!            zval = exp(-((z - 2._rk)**2._rk)/0.5_rk)

            xval = exp(-((x - 0._rk)**2._rk)/3.0_rk)
            yval = exp(-((y - 0._rk)**2._rk)/3.0_rk)
            zval = exp(-((z - 0._rk)**2._rk)/3.0_rk)

!             xval = 0._rk
!             yval = 0._rk
!             zval = 0._rk

            init_vals(inode) = xval*yval*zval
        end do

!        avals = 1.4_rk
!!        bvals = 10.0
!        bvals = ZERO
!        cvals = ZERO
!        dvals = ZERO
!        evals = 300000._rk

        avals = 1.23_rk
        bvals = 208._rk
        cvals = ZERO
        dvals = ZERO
        evals = 270470.0_rk


        ! Invert polynomial quadrature matrix and multiply by array of function values
!        self%q%vals(:,1) = matmul(inv(self%gq_coll%vol%val),init_vals)
        self%q%vals(:,1) = matmul(inv(self%gq_coll%vol%val),avals)
        self%q%vals(:,2) = matmul(inv(self%gq_coll%vol%val),bvals)
        self%q%vals(:,3) = matmul(inv(self%gq_coll%vol%val),cvals)
        self%q%vals(:,4) = matmul(inv(self%gq_coll%vol%val),dvals)
        self%q%vals(:,5) = matmul(inv(self%gq_coll%vol%val),evals)


    end subroutine


    !============================================================
    !
    !  Compute mass matrix
    !
    !============================================================
    subroutine compute_mass_matrix(self)
        class(element_t), intent(inout) :: self

        type(point_t)                      :: node
        real(kind=rk)                      :: weight
        type(elementQuadrature_t), pointer :: gq

        integer(kind=ik)                   :: iterm,i,j
        real(kind=rk)                      :: xi_val,eta_val,polyval_i,polyval_j,val,invjac

        real(kind=rk)                      :: temp(self%nterms_sol,self%gq%vol%nnodes)

        ! Compute mass matrix for current element
        self%invmass = 0._rk

        gq => self%gq


        ! First, just the polynomials and no jacobians, for projection routines
        temp = transpose(self%gq%vol%val)

        ! Multiply rows by quadrature weights and cell jacobians
        do iterm = 1,self%nterms_sol
            temp(iterm,:) = temp(iterm,:)*gq%vol%weights
        end do


        ! Perform the matrix multiplication of the transpose val matrix by
        ! the standard matrix. This produces the mass matrix. I think...
        self%mass = matmul(temp,gq%vol%val)


        ! Store the diagonal of the mass matrix
        do j = 1,self%nterms_sol
            do i = 1,self%nterms_sol

                if (i==j) then
                    self%diagmass(i) = self%mass(i,j)
                end if

            end do
        end do



        ! Redo, with jacobians for actuall mass matrix
        self%mass = 0._rk
        temp = transpose(self%gq%vol%val)

        ! Multiply rows by quadrature weights and cell jacobians
        do iterm = 1,self%nterms_sol
            temp(iterm,:) = temp(iterm,:)*gq%vol%weights*self%jinv
        end do

        ! Perform the matrix multiplication of the transpose val matrix by
        ! the standard matrix. This produces the mass matrix. I think...
        self%mass = matmul(temp,gq%vol%val)


        ! Store the inverted mass matrix
        self%invmass = inv(self%mass)



    end subroutine ! compute_mass_matrix



    !============================================================
    !
    !  Compute stiffness matrices
    !
    !============================================================
    subroutine compute_stiffness_matrices(self)
        class(element_t), intent(inout) :: self

        type(point_t)                      :: node
        real(kind=rk)                      :: weight
        type(elementQuadrature_t), pointer :: gq

        integer(kind=ik)                   :: iterm
        real(kind=rk)                      :: xi_val,eta_val,polyval_i,polyval_j,val,invjac

        real(kind=rk), dimension(self%nterms_sol,self%gq%vol%nnodes) :: &
                            tempx, tempy, tempz


        gq => self%gq


        ! Compute stiffness matrices for current element
        self%Sxi   = 0._rk
        self%Seta  = 0._rk
        self%Szeta = 0._rk



        tempx = transpose(self%dtdx)
        tempy = transpose(self%dtdy)
        tempz = transpose(self%dtdz)

        ! Multiply rows by quadrature weights and cell jacobians
        do iterm = 1,self%nterms_sol
            tempx(iterm,:) = tempx(iterm,:) * gq%vol%weights * self%jinv
            tempy(iterm,:) = tempy(iterm,:) * gq%vol%weights * self%jinv
            tempz(iterm,:) = tempz(iterm,:) * gq%vol%weights * self%jinv
        end do

        ! Perform the matrix multiplication of the transposed derivative matrix by
        ! the standard polynomial matrix. This produces the stiffness matrix.
        self%Sxi   = matmul(tempx,gq%vol%val)
        self%Seta  = matmul(tempy,gq%vol%val)
        self%Szeta = matmul(tempz,gq%vol%val)


    end subroutine ! compute_stiffness_matrices










    !=============================================================================
    !
    !   Integrate volume flux.
    !
    !   Multiplies volume fluxes by gradient of test functions and integrates.
    !
    !=============================================================================
    subroutine integrate_volume_flux(self,fluxx,fluxy,fluxz,varindex)
        class(element_t),   intent(inout)  :: self
        real(kind=rk),      intent(inout)  :: fluxx(:),fluxy(:),fluxz(:)
        integer(kind=ik),   intent(in)     :: varindex

        real(kind=rk),  dimension(self%nterms_sol)  ::   &
                        integralx, integraly, integralz, &
                        modesx,    modesy,    modesz
        integer(4)                            :: numrows,numcols




        ! So, we have the flux functions at quadrature points. Now, project
        ! to modal space.

        ! First, multiply by jacobian and quadrature weights

        ! I DON"T THINK WE NEED THE JACOBIANS IN THE PROJECTION
!        fluxx = fluxx * self%gq%vol%weights * self%jinv
!        fluxy = fluxy * self%gq%vol%weights * self%jinv
!        fluxz = fluxz * self%gq%vol%weights * self%jinv

        fluxx = fluxx * self%gq%vol%weights
        fluxy = fluxy * self%gq%vol%weights
        fluxz = fluxz * self%gq%vol%weights



        ! Project by integrating the function multiplied by a
        ! given mode and divide by the integral of the mode squared(diagonal of mass matrix).
        modesx = matmul(transpose(self%gq%vol%val),fluxx) / self%diagmass
        modesy = matmul(transpose(self%gq%vol%val),fluxy) / self%diagmass
        modesz = matmul(transpose(self%gq%vol%val),fluxz) / self%diagmass


        ! Contribute to rhs vector
        self%rhs%vals(:,varindex) = self%rhs%vals(:,varindex) + &
                                    matmul(self%Sxi,modesx)   + &
                                    matmul(self%Seta,modesy)  + &
                                    matmul(self%Szeta,modesz)




!        ! FROM HERE DOWN WAS PREVIOUS IMPLEMENTATION
!
!        numrows = size(self%dtdx,1)
!        numcols = size(self%dtdx,2)
!
!        ! Multiply each component by quadrature weights and jacobian
!        fluxx = fluxx * self%gq%vol%weights * self%jinv
!        fluxy = fluxy * self%gq%vol%weights * self%jinv
!        fluxz = fluxz * self%gq%vol%weights * self%jinv
!
!
!
!        ! Integrate flux
!        integralx = matmul(transpose(self%dtdx),fluxx)
!        integraly = matmul(transpose(self%dtdy),fluxy)
!        integralz = matmul(transpose(self%dtdz),fluxz)
!
!
!
!!        call dgemv('T',numrows,numcols,1.0_rk,self%dtdx,numrows,fluxx,1,0.0_rk,integralx,1)
!!        call dgemv('T',numrows,numcols,1.0_rk,self%dtdy,numrows,fluxy,1,0.0_rk,integraly,1)
!!        call dgemv('T',numrows,numcols,1.0_rk,self%dtdz,numrows,fluxz,1,0.0_rk,integralz,1)
!
!        ! Contribute to rhs vector
!        self%rhs%vals(:,varindex) = self%rhs%vals(:,varindex) + integralx + integraly + integralz

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
        real(kind=rk),      intent(inout)  :: source(:)
        integer(kind=ik),   intent(in)     :: varindex

        real(kind=rk),  dimension(self%nterms_sol)  :: integral
        integer(4) :: numrows,numcols

        numrows = size(self%gq%vol%val,1)
        numcols = size(self%gq%vol%val,2)

        ! Multiply by weights
        source = source * self%gq%vol%weights * self%jinv

        ! Integrate
        integral = matmul(transpose(self%gq%vol%val),source)


!        call dgemv('T',numrows,numcols,1.0_rk,self%gq%vol%val,numrows,source,1,0.0_rk,integral,1)

        ! Contribute to rhs vector
        self%rhs%vals(:,varindex) = self%rhs%vals(:,varindex) - integral


    end subroutine


    !=============================================================================
    !
    !
    !   Compute variable at quadrature nodes.
    !
    !
    !=============================================================================
    subroutine compute_var(self,varindex,vargq)
        class(element_t),   intent(in)      :: self
        integer(kind=ik),   intent(in)      :: varindex
        real(kind=rk),      intent(inout)   :: vargq(:)

        integer :: numrows,numcols

        numrows = size(self%gq%vol%val,1)
        numcols = size(self%gq%vol%val,2)

        ! Compute variables at volume GQ nodes
        vargq = matmul(self%gq%vol%val,self%q%vals(:,varindex))

!        call dgemv('N',numrows,numcols,1.0_rk,self%gq%vol%val,numrows,self%q%vals(:,varindex),1,0.0_rk,vargq,1)


    end subroutine


    subroutine destructor(self)
        type(element_t), intent(in) :: self
    end subroutine

end module type_element
