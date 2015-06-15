module type_face
    use mod_types,              only: rk,ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, &
                                      ZETA_MIN, ZETA_MAX, XI_DIR, ETA_DIR, ZETA_DIR, &
                                      SPACEDIM, NFACES
    use type_point,             only: point_t
    use type_element,           only: element_t
    use type_elementQuadrature, only: elementQuadrature_t
    use type_variables,         only: variables_t
    use type_variablesVector,   only: variablesVector_t

    implicit none

    ! Declare BLAS routines
    EXTERNAL    dgemv


    !------------------------------
    type, public :: face_t
        integer(kind=ik)                :: neqns
        integer(kind=ik)                :: nterms_face
        integer(kind=ik)                :: nterms_sol

        integer(kind=ik)                :: ftype           ! interior or boundary face
        integer(kind=ik)                :: iface           ! XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, etc
        type(point_t),    allocatable   :: mesh_pts(:)     ! points defining the face geometry

        type(variables_t), pointer      :: mesh_modes

        ! cartesian coordinates at quadrature values,
        ! mostly for boundary conditions
        type(point_t),    allocatable   :: quad_pts(:)


        real(kind=rk),    allocatable   :: jinv(:)
        real(kind=rk),    allocatable   :: metric(:,:,:)
        real(kind=rk),    allocatable   :: norm(:,:)

        real(kind=rk),    allocatable   :: massfaceref(:,:)
        real(kind=rk),    allocatable   :: massface(:,:,:)
        real(kind=rk),    allocatable   :: diagmassface(:)

        ! array containins lagrange polys and derivatives evaluated at the
        ! quadrature points. Used to resontstruct the solution by multiplying
        ! by the solution array
        type(elementQuadrature_t),  pointer     :: gq
        type(elementQuadrature_t),  pointer     :: gq_over
        type(elementQuadrature_t),  pointer     :: gqmesh
        type(elementQuadrature_t),  pointer     :: gqmesh_over

        ! pointers to the correct locations in the rhs and q vectors
        type(variables_t),          pointer     :: q
        type(variables_t),          pointer     :: q_ref
        type(variables_t),          pointer     :: q_next

        type(variables_t),          pointer     :: rhs
        type(variables_t),          pointer     :: rhs_ref

        type(variablesVector_t),    pointer     :: lift_g
        type(variablesVector_t),    pointer     :: lift_l(:)


        real(kind=rk),              pointer        :: invmass(:,:)  ! For BR2 lifting modes
    contains
        procedure   :: init
        procedure   :: compute_metrics
        procedure   :: compute_face_mass_matrix

        ! Interpolation to quadrature nodes
        procedure, public   :: compute_var

        ! Integration procedures
        procedure, public   :: integrate_flux
        procedure, public   :: integrate_scalar



        final       :: destructor
    end type face_t
    !------------------------------

    private
contains


    subroutine init(self,element,nterms_face,nterms_sol)
        class(face_t),    intent(inout)       :: self
        type(element_t),  intent(in), target  :: element
        integer(kind=ik), intent(in)          :: nterms_face,nterms_sol



        self%neqns       = element%neqns
        self%nterms_face = nterms_face
        self%nterms_sol  = nterms_sol
        self%mesh_modes  => element%mesh_modes

        self%q           => element%q
        self%q_ref       => element%q_ref
        self%q_next      => element%q_next

        self%rhs         => element%rhs
        self%rhs_ref     => element%rhs_ref

        self%lift_g      => element%lift_g
        self%lift_l      => element%lift_l

        self%invmass     => element%invmass

!        self%gq          => element%gq_coll
!        self%gqmesh      => element%gqmesh_coll
        self%gq          => element%gq
        self%gqmesh      => element%gqmesh
        self%gq_over     => element%gq_over
        self%gqmesh_over => element%gqmesh_over

        ! Allocate face storage
        allocate(self%quad_pts(self%gq%face%nnodes))
        allocate(self%jinv(self%gq%face%nnodes))
        allocate(self%metric(SPACEDIM,SPACEDIM,self%gq%face%nnodes))
        allocate(self%norm(self%gq%face%nnodes,SPACEDIM))

        allocate(self%massfaceref(nterms_sol,nterms_face))
        allocate(self%massface(nterms_sol,nterms_face,SPACEDIM))
        allocate(self%diagmassface(nterms_face))

        call self%compute_metrics()
        call self%compute_face_mass_matrix()

    end subroutine



    subroutine compute_metrics(self)
        class(face_t),  intent(inout)   :: self

        integer(kind=ik) :: inode, iface, nnodes

        real(kind=rk)    :: dxdxi(self%gq%face%nnodes), dxdeta(self%gq%face%nnodes), dxdzeta(self%gq%face%nnodes)
        real(kind=rk)    :: dydxi(self%gq%face%nnodes), dydeta(self%gq%face%nnodes), dydzeta(self%gq%face%nnodes)
        real(kind=rk)    :: dzdxi(self%gq%face%nnodes), dzdeta(self%gq%face%nnodes), dzdzeta(self%gq%face%nnodes)
        real(kind=rk)    :: invjac(self%gq%face%nnodes)
        real(kind=rk)    :: x(self%gq%face%nnodes),y(self%gq%face%nnodes),z(self%gq%face%nnodes)

        iface = self%iface
        nnodes = self%gq%face%nnodes

        associate (gq_f => self%gqmesh%face)
            dxdxi   = matmul(gq_f%ddxi(:,:,iface),  self%mesh_modes%vals(:,1))
            dxdeta  = matmul(gq_f%ddeta(:,:,iface), self%mesh_modes%vals(:,1))
            dxdzeta = matmul(gq_f%ddzeta(:,:,iface),self%mesh_modes%vals(:,1))

            dydxi   = matmul(gq_f%ddxi(:,:,iface),  self%mesh_modes%vals(:,2))
            dydeta  = matmul(gq_f%ddeta(:,:,iface), self%mesh_modes%vals(:,2))
            dydzeta = matmul(gq_f%ddzeta(:,:,iface),self%mesh_modes%vals(:,2))

            dzdxi   = matmul(gq_f%ddxi(:,:,iface),  self%mesh_modes%vals(:,3))
            dzdeta  = matmul(gq_f%ddeta(:,:,iface), self%mesh_modes%vals(:,3))
            dzdzeta = matmul(gq_f%ddzeta(:,:,iface),self%mesh_modes%vals(:,3))

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



        ! compute cartesian coordinates associated with quadrature points
        associate(gq_f => self%gqmesh%face)
            x = matmul(gq_f%val(:,:,iface),self%mesh_modes%vals(:,1))
            y = matmul(gq_f%val(:,:,iface),self%mesh_modes%vals(:,2))
            z = matmul(gq_f%val(:,:,iface),self%mesh_modes%vals(:,3))
        end associate

        do inode = 1,nnodes
            call self%quad_pts(inode)%set_coord(x(inode),y(inode),z(inode))
        end do

        ! Compute normal vectors
        select case (iface)
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
        if (iface == XI_MIN .or. iface == ETA_MIN .or. iface == ZETA_MIN) then
            self%norm(:,XI_DIR)   = -self%norm(:,XI_DIR)
            self%norm(:,ETA_DIR)  = -self%norm(:,ETA_DIR)
            self%norm(:,ZETA_DIR) = -self%norm(:,ZETA_DIR)
        end if


    end subroutine


    !============================================================================
    !
    !
    !   Compute face mass matrix
    !
    !
    !============================================================================
    subroutine compute_face_mass_matrix(self)
        class(face_t), intent(inout)       :: self

        type(elementQuadrature_t), pointer :: gq

        integer(kind=ik)                   :: iterm,iface,i,j
        real(kind=rk), dimension(self%nterms_sol,self%gq%face%nnodes) ::  &
                            temp_xi, temp_eta, temp_zeta, temp_ref

        real(kind=rk), dimension(self%nterms_face,self%gq%face%nnodes) :: temp

        real(kind=rk), dimension(self%nterms_face,self%nterms_face) :: tempmass

        iface = self%iface

        gq => self%gq


!        temp      = transpose(gq%face%faceval(:,:,iface))
!        temp_xi   = transpose(gq%face%faceval(:,:,iface))
!        temp_eta  = transpose(gq%face%faceval(:,:,iface))
!        temp_zeta = transpose(gq%face%faceval(:,:,iface))

        temp      = transpose(gq%face%faceval(:,:,iface))

        temp_xi   = transpose(gq%face%val(:,:,iface))
        temp_eta  = transpose(gq%face%val(:,:,iface))
        temp_zeta = transpose(gq%face%val(:,:,iface))
        temp_ref  = transpose(gq%face%val(:,:,iface))


        ! Multiply rows by quadrature weights and face jacobian/normals
        do iterm = 1,self%nterms_sol
            temp_ref(iterm,:)  = temp_ref(iterm,:)  * gq%face%weights(:,iface)
            temp_xi(iterm,:)   = temp_xi(iterm,:)   * gq%face%weights(:,iface) * self%norm(:,1)
            temp_eta(iterm,:)  = temp_eta(iterm,:)  * gq%face%weights(:,iface) * self%norm(:,2)
            temp_zeta(iterm,:) = temp_zeta(iterm,:) * gq%face%weights(:,iface) * self%norm(:,3)
        end do

        do iterm = 1,self%nterms_face
            temp(iterm,:)      = temp(iterm,:)      * gq%face%weights(:,iface)
        end do
        tempmass = matmul(temp,gq%face%faceval(:,:,iface))

        ! Perform the matrix multiplication of the transpose of the face expansion
        ! with the volume expansion at the face. This produces the face mass matrix.
        self%massfaceref     = matmul(temp_ref, gq%face%faceval(:,:,iface))
        self%massface(:,:,1) = matmul(temp_xi,  gq%face%faceval(:,:,iface))
        self%massface(:,:,2) = matmul(temp_eta, gq%face%faceval(:,:,iface))
        self%massface(:,:,3) = matmul(temp_zeta,gq%face%faceval(:,:,iface))

        ! Store the diagonal of the mass matrix for projections
        do j = 1,self%nterms_face
            do i = 1,self%nterms_face

                if (i==j) then
                    self%diagmassface(i) = 1._rk/tempmass(i,j)
                end if

            end do
        end do

    end subroutine

    !=============================================================================
    !
    !   Integrate boundary flux.
    !
    !   Multiplies volume fluxes by gradient of test functions and integrates.
    !
    !=============================================================================
    subroutine integrate_flux(self,fluxx,fluxy,fluxz,varindex)
        class(face_t),      intent(inout)  :: self
        real(kind=rk),      intent(inout)  :: fluxx(:), fluxy(:), fluxz(:)
        integer(kind=ik),   intent(in)     :: varindex

        real(kind=rk),  dimension(self%nterms_face)  :: &
                          integral, modesx, modesy, modesz

        integer(4)                            :: iface,numrows,numcols


        iface = self%iface

        ! So, we have the numerical flux function at quadrature nodes. Now, need to
        ! project that information to face polynomial expansion.

        ! Multiply by quadrature weights on the face
        fluxx = fluxx * self%gq%face%weights(:,iface)
        fluxy = fluxy * self%gq%face%weights(:,iface)
        fluxz = fluxz * self%gq%face%weights(:,iface)

        ! Project by integrating the function multiplied by
        ! given mode and divied by the integral of the mode squared.
!        modesx = matmul(transpose(self%gq%face%faceval(:,:,iface)),fluxx) / self%diagmassface
!        modesy = matmul(transpose(self%gq%face%faceval(:,:,iface)),fluxy) / self%diagmassface
!        modesz = matmul(transpose(self%gq%face%faceval(:,:,iface)),fluxz) / self%diagmassface

        modesx = matmul(transpose(self%gq%face%faceval(:,:,iface)),fluxx) * self%diagmassface
        modesy = matmul(transpose(self%gq%face%faceval(:,:,iface)),fluxy) * self%diagmassface
        modesz = matmul(transpose(self%gq%face%faceval(:,:,iface)),fluxz) * self%diagmassface


        ! Contribute to rhs vector
        self%rhs%vals(:,varindex) = self%rhs%vals(:,varindex) - &
                                        matmul(self%massface(:,:,1),modesx) - &
                                        matmul(self%massface(:,:,2),modesy) - &
                                        matmul(self%massface(:,:,3),modesz)


!        ! BELOW HERE IS PREVIOUS IMPLEMENTATION
!
!        iface = self%iface
!
!        numrows = size(self%gq%face%val(:,:,iface),1)
!        numcols = size(self%gq%face%val(:,:,iface),2)
!
!
!        ! Multiply by quadrature weights
!        flux = flux  * self%gq%face%weights(:,iface)
!
!
!        ! Multiply by column of test functions and integrate
!        integral = matmul(transpose(self%gq%face%val(:,:,iface)),flux)
!
!!        call dgemv('T',numrows,numcols,1.0_rk,self%gq%face%val(:,:,iface),numrows,flux,1,0.0_rk,integral,1)
!
!
!        ! Set face RHS
!        self%rhs%vals(:,varindex) = self%rhs%vals(:,varindex) - integral

    end subroutine








    !=============================================================================
    !
    !   Integrate scalar flux.
    !
    !
    !
    !=============================================================================
    subroutine integrate_scalar(self,scalar,varindex)
        class(face_t),      intent(inout)  :: self
        real(kind=rk),      intent(inout)  :: scalar(:)
        integer(kind=ik),   intent(in)     :: varindex

        real(kind=rk),  dimension(self%nterms_face)  :: modes

        integer(4)                            :: iface,numrows,numcols


        iface = self%iface
        numrows = size(self%gq%face%faceval(:,:,iface),1)
        numcols = size(self%gq%face%faceval(:,:,iface),2)

        ! So, we have the numerical flux function at quadrature nodes. Now, need to
        ! project that information to face polynomial expansion.

        ! Multiply by quadrature weights on the face
        scalar = scalar * self%gq%face%weights(:,iface)

        ! Project
!        modes = matmul(transpose(self%gq%face%faceval(:,:,iface)),scalar) * self%diagmassface
        call dgemv('T',numrows,numcols,self%diagmassface,self%gq%face%faceval(:,:,iface),numrows,scalar,1,0.0_rk,modes,1)


        ! Contribute to rhs vector
        self%rhs%vals(:,varindex) = self%rhs%vals(:,varindex) - &
                                        matmul(self%massfaceref,modes)


!        ! BELOW HERE IS PREVIOUS IMPLEMENTATION
!
!        iface = self%iface
!
!        numrows = size(self%gq%face%val(:,:,iface),1)
!        numcols = size(self%gq%face%val(:,:,iface),2)
!
!
!        ! Multiply by quadrature weights
!        flux = flux  * self%gq%face%weights(:,iface)
!
!
!        ! Multiply by column of test functions and integrate
!        integral = matmul(transpose(self%gq%face%val(:,:,iface)),flux)
!
!!        call dgemv('T',numrows,numcols,1.0_rk,self%gq%face%val(:,:,iface),numrows,flux,1,0.0_rk,integral,1)
!
!
!        ! Set face RHS
!        self%rhs%vals(:,varindex) = self%rhs%vals(:,varindex) - integral

    end subroutine
















    !=============================================================================
    !
    !
    !   Compute variable at quadrature nodes.
    !
    !
    !=============================================================================
    subroutine compute_var(self,varindex,vargq)
        class(face_t),      intent(in)      :: self
        integer(kind=ik),   intent(in)      :: varindex
        real(kind=rk),      intent(inout)   :: vargq(:)

        integer(kind=ik)                    :: iface
        integer(4)                          :: numrows, numcols

        iface = self%iface


        ! Compute variables at volume GQ nodes
!        vargq = matmul(self%gq%face%val(:,:,iface),self%q%vals(:,varindex))



        numrows = size(self%gq%face%val(:,:,iface),1)
        numcols = size(self%gq%face%val(:,:,iface),2)
        call dgemv('N',numrows,numcols,1.0_rk,self%gq%face%val(:,:,iface),numrows,self%q%vals(:,varindex),1,0.0_rk,vargq,1)


    end subroutine










    subroutine destructor(self)
        type(face_t), intent(in) :: self
    end subroutine

end module type_face
