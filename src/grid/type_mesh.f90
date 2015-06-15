module type_mesh
    use mod_types,          only: rk,ik
    use mod_constants,      only: NFACES,XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX
    use mod_io,             only: ntime_instances, mesh_poly_order,nterms_mesh3d

    use type_element,       only: element_t
    use type_face,          only: face_t
    use type_point,         only: point_t

    implicit none
    private


    !> Data type for mesh information
    type, public :: mesh_t
        ! Integer parameters
        integer(kind=ik)                    :: neqns
        integer(kind=ik)                    :: ntime_instances

        integer(kind=ik)                    :: nelem_xi, nelem_eta, nelem_zeta, nelem

        ! Grid data
        type(element_t), dimension(:,:,:),      allocatable     :: elements
        type(face_t),    dimension(:,:,:,:),    allocatable     :: faces


    contains
        procedure, public   :: init
        procedure, private  :: init_elements
        procedure, private  :: init_faces

        final :: destructor
    end type mesh_t

contains

    !============================================================
    !
    ! Mesh Initialization
    !
    ! Sets integer values and calls sub initialization routines
    ! for elements, faces, and boundary conditions
    !============================================================
    subroutine init(self,neqns,nterms_sol,nterms_mesh,points_g)
        class(mesh_t),      intent(inout)  :: self
        integer(kind=ik),   intent(in)     :: neqns
        integer(kind=ik),   intent(in)     :: nterms_sol,nterms_mesh
        type(point_t),      intent(in)     :: points_g(:,:,:)
        integer(kind=ik)                   :: nterms

        self%neqns           = neqns
        self%nterms_sol3d    = nterms_sol
        self%nterms_mesh3d   = nterms_mesh
        self%ntime_instances = ntime_instances

        ! compute number of 1d terms
        nterms = 0
        do while (nterms*nterms*nterms /= self%nterms_sol3d)
            nterms = nterms+1
        end do
        self%nterms_sol1d = nterms
        self%nterms_sol2d = nterms*nterms

        nterms = 0
        do while (nterms*nterms*nterms /= self%nterms_mesh3d)
            nterms = nterms+1
        end do
        self%nterms_mesh1d = nterms


        call self%init_elements(points_g)
        call self%init_faces()



        print*, "Mesh initialization completed!"
    end subroutine
    

    !============================================================
    !
    !  Element Initialization
    !
    !============================================================
    subroutine init_elements(self,points_g)
        class(mesh_t),  intent(inout)   :: self
        type(point_t),  intent(in)      :: points_g(:,:,:)
        type(point_t),  allocatable     :: points_l(:)

        integer(kind=ik)                :: AllocateStatus
        integer(kind=ik)                :: npts_xi,npts_eta,npts_zeta
        integer(kind=ik)                :: ipt_xi,ipt_eta,ipt_zeta,ipt
        integer(kind=ik)                :: xi_start,eta_start,zeta_start
        integer(kind=ik)                :: nelem_xi,nelem_eta,nelem_zeta
        integer(kind=ik)                :: neqns,nterms_sol,nquad_nodes,nquad_nodes_over,nterms_mesh
        integer(kind=ik)                :: ixi,ieta,izeta

        npts_xi   = size(points_g,1)
        npts_eta  = size(points_g,2)
        npts_zeta = size(points_g,3)

        ! Count number of elements in each direction
        nelem_xi = 0
        ipt_xi = 1
        do while (ipt_xi < npts_xi)
            nelem_xi = nelem_xi + 1
            ipt_xi = ipt_xi + (self%nterms_mesh1d-1)
        end do
        if (ipt_xi > npts_xi) stop "Mesh does not conform to agglomeration routine in xi"

        nelem_eta = 0
        ipt_eta = 1
        do while (ipt_eta < npts_eta)
            nelem_eta = nelem_eta + 1
            ipt_eta = ipt_eta + (self%nterms_mesh1d-1)
        end do
        if (ipt_eta > npts_eta) stop "Mesh does not conform to agglomeration routine in eta"

        nelem_zeta = 0
        ipt_zeta = 1
        do while (ipt_zeta < npts_zeta)
            nelem_zeta = nelem_zeta + 1
            ipt_zeta = ipt_zeta + (self%nterms_mesh1d-1)
        end do
        if (ipt_zeta > npts_zeta) stop "Mesh does not conform to agglomeration routine in zeta"

        print*, "       ... npts_xi, npts_eta, npts_zeta"
        print*, npts_xi, npts_eta, npts_zeta

        print*, "       ... mesh polynomial order"
        print*, mesh_poly_order

        print*, "       ... nelem_xi, nelem_eta, nelem_zeta"
        print*, nelem_xi, nelem_eta, nelem_zeta

        self%nelem_xi   = nelem_xi
        self%nelem_eta  = nelem_eta
        self%nelem_zeta = nelem_zeta

        ! Allocate element storage
        allocate(self%elements(nelem_xi,nelem_eta,nelem_zeta),stat=AllocateStatus)
        allocate(points_l(self%nterms_mesh3d))
        if(AllocateStatus /= 0) stop "Memory allocation error: init_elements"

        ! Initialize elements
        do izeta = 1,nelem_zeta
            do ieta = 1,nelem_eta
                do ixi = 1,nelem_xi

                    xi_start   = 1 + (ixi  -1)*(self%nterms_mesh1d-1)
                    eta_start  = 1 + (ieta -1)*(self%nterms_mesh1d-1)
                    zeta_start = 1 + (izeta-1)*(self%nterms_mesh1d-1)
                    ! For this element, collect the necessary points from the global points
                    ! array into a local points array for initialization
                    ipt = 1
                    do ipt_zeta = 1,self%nterms_mesh1d
                        do ipt_eta = 1,self%nterms_mesh1d
                            do ipt_xi = 1,self%nterms_mesh1d
                                points_l(ipt) = points_g(xi_start+(ipt_xi-1),eta_start+(ipt_eta-1),zeta_start+(ipt_zeta-1))
                                ipt = ipt + 1
                            end do
                        end do
                    end do

                    call self%elements(ixi,ieta,izeta)%init(self%neqns,self%nterms_sol3d,self%nterms_mesh3d,points_l)
                end do
            end do
        end do
    end subroutine


    !============================================================
    !
    !  Face Initialization
    !
    !============================================================
    subroutine init_faces(self)
        class(mesh_t), intent(inout)  :: self

        integer(kind=ik)              :: nterms_face,nterms_sol
        integer(kind=ik)              :: ixi,ieta,izeta,iface

        ! Allocate face storage
        allocate(self%faces(self%nelem_xi,self%nelem_eta,self%nelem_zeta,NFACES))
        nterms_face = self%nterms_sol2d
        nterms_sol  = self%nterms_sol3d

        do izeta = 1,self%nelem_zeta
            do ieta = 1,self%nelem_eta
                do ixi = 1,self%nelem_xi

                    do iface = 1,NFACES

                        ! Set ftype to designate interior and boundary faces
                        if ( (ixi == 1                 .and. iface == XI_MIN)   .or. &
                             (ixi == self%nelem_xi     .and. iface == XI_MAX)   .or. &
                             (ieta == 1                .and. iface == ETA_MIN)  .or. &
                             (ieta == self%nelem_eta   .and. iface == ETA_MAX)  .or. &
                             (izeta == 1               .and. iface == ZETA_MIN) .or. &
                             (izeta == self%nelem_zeta .and. iface == ZETA_MAX) ) then
                            self%faces(ixi,ieta,izeta,iface)%ftype = 1  ! boundary face
                        else
                            self%faces(ixi,ieta,izeta,iface)%ftype = 0  ! interior face
                        end if

                        ! Set face index - XI_MIN,XI_MAX,ETA_MIN, eta.
                        self%faces(ixi,ieta,izeta,iface)%iface = iface

                        ! Call face initialization routine
                        call self%faces(ixi,ieta,izeta,iface)%init(self%elements(ixi,ieta,izeta),nterms_face,nterms_sol)

                    end do

                end do !ixi
            end do ! ieta
        end do ! izeta

    end subroutine






    subroutine destructor(self)
        type(mesh_t), intent(in) :: self
    end subroutine

end module type_mesh
