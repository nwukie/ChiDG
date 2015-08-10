module type_mesh
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: NFACES,XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX

    use type_element,       only: element_t
    use type_face,          only: face_t
    use type_point,         only: point_t

    implicit none
    private


    !> Data type for mesh information
    type, public :: mesh_t
        ! Integer parameters
        integer(ik)         :: neqns    = 0     !> Number of equations being solved
        integer(ik)         :: nterms_s = 0     !> Number of terms in the solution expansion
        integer(ik)         :: nterms_c = 0     !> Number of terms in the grid coordinate expansion
        integer(ik)         :: nelem_xi, nelem_eta, nelem_zeta, nelem

        ! Grid data
        type(element_t),  allocatable  :: elems(:)                  !> Element storage (1:nelem)
        type(element_t),  pointer      :: elems_m(:,:,:) => null()  !> Matrix view of element storage (1:nelem_xi, 1:nelem_eta, 1:nelem_zeta)
        type(face_t),     allocatable  :: faces(:,:)                !> Face storage    (1:nelem,1:nfaces)

        ! TODO: Needs tested
!        type(face_t),     pointer      :: faces_c(:,:,:,:) => null()    !> Matrix view of face storage

    contains
        procedure           :: init_geom
        procedure           :: init_sol

        procedure, private  :: init_elems_geom
        procedure, private  :: init_elems_sol
        procedure, private  :: init_faces_geom
        procedure, private  :: init_faces_sol

        procedure, private  :: init_boundary_conditions

        final :: destructor
    end type mesh_t

contains

    !> Mesh geometry initialization procedure
    !!
    !!  Sets number of terms in coordinate expansion for the entire domain
    !!  and calls sub-initialization routines for individual element and face geometry
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]  nterms_c    Number of terms in the coordinate expansion
    !!  @param[in]  points_g    Rank-3 matrix of coordinate points defining a block mesh
    !---------------------------------------------------------------------------------------
    subroutine init_geom(self,nterms_c,points_g)
        class(mesh_t),  intent(inout), target   :: self
        integer(ik),    intent(in)              :: nterms_c
        type(point_t),  intent(in)              :: points_g(:,:,:)
        type(element_t), pointer                :: temp(:)
!        type(face_t),    pointer                :: ftemp(:,:)

        self%nterms_c = nterms_c

        call self%init_elems_geom(points_g)
        call self%init_faces_geom()

        ! Initialize element matrix view
        temp => self%elems
        self%elems_m(1:self%nelem_xi,1:self%nelem_eta,1:self%nelem_zeta) => temp(1:self%nelem)

!        TODO: NEEDS TESTED
!        ftemp => self%faces
!        self%faces_m(1:self%nelem_xi,1:self%nelem_eta,1:self%nelem_zeta,NFACES) => ftemp(1:self%nelem,NFACES)

        !> Initialize boundary conditions after geometry is set up
        call self%init_boundary_conditions()
    end subroutine


    !> Mesh numerics initialization procedure
    !!
    !!  Sets number of equations being solved, number of terms in the solution expansion and
    !!  calls sub-initialization routines for individual element and face numerics
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]  neqns       Number of equations being solved in the current domain
    !!  @param[in]  nterms_s    Number of terms in the solution expansion
    !---------------------------------------------------------------------------------------
    subroutine init_sol(self,neqns,nterms_s)
        class(mesh_t),  intent(inout)   :: self
        integer(ik),    intent(in)      :: neqns
        integer(ik),    intent(in)      :: nterms_s

        self%neqns    = neqns
        self%nterms_s = nterms_s

        call self%init_elems_sol(neqns,nterms_s)
        call self%init_faces_sol()

    end subroutine





    !> Mesh - element initialization procedure
    !!
    !!  Computes the number of elements based on the element mapping selected and
    !!  calls the element initialization procedure on individual elements.
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]  points_g    Rank-3 matrix of coordinate points defining a block mesh
    !-------------------------------------------------------------------------------------------
    subroutine init_elems_geom(self,points_g)
        class(mesh_t),  intent(inout)   :: self
        type(point_t),  intent(in)      :: points_g(:,:,:)
        type(point_t),  allocatable     :: points_l(:)

        integer(ik)                ::   ierr,     ipt,       ielem,             &
                                        ipt_xi,   ipt_eta,   ipt_zeta,          &
                                        ixi,      ieta,      izeta,             &
                                        npts_xi,  npts_eta,  npts_zeta,         &
                                        xi_start, eta_start, zeta_start,        &
                                        nelem_xi, nelem_eta, nelem_zeta, nelem, &
                                        neqns,    nterms_s,  nnodes, nterms_c, npts_1d, mapping

        npts_xi   = size(points_g,1)    !> Number of points in the xi-direction
        npts_eta  = size(points_g,2)    !> Number of points in the eta-direction
        npts_zeta = size(points_g,3)    !> Number of points in the zeta-direction

        ! Compute number of 1d points for a single element
        npts_1d = 0
        do while (npts_1d*npts_1d*npts_1d < self%nterms_c)
            npts_1d = npts_1d + 1       ! really just computing the cubed root of nterms_c, the number of terms in the coordinate expansion
        end do



        ! Count number of elements in each direction and check mesh conforms to
        ! the agglomeration rule for higher-order elements
        nelem_xi = 0
        ipt = 1
        do while (ipt < npts_xi)
            nelem_xi = nelem_xi + 1
            ipt = ipt + (npts_1d-1)
        end do
        if (ipt > npts_xi) stop "Mesh does not conform to agglomeration routine in xi"

        nelem_eta = 0
        ipt = 1
        do while (ipt < npts_eta)
            nelem_eta = nelem_eta + 1
            ipt = ipt + (npts_1d-1)
        end do
        if (ipt > npts_eta) stop "Mesh does not conform to agglomeration routine in eta"

        nelem_zeta = 0
        ipt = 1
        do while (ipt < npts_zeta)
            nelem_zeta = nelem_zeta + 1
            ipt = ipt + (npts_1d-1)
        end do
        if (ipt > npts_zeta) stop "Mesh does not conform to agglomeration routine in zeta"


!        ! Print mesh characteristics
!        print*, "       ... npts_xi, npts_eta, npts_zeta"
!        print*, npts_xi, npts_eta, npts_zeta
!
!        print*, "       ... mesh mapping"
!        select case (npts_1d-1)
!            case (1)
!                print*, 'Linear'
!            case (2)
!                print*, 'Quadratic'
!            case (3)
!                print*, 'Cubic'
!            case (4)
!                print*, 'Quartic'
!            case default
!                stop "Error: mesh%init - Invalid element mapping"
!        end select
!
!        print*, "       ... nelem_xi, nelem_eta, nelem_zeta"
!        print*, nelem_xi, nelem_eta, nelem_zeta


        self%nelem_xi   = nelem_xi
        self%nelem_eta  = nelem_eta
        self%nelem_zeta = nelem_zeta
        nelem           = nelem_xi * nelem_eta * nelem_zeta
        self%nelem      = nelem
        mapping         = (npts_1d - 1)     !> 1 - linear, 2 - quadratic, 3 - cubic, etc.

        ! Allocate element storage
        allocate(self%elems(nelem),stat=ierr)
        allocate(points_l(self%nterms_c))
        if(ierr /= 0) stop "Memory allocation error: init_elements"

        ielem = 1
        ! Initialize elements
        do izeta = 1,nelem_zeta
            do ieta = 1,nelem_eta
                do ixi = 1,nelem_xi

                    xi_start   = 1 + (ixi  -1)*(npts_1d-1)
                    eta_start  = 1 + (ieta -1)*(npts_1d-1)
                    zeta_start = 1 + (izeta-1)*(npts_1d-1)
                    ! For this element, collect the necessary points from the global points
                    ! array into a local points array for initializing an individual element
                    ipt = 1
                    do ipt_zeta = 1,npts_1d
                        do ipt_eta = 1,npts_1d
                            do ipt_xi = 1,npts_1d
                                points_l(ipt) = points_g(xi_start+(ipt_xi-1),eta_start+(ipt_eta-1),zeta_start+(ipt_zeta-1))
                                ipt = ipt + 1
                            end do
                        end do
                    end do

                    call self%elems(ielem)%init_geom(mapping,points_l,ielem)
                    ielem = ielem + 1
                end do
            end do
        end do
    end subroutine




    subroutine init_elems_sol(self,neqns,nterms_s)
        class(mesh_t),  intent(inout)   :: self
        integer(ik),    intent(in)      :: neqns
        integer(ik),    intent(in)      :: nterms_s
        integer(ik) :: ielem

        self%neqns    = neqns
        self%nterms_s = nterms_s

        do ielem = 1,self%nelem
            call self%elems(ielem)%init_sol(self%neqns,self%nterms_s)
        end do
    end subroutine


    !> Mesh - face initialization procedure
    !!
    !!  @author Nathan A. Wukie
    !!
    !-----------------------------------------------------------------------------------
    subroutine init_faces_geom(self)
        class(mesh_t), intent(inout)  :: self

        integer(ik)              :: ixi,ieta,izeta,iface,ftype,ineighbor,ielem,ierr

        ! Allocate face storage
        allocate(self%faces(self%nelem,NFACES),stat=ierr)
        if (ierr /= 0) stop "Error: mesh%init -- face allocation error"

        ielem = 1
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
                            ftype = 1       ! boundary face
                            ineighbor = 0   ! No neighbor
                        else
                            ftype = 0  ! interior face

                            select case (iface)
                                case (XI_MIN)
                                    ineighbor = ielem - 1
                                case (XI_MAX)
                                    ineighbor = ielem + 1
                                case (ETA_MIN)
                                    ineighbor = ielem - self%nelem_xi
                                case (ETA_MAX)
                                    ineighbor = ielem + self%nelem_xi
                                case (ZETA_MIN)
                                    ineighbor = ielem - self%nelem_xi*self%nelem_eta
                                case (ZETA_MAX)
                                    ineighbor = ielem + self%nelem_xi*self%nelem_eta
                            end select

                        end if


                        ! Call face initialization routine
                        call self%faces(ielem,iface)%init_geom(iface,ftype,self%elems(ielem),ineighbor)

                    end do

                    ielem = ielem + 1
                end do !ixi
            end do ! ieta
        end do ! izeta

    end subroutine




    !> Mesh - face initialization procedure
    !!
    !!  @author Nathan A. Wukie
    !!
    !-----------------------------------------------------------------------------------
    subroutine init_faces_sol(self)
        class(mesh_t), intent(inout)  :: self

        integer(ik) :: ielem, iface

        do ielem = 1,self%nelem
            do iface = 1,NFACES
                call self%faces(ielem,iface)%init_sol(self%elems(ielem))
            end do
        end do
    end subroutine








    !> Initialize boundary conditions
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------
    subroutine init_boundary_conditions(self)
        class(mesh_t),  intent(inout)   :: self

        integer(ik) :: ielem, ielem_p, ixi, ieta, izeta, ixi_p, ieta_p, izeta_p



        ! Apply periodic XI
        ixi = 1
        ixi_p = self%nelem_xi
        do izeta = 1,self%nelem_zeta
            do ieta = 1,self%nelem_eta
                ielem = self%elems_m(ixi,ieta,izeta)%ielem
                ielem_p = self%elems_m(ixi_p,ieta,izeta)%ielem

                self%faces(ielem,XI_MIN)%ftype = 0              !> Interior face
                self%faces(ielem,XI_MIN)%ineighbor = ielem_p    !> Set neighbor face to be periodic

                self%faces(ielem_p,XI_MAX)%ftype = 0
                self%faces(ielem_p,XI_MAX)%ineighbor = ielem
            end do
        end do


        ! Apply periodic ETA
        ieta = 1
        ieta_p = self%nelem_eta
        do izeta = 1,self%nelem_zeta
            do ixi = 1,self%nelem_xi
                ielem = self%elems_m(ixi,ieta,izeta)%ielem
                ielem_p = self%elems_m(ixi,ieta_p,izeta)%ielem

                self%faces(ielem,ETA_MIN)%ftype = 0              !> Interior face
                self%faces(ielem,ETA_MIN)%ineighbor = ielem_p    !> Set neighbor face to be periodic

                self%faces(ielem_p,ETA_MAX)%ftype = 0
                self%faces(ielem_p,ETA_MAX)%ineighbor = ielem
            end do
        end do


        ! Apply periodic ZETA
        izeta = 1
        izeta_p = self%nelem_zeta
        do ieta = 1,self%nelem_zeta
            do ixi = 1,self%nelem_xi
                ielem = self%elems_m(ixi,ieta,izeta)%ielem
                ielem_p = self%elems_m(ixi,ieta,izeta_p)%ielem

                self%faces(ielem,ZETA_MIN)%ftype = 0              !> Interior face
                self%faces(ielem,ZETA_MIN)%ineighbor = ielem_p    !> Set neighbor face to be periodic

                self%faces(ielem_p,ZETA_MAX)%ftype = 0
                self%faces(ielem_p,ZETA_MAX)%ineighbor = ielem
            end do
        end do





    end subroutine



















    subroutine destructor(self)
        type(mesh_t), intent(inout) :: self

!> Shouldn't need to deallocate an 'allocatable'. The compiler is supposed to do that for you
!        if (allocated(self%elems)) deallocate(self%elems)
!        if (allocated(self%faces)) deallocate(self%faces)
    end subroutine

end module type_mesh
