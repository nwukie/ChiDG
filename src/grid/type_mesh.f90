module type_mesh
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: NFACES,XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX, &
                                  ORPHAN, INTERIOR, BOUNDARY, CHIMERA

    use type_element,       only: element_t
    use type_face,          only: face_t
    use type_point,         only: point_t

    use type_chimera,       only: chimera_t
    implicit none
    private


    !> Data type for mesh information
    !!      - contains array of elements, array of faces for each element
    !!      - calls initialization procedure for elements and faces
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !------------------------------------------------------------------------------------------------------------
    type, public :: mesh_t
        ! Integer parameters
        integer(ik)         :: neqns      = 0                           !< Number of equations being solved
        integer(ik)         :: nterms_s   = 0                           !< Number of terms in the solution expansion
        integer(ik)         :: nterms_c   = 0                           !< Number of terms in the grid coordinate expansion
        integer(ik)         :: nelem_xi   = 0                           !< Number of elements in the xi-direction
        integer(ik)         :: nelem_eta  = 0                           !< Number of elements in the eta-direction
        integer(ik)         :: nelem_zeta = 0                           !< Number of elements in the zeta-direction
        integer(ik)         :: nelem      = 0                           !< Number of total elements

        ! Grid data
        integer(ik)                     :: idomain
        type(element_t),  allocatable   :: elems(:)                     !< Element storage (1:nelem)
        type(face_t),     allocatable   :: faces(:,:)                   !< Face storage    (1:nelem,1:nfaces)
        type(chimera_t)                 :: chimera                      !< Chimera interface data

        ! Initialization flags
        logical             :: geomInitialized = .false.                !< Status of geometry initialization
        logical             :: solInitialized  = .false.                !< Status of numerics initialization


    contains

        procedure           :: init_geom
        procedure           :: init_sol

        procedure, private  :: init_elems_geom
        procedure, private  :: init_elems_sol
        procedure, private  :: init_faces_geom
        procedure, private  :: init_faces_sol

        final               :: destructor

    end type mesh_t
    !************************************************************************************************************





contains






    !> Mesh geometry initialization procedure
    !!
    !!  Sets number of terms in coordinate expansion for the entire domain
    !!  and calls sub-initialization routines for individual element and face geometry
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]  nterms_c    Number of terms in the coordinate expansion
    !!  @param[in]  points_g    Rank-3 matrix of coordinate points defining a block mesh
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine init_geom(self,idomain,nterms_c,points_g)
        class(mesh_t),  intent(inout), target   :: self
        integer(ik),    intent(in)              :: idomain
        integer(ik),    intent(in)              :: nterms_c
        type(point_t),  intent(in)              :: points_g(:,:,:)
        type(element_t), pointer                :: temp(:)


        !
        ! Store number of terms in coordinate expansion and domain index
        !
        self%nterms_c = nterms_c
        self%idomain  = idomain


        !
        ! Call geometry initialization for elements and faces
        !
        call self%init_elems_geom(points_g)
        call self%init_faces_geom()


        !
        ! Confirm initialization
        !
        self%geomInitialized = .true.


    end subroutine init_geom
    !************************************************************************************************************









    !> Mesh numerics initialization procedure
    !!
    !!  Sets number of equations being solved, number of terms in the solution expansion and
    !!  calls sub-initialization routines for individual element and face numerics
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  neqns       Number of equations being solved in the current domain
    !!  @param[in]  nterms_s    Number of terms in the solution expansion
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine init_sol(self,neqns,nterms_s)
        class(mesh_t),  intent(inout)   :: self
        integer(ik),    intent(in)      :: neqns
        integer(ik),    intent(in)      :: nterms_s


        !
        ! Store number of equations and number of terms in solution expansion
        !
        self%neqns    = neqns
        self%nterms_s = nterms_s


        !
        ! Call numerics initialization for elements and faces
        !
        call self%init_elems_sol(neqns,nterms_s) 
        call self%init_faces_sol()               

        
        !
        ! Confirm initialization
        !
        self%solInitialized = .true.


    end subroutine init_sol
    !************************************************************************************************************







    !> Mesh - element initialization procedure
    !!
    !!  Computes the number of elements based on the element mapping selected and
    !!  calls the element initialization procedure on individual elements.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  points_g    Rank-3 matrix of coordinate points defining a block mesh
    !!
    !------------------------------------------------------------------------------------------------------------
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
                                        neqns,    nterms_s,  nnodes, nterms_c,  &
                                        npts_1d, mapping, idomain


        !
        ! Get number of points in each direction
        !
        npts_xi   = size(points_g,1)    ! Number of points in the xi-direction
        npts_eta  = size(points_g,2)    ! Number of points in the eta-direction
        npts_zeta = size(points_g,3)    ! Number of points in the zeta-direction


        !
        ! Compute number of 1d points for a single element
        !
        npts_1d = 0
        do while (npts_1d*npts_1d*npts_1d < self%nterms_c)
            npts_1d = npts_1d + 1       ! really just computing the cubed root of nterms_c, the number of terms in the coordinate expansion
        end do


        !
        ! Count number of elements in each direction and check mesh conforms to
        ! the agglomeration rule for higher-order elements
        !
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



        !
        ! Store number of elements in each direction along with total number of elements
        !
        self%nelem_xi   = nelem_xi
        self%nelem_eta  = nelem_eta
        self%nelem_zeta = nelem_zeta
        nelem           = nelem_xi * nelem_eta * nelem_zeta
        self%nelem      = nelem
        mapping         = (npts_1d - 1)     !> 1 - linear, 2 - quadratic, 3 - cubic, etc.


        !
        ! Allocate element storage
        !
        allocate(self%elems(nelem),       &
                 points_l(self%nterms_c), stat=ierr)
        if(ierr /= 0) stop "Memory allocation error: init_elements"



        idomain = self%idomain
        ielem = 1

        !
        ! Accumulate points for each element and call element initialization procedure
        !
        do izeta = 1,nelem_zeta
            do ieta = 1,nelem_eta
                do ixi = 1,nelem_xi

                    xi_start   = 1 + (ixi  -1)*(npts_1d-1)
                    eta_start  = 1 + (ieta -1)*(npts_1d-1)
                    zeta_start = 1 + (izeta-1)*(npts_1d-1)

                    !
                    ! For this element, collect the necessary points from the global points
                    ! array into a local points array for initializing an individual element
                    !
                    ipt = 1
                    do ipt_zeta = 1,npts_1d
                        do ipt_eta = 1,npts_1d
                            do ipt_xi = 1,npts_1d
                                points_l(ipt) = points_g(xi_start+(ipt_xi-1),eta_start+(ipt_eta-1),zeta_start+(ipt_zeta-1))
                                ipt = ipt + 1
                            end do
                        end do
                    end do


                    !
                    ! Element geometry initialization
                    !
                    call self%elems(ielem)%init_geom(mapping,points_l,idomain,ielem)
                    ielem = ielem + 1
                end do
            end do
        end do

    end subroutine init_elems_geom
    !**************************************************************************************************************









    !> Mesh - element solution data initialization
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  neqns       Number of equations in the domain equation set
    !!  @param[in]  nterms_s    Number of terms in the solution expansion
    !!
    !--------------------------------------------------------------------------------------------------------------
    subroutine init_elems_sol(self,neqns,nterms_s)
        class(mesh_t),  intent(inout)   :: self
        integer(ik),    intent(in)      :: neqns
        integer(ik),    intent(in)      :: nterms_s
        integer(ik) :: ielem


        !
        ! Store number of equations and number of terms in the solution expansion
        !
        self%neqns    = neqns
        self%nterms_s = nterms_s


        !
        ! Call the numerics initialization procedure for each element
        !
        do ielem = 1,self%nelem
            call self%elems(ielem)%init_sol(self%neqns,self%nterms_s)
        end do


    end subroutine init_elems_sol
    !***************************************************************************************************************












    !> Mesh - face initialization procedure
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------------------
    subroutine init_faces_geom(self)
        class(mesh_t), intent(inout)  :: self

        integer(ik)             :: ixi,ieta,izeta,iface,ftype,ineighbor,ielem,ierr
        logical                 :: boundary_face = .false.

        !
        ! Allocate face storage
        !
        allocate(self%faces(self%nelem,NFACES),stat=ierr)
        if (ierr /= 0) stop "Error: mesh%init -- face allocation error"


        !
        ! Loop through each element and call initialization for each face
        !
        ielem = 1
        do izeta = 1,self%nelem_zeta
            do ieta = 1,self%nelem_eta
                do ixi = 1,self%nelem_xi

                    !
                    ! For each face of the current element, call initialization procedure
                    !
                    do iface = 1,NFACES

                        !
                        ! Set ftype to designate interior and boundary faces
                        !
                        boundary_face = ( (ixi == 1                 .and. iface == XI_MIN)   .or. &
                                          (ixi == self%nelem_xi     .and. iface == XI_MAX)   .or. &
                                          (ieta == 1                .and. iface == ETA_MIN)  .or. &
                                          (ieta == self%nelem_eta   .and. iface == ETA_MAX)  .or. &
                                          (izeta == 1               .and. iface == ZETA_MIN) .or. &
                                          (izeta == self%nelem_zeta .and. iface == ZETA_MAX) )



                        if (boundary_face) then

                            !
                            ! Default ftype to ORPHAN face
                            !
                            ftype = ORPHAN      ! This should be processed later; either by a boundary condition(ftype=1), or a chimera boundary(ftype=2)
                            ineighbor = 0       ! No neighbor


                        else
                            ftype = INTERIOR    ! interior face

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


                        !
                        ! Call face initialization routine
                        !
                        call self%faces(ielem,iface)%init_geom(iface,ftype,self%elems(ielem),ineighbor)

                    end do

                    !
                    ! Increment element index
                    !
                    ielem = ielem + 1

                end do !ixi
            end do ! ieta
        end do ! izeta

    end subroutine init_faces_geom
    !**************************************************************************************************************






    !> Mesh - face initialization procedure
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------------------
    subroutine init_faces_sol(self)
        class(mesh_t), intent(inout)  :: self

        integer(ik) :: ielem, iface

        !
        ! Loop through elements
        !
        do ielem = 1,self%nelem

            !
            ! Loop through faces and call numerics initialization routine
            !
            do iface = 1,NFACES

                call self%faces(ielem,iface)%init_sol(self%elems(ielem))

            end do ! iface

        end do ! ielem

    end subroutine init_faces_sol
    !***************************************************************************************************************









    subroutine destructor(self)
        type(mesh_t), intent(inout) :: self

    
    end subroutine




end module type_mesh
