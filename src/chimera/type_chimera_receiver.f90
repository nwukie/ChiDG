module type_chimera_receiver
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: NO_PROC, NO_ID
    use type_ivector,       only: ivector_t
    use type_mvector,       only: mvector_t
    use type_pvector,       only: pvector_t
    use type_rvector,       only: rvector_t
    use type_chimera_donor, only: chimera_donor_t, chimera_donor
    implicit none




    !> Chimera receiver data container. One instance per Chimera face receiving data.
    !! Holds donor domain indices and donor element indices corresponding to their
    !! location in the donor domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------------------------
    type, public :: chimera_receiver_t

        integer(ik)                     :: idomain_g    ! ChiDG-global domain index of receiver
        integer(ik)                     :: idomain_l    ! Proc-local domain index of receiver
        integer(ik)                     :: ielement_g   ! Domain-global element index of receiver
        integer(ik)                     :: ielement_l   ! Proc-local element index of receiver
        integer(ik)                     :: iface        ! Face index of receiver


!        ! Access via data%donor_domain%at(idonor)
!        type(ivector_t)                 :: donor_eqn_ID             ! Equation set identifier
!        type(ivector_t)                 :: donor_neqns              ! Number of equations in donor element
!        type(ivector_t)                 :: donor_nterms_s           ! Number of terms in donor expansion
!        type(ivector_t)                 :: donor_proc               ! Vector of processor ranks
!        type(ivector_t)                 :: donor_domain_g           ! Vector of domain indices
!        type(ivector_t)                 :: donor_domain_l           ! Vector of domain indices
!        type(ivector_t)                 :: donor_element_g          ! Vector of element indices for the location in the corresponding domain
!        type(ivector_t)                 :: donor_element_l          ! Vector of element indices for the location in the corresponding domain
!        type(mvector_t)                 :: donor_interpolator       ! Vector of matrices defining the Chimera interpolation
!        type(mvector_t)                 :: donor_interpolator_grad1
!        type(mvector_t)                 :: donor_interpolator_grad2
!        type(mvector_t)                 :: donor_interpolator_grad3
!
!        ! For donor elements being comminucated from off-processor, their location in the recv container for accessing components
!        type(ivector_t)                 :: donor_recv_comm          ! Indices of comm container
!        type(ivector_t)                 :: donor_recv_domain        ! Domain index within comm
!        type(ivector_t)                 :: donor_recv_element       ! Element index within domain
!
!        ! The access for this component is slightly different than the above components
!        ! Access via data%donor_gq_indices(idonor)%data()
!        type(ivector_t), allocatable    :: donor_gq_indices(:)      ! Array of integer vectors defining the GQ node indices associated with a given donor
!        type(pvector_t), allocatable    :: donor_coords(:)          ! Array of points defining the local coordinates of the GQ nodes
!        type(mvector_t), allocatable    :: donor_metrics(:)         ! For each donor, matrices of metric terms for each donor GQ node
!        type(rvector_t), allocatable    :: donor_jinv(:)            ! For each donor, inverse element jacobian term for each donor GQ node


        !
        ! NEW
        !
        type(chimera_donor_t), allocatable :: donor(:)


        real(rk),   allocatable :: jinv(:)
        real(rk),   allocatable :: det_jacobian_grid(:)
        real(rk),   allocatable :: inv_jacobian_grid(:)

    contains

        procedure   :: add_donor
        procedure   :: new_donor
        procedure   :: find_donor
        procedure   :: ndonors

        procedure   :: clear

    end type chimera_receiver_t
    !******************************************************************************************




contains





    !>  Add chimera donor to the receiver face object.
    !!
    !!  If donor matching the incoming donor already exists, return its ID instead of
    !!  creating a new instance.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !----------------------------------------------------------------------------------------
    function add_donor(self, idomain_g, idomain_l, ielement_g, ielement_l, iproc) result(donor_ID)
        class(chimera_receiver_t),  intent(inout)   :: self
        integer(ik),                intent(in)      :: idomain_g
        integer(ik),                intent(in)      :: idomain_l
        integer(ik),                intent(in)      :: ielement_g
        integer(ik),                intent(in)      :: ielement_l
        integer(ik),                intent(in)      :: iproc

        integer(ik) :: donor_ID

        ! Check if receiver matching the incoming face already exists.
        ! If not, call new and construct new object.
        donor_ID = self%find_donor(idomain_g,ielement_g)
        if ( donor_ID == NO_ID ) then
            donor_ID = self%new_donor()
            self%donor(donor_ID) = chimera_donor(idomain_g, idomain_l, ielement_g, ielement_l, iproc)
        end if

    end function add_donor
    !****************************************************************************************






    !>  Extend allocation of donor instances, return identifier for new instance.
    !!
    !!  Returns donor_ID for chimera_donor_t in domain%chimera%recv(Chi_ID)%donor(donor_ID)
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !-----------------------------------------------------------------------------------------
    function new_donor(self) result(donor_ID)
        class(chimera_receiver_t),  intent(inout)   :: self

        type(chimera_donor_t),  allocatable :: temp(:)
        integer(ik)                         :: ierr, donor_ID



        ! Resize array storage
        allocate(temp(self%ndonors() + 1), stat=ierr)


        ! Copy previously initialized instances to new array.
        if (self%ndonors() > 0) then
            temp(1:size(self%donor)) = self%donor(1:size(self%donor))
        end if


        ! Move resized temp allocation back to parent. 
        call move_alloc(temp,self%donor)


        ! Set domain identifier of newly allocated domain that will be returned
        donor_ID = self%ndonors()


    end function new_donor
    !*****************************************************************************************






    !>  Given a global element location, return the index of a chimera_donor instance 
    !!  corresponding to the element, if one exists in the chimera collection.
    !!
    !!  If no chimera_donor instance can be found that matches the incoming description
    !!  the return NO_ID.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !------------------------------------------------------------------------------------------
    function find_donor(self,idomain_g,ielement_g) result(donor_ID)
        class(chimera_receiver_t),  intent(in)  :: self
        integer(ik),                intent(in)  :: idomain_g
        integer(ik),                intent(in)  :: ielement_g

        integer(ik) :: donor_ID, idonor

        donor_ID = NO_ID
        do idonor = 1,self%ndonors()
            if ( (self%donor(idonor)%idomain_g  == idomain_g ) .and. &
                 (self%donor(idonor)%ielement_g == ielement_g) ) then
                 donor_ID = idonor
                 exit
            end if
        end do !idonor

    end function find_donor
    !******************************************************************************************











!    !>  Return the number of donor elements providing the receiver with data.
!    !!
!    !!  @author Nathan A. Wukie (AFRL)
!    !!  @date   7/13/2016
!    !!
!    !!
!    !!
!    !-------------------------------------------------------------------------------------------
!    function ndonors(self) result(res)
!        class(chimera_receiver_t), intent(in)  :: self
!
!        integer(ik) :: res
!
!        res = self%donor_neqns%size()
!
!    end function ndonors
!    !*******************************************************************************************


    !>  Return the number of donor elements providing the receiver with data.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/13/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function ndonors(self) result(n)
        class(chimera_receiver_t), intent(in)  :: self

        integer(ik) :: n

        if (allocated(self%donor)) then
            n = size(self%donor)
        else
            n = 0
        end if

    end function ndonors
    !*******************************************************************************************





    !>  Clear the data for the chimera receiver face.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/3/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(chimera_receiver_t), intent(inout)   :: self

        self%idomain_g    = 0
        self%idomain_l    = 0
        self%ielement_g   = 0
        self%ielement_l   = 0
        self%iface        = 0


!        call self%donor_eqn_ID%clear()
!        call self%donor_neqns%clear()
!        call self%donor_nterms_s%clear()
!        call self%donor_proc%clear()               
!        call self%donor_domain_g%clear()           
!        call self%donor_domain_l%clear()           
!        call self%donor_element_g%clear()         
!        call self%donor_element_l%clear()         
!        call self%donor_interpolator%clear()       
!        call self%donor_interpolator_grad1%clear()
!        call self%donor_interpolator_grad2%clear()
!        call self%donor_interpolator_grad3%clear()
!
!        call self%donor_recv_comm%clear()
!        call self%donor_recv_domain%clear()
!        call self%donor_recv_element%clear()
!
!        if (allocated(self%donor_gq_indices))   deallocate(self%donor_gq_indices) 
!        if (allocated(self%donor_coords))       deallocate(self%donor_coords) 
!        if (allocated(self%donor_metrics))      deallocate(self%donor_metrics) 
!        if (allocated(self%donor_jinv))         deallocate(self%donor_jinv) 
        if (allocated(self%donor))             deallocate(self%donor)
        if (allocated(self%jinv))              deallocate(self%jinv)
        if (allocated(self%det_jacobian_grid)) deallocate(self%det_jacobian_grid)
        if (allocated(self%inv_jacobian_grid)) deallocate(self%inv_jacobian_grid)

    end subroutine clear
    !******************************************************************************************




end module type_chimera_receiver
