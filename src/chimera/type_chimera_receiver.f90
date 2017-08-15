module type_chimera_receiver
#include <messenger.h>
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

        !
        ! Array of donors
        !
        type(chimera_donor_t), allocatable :: donor(:)


        !
        ! Data assembled from all donors to define the complete node set
        !
        real(rk),   allocatable :: ale_g(:)
        real(rk),   allocatable :: ale_g_grad1(:)
        real(rk),   allocatable :: ale_g_grad2(:)
        real(rk),   allocatable :: ale_g_grad3(:)
        real(rk),   allocatable :: ale_Dinv(:,:,:)
        real(rk),   allocatable :: interp_coords_vel(:,:)

    contains


        procedure   :: add_donor
        procedure   :: new_donor
        procedure   :: find_donor
        procedure   :: ndonors

        procedure   :: clear

    end type chimera_receiver_t
    !******************************************************************************************



    interface chimera_receiver
        module procedure chimera_receiver
    end interface chimera_receiver




contains



    !>  Constructor for chimera_receiver_t.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/27/2017
    !!
    !---------------------------------------------------------------------------------------
    function chimera_receiver(idomain_g, idomain_l, ielement_g, ielement_l, iface, nnodes) result(instance)
        integer(ik),    intent(in)  :: idomain_g
        integer(ik),    intent(in)  :: idomain_l
        integer(ik),    intent(in)  :: ielement_g
        integer(ik),    intent(in)  :: ielement_l
        integer(ik),    intent(in)  :: iface
        integer(ik),    intent(in)  :: nnodes

        integer(ik)                 :: ierr
        type(chimera_receiver_t)    :: instance


        instance%idomain_g  = idomain_g
        instance%idomain_l  = idomain_l
        instance%ielement_g = ielement_g
        instance%ielement_l = ielement_l
        instance%iface      = iface


        allocate(instance%ale_g(nnodes),        &
                 instance%ale_g_grad1(nnodes),  &
                 instance%ale_g_grad2(nnodes),  &
                 instance%ale_g_grad3(nnodes),  &
                 instance%ale_Dinv(3,3,nnodes), &
                 instance%interp_coords_vel(nnodes,3), stat=ierr)
        if (ierr /= 0) call AllocationError


    end function chimera_receiver
    !***************************************************************************************







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

        if (allocated(self%donor))             deallocate(self%donor)
        if (allocated(self%ale_g))             deallocate(self%ale_g)
        if (allocated(self%ale_g_grad1))       deallocate(self%ale_g_grad1)
        if (allocated(self%ale_g_grad2))       deallocate(self%ale_g_grad2)
        if (allocated(self%ale_g_grad3))       deallocate(self%ale_g_grad3)
        if (allocated(self%ale_Dinv))          deallocate(self%ale_Dinv)
        if (allocated(self%interp_coords_vel)) deallocate(self%interp_coords_vel)

    end subroutine clear
    !******************************************************************************************




end module type_chimera_receiver
