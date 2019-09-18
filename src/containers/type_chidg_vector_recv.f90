module type_chidg_vector_recv
#include <messenger.h>
    use mod_kinds,                   only: ik
    use mod_constants,              only: NO_ID
    use type_ivector,                only: ivector_t
    use type_mesh,                   only: mesh_t
    use type_chidg_vector_recv_comm, only: chidg_vector_recv_comm_t
    implicit none



    !>  Container for receiving vector data from other chidg_vectors. Holds a container for
    !!  each process that is sending data here.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    type, public :: chidg_vector_recv_t

        type(chidg_vector_recv_comm_t),   allocatable :: comm(:)


    contains

        procedure, public   :: init
        procedure, public   :: clear
        procedure, public   :: ncomm
        procedure, public   :: restrict
        procedure, public   :: prolong

    end type chidg_vector_recv_t
    !*****************************************************************************************





contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self,mesh)
        class(chidg_vector_recv_t), intent(inout)   :: self
        type(mesh_t),               intent(inout)   :: mesh

        integer(ik)                 :: idom, iproc, icomm, ncomm, ierr, pelem_ID
        integer(ik),    allocatable :: comm_procs_array(:)


        !
        ! Deallocate storage if necessary in case this is being called as a 
        ! reinitialization routine.
        !
        if (allocated(self%comm)) deallocate(self%comm)


        !
        ! Get processor ranks that we are receiving from: mesh
        !
        comm_procs_array = mesh%get_recv_procs()


        !
        ! Allocate recv communication for each processor sending data here.
        !
        ncomm = size(comm_procs_array)
        if (allocated(self%comm)) deallocate(self%comm)
        allocate(self%comm(ncomm), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Call initialization for each communicating process
        !
        do icomm = 1,ncomm
            iproc = comm_procs_array(icomm)
            call self%comm(icomm)%init(mesh,iproc,icomm)
        end do


    end subroutine init
    !*******************************************************************************************







    !>  Clear the received data, but do NOT deallocated.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/7/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(chidg_vector_recv_t),  intent(inout)   :: self
        
        integer(ik) :: icomm

    
        do icomm = 1,size(self%comm)
            call self%comm(icomm)%clear()
        end do !icomm


    end subroutine clear
    !********************************************************************************************




    !>
    !!
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function ncomm(self) result(ncomm_)
        class(chidg_vector_recv_t), intent(in)  :: self

        integer(ik) :: ncomm_

        if (allocated(self%comm)) then
            ncomm_ = size(self%comm)
        else 
            ncomm_ = 0
        end if

    end function ncomm
    !*******************************************************************************************



    !>
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/21/2017
    !!
    !------------------------------------------------------------------------------------------
    function restrict(self, nterms_r) result(restricted)
        class(chidg_vector_recv_t), intent(in)  :: self
        integer(ik),                intent(in)  :: nterms_r

        type(chidg_vector_recv_t)   :: restricted
        integer(ik)                 :: ierr, icomm

        allocate(restricted%comm(self%ncomm()), stat=ierr)
        if (ierr /= 0) call AllocationError

        do icomm = 1,self%ncomm()
            restricted%comm(icomm) = self%comm(icomm)%restrict(nterms_r)
        end do !icomm

    end function restrict
    !******************************************************************************************




    !>
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/21/2017
    !!
    !------------------------------------------------------------------------------------------
    function prolong(self, nterms_p) result(prolonged)
        class(chidg_vector_recv_t), intent(in)  :: self
        integer(ik),                intent(in)  :: nterms_p

        type(chidg_vector_recv_t)   :: prolonged
        integer(ik)                 :: ierr, icomm

        allocate(prolonged%comm(self%ncomm()), stat=ierr)
        if (ierr /= 0) call AllocationError

        do icomm = 1,self%ncomm()
            prolonged%comm(icomm) = self%comm(icomm)%prolong(nterms_p)
        end do !icomm

    end function prolong
    !******************************************************************************************










end module type_chidg_vector_recv
