module type_chidgVector
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: ZERO, TWO
    use mod_chidg_mpi,              only: GROUP_MASTER, ChiDG_COMM, IRANK
    use type_mesh,                  only: mesh_t
    use type_chidgVector_send,      only: chidgVector_send_t
    use type_chidgVector_recv,      only: chidgVector_recv_t
    use type_blockvector
    use mpi_f08,                    only: MPI_AllReduce, MPI_Reduce, MPI_COMM, MPI_REAL8, MPI_SUM, MPI_STATUS_IGNORE, MPI_Recv, MPI_Request, MPI_STATUSES_IGNORE
    implicit none





    !>  High-level ChiDG vector container.
    !! 
    !!  Container stores a blockvector_t for each domain_t
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    type, public :: chidgVector_t

        type(blockvector_t),    allocatable :: dom(:)       !< Local block vector storage

        type(chidgVector_send_t)            :: send         !< Information on what to send to other processors
        type(chidgVector_recv_t)            :: recv         !< Storage to receive data from other processors

    contains

        generic,    public  :: init => initialize
        procedure,  private :: initialize

        procedure,  public  :: clear                            !< Zero the densevector data nested in the container
        generic,    public  :: norm => norm_local, norm_comm    !< Compute the L2 vector norm
        procedure,  public  :: norm_local                       !< Return the processor-local L2 vector norm of the chidgVector 
        procedure,  public  :: norm_comm                        !< Return the MPI group L2 vector norm of the chidgVector 
        procedure,  public  :: sumsqr                           !< Return the sum of the squared processor-local chidgVector entries 
        procedure,  public  :: dump

        procedure,  public  :: comm_send                        !< Execute non-blocking send of data to communicating processors
        procedure,  public  :: comm_recv                        !< Execute blocking receives of incomming data
        procedure,  public  :: comm_wait                        !< Execute a wait on outstanding non-blocking send data

    end type chidgVector_t
    !****************************************************************************************************










    !------------------------       OPERATORS       --------------------------------------

    public operator (*)
    interface operator(*)
        module procedure mult_real_chidgVector          ! real * chidgVector
        module procedure mult_chidgVector_real          ! chidgVector * real
    end interface


    public operator (/)
    interface operator (/)
        module procedure div_real_chidgVector           ! real / chidgVector
        module procedure div_chidgVector_real           ! chidgVector / real
    end interface


    public operator (-)
    interface operator (-)
        module procedure sub_chidgVector_chidgVector    ! chidgVector - chidgVector
    end interface

    public operator (+)
    interface operator (+)
        module procedure add_chidgVector_chidgVector    ! chidgVector + chidgVector
    end interface
















contains





    !>  Allocate and initialize chidgVector_t storage and data.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    Array of mesh_t instances used to initialize each blockvector_t subcomponent.
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine initialize(self,mesh)
        class(chidgVector_t),   intent(inout)   :: self
        type(mesh_t),           intent(inout)   :: mesh(:)

        integer(ik) :: ierr, ndomains, idom


        ! Allocate blockvector_t for each mesh
        ndomains = size(mesh)
        allocate(self%dom(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError



        ! Call initialization procedure for each blockvector_t
        do idom = 1,ndomains
            call self%dom(idom)%init(mesh(idom))
        end do


        ! Call initialization for determining what data to send and where
        call self%send%init(mesh)
        ! Call initialization for determining what data to receive and allocate storage for it
        call self%recv%init(mesh)


    end subroutine initialize
    !****************************************************************************************************













    !> Set all floating-point vector entries to zero.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(chidgVector_t),   intent(inout)   :: self

        integer :: idom


        ! Call clear procedure for each blockvector_t
        do idom = 1,size(self%dom)
            call self%dom(idom)%clear()
        end do

        ! Call clear on recv storage
        call self%recv%clear()

    end subroutine clear
    !*****************************************************************************************************










    !>  Compute the process-local L2-Norm of the vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @return res     L2-norm of the vector
    !!
    !-----------------------------------------------------------------------------------------------------
    function norm_local(self) result(res)
        class(chidgVector_t),   intent(in)   :: self

        real(rk)    :: res
        integer(ik) :: idom, ielem


        res = ZERO

        ! Loop through domain vectors and compute contribution to vector sum of the squared elements
        do idom = 1,size(self%dom)
            res = res + self%dom(idom)%sumsqr()
        end do ! idom


        ! Take the square root of the result
        res = sqrt(res)

    end function norm_local
    !*****************************************************************************************************











    !>  Compute the L2-Norm of the vector within the space of processors given by the MPI communicator
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   6/23/2016
    !!
    !!  @return res     L2-norm of the vector
    !!
    !-----------------------------------------------------------------------------------------------------
    function norm_comm(self,comm) result(norm)
        class(chidgVector_t),   intent(in)  :: self
        type(mpi_comm),         intent(in)  :: comm

        real(rk)    :: sumsqr, norm
        integer     :: ierr


        norm = ZERO

        ! Compute sum of the squared elements of the processor-local vector
        sumsqr = self%sumsqr()


        ! Reduce sumsqr across all procs, distribute result back to all
        call MPI_AllReduce(sumsqr,norm,1,MPI_REAL8,MPI_SUM,comm,ierr)


        norm = sqrt(norm)

    end function norm_comm
    !*****************************************************************************************************














    !< Return the sum of the squared chidgVector entries 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/23/2016
    !!
    !!  @return res     sum of the squared chidgVector entries
    !!
    !-----------------------------------------------------------------------------------------------------
    function sumsqr(self) result(res)
        class(chidgVector_t),   intent(in)   :: self

        real(rk)    :: res
        integer(ik) :: idom, ielem


        res = ZERO

        ! Loop through domain vectors and compute contribution to vector sum of the squared elements
        do idom = 1,size(self%dom)
            res = res + self%dom(idom)%sumsqr()
        end do ! idom


    end function sumsqr
    !*****************************************************************************************************















    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------------
    subroutine comm_send(self)
        class(chidgVector_t),   intent(inout)   :: self

        integer(ik)         :: icomm, iproc_send, idom_send, idom, ielem_send, ielem, ierr, data_size, isend
        type(mpi_request)   :: isend_handle


        !
        ! Loop through comms to send
        !
        isend = 1
        do icomm = 1,size(self%send%comm)

            ! Get processor rank we are sending to
            iproc_send = self%send%comm(icomm)%proc

            ! Loop through domains to send
            do idom_send = 1,self%send%comm(icomm)%dom_send%size()

                ! Get domain index to be send
                idom = self%send%comm(icomm)%dom_send%at(idom_send)

                do ielem_send = 1,self%send%comm(icomm)%elems_send(idom_send)%size()
                    ielem = self%send%comm(icomm)%elems_send(idom_send)%at(ielem_send)

                    ! Post non-blocking send message for the vector data
                    data_size = size(self%dom(idom)%vecs(ielem)%vec)
                    call MPI_ISend(self%dom(idom)%vecs(ielem)%vec, data_size, MPI_REAL8, iproc_send, 0, ChiDG_COMM, isend_handle, ierr)

                    ! Add non-blocking send handle to list of things to wait on
                    self%send%isend_handles(isend) = isend_handle

                    ! Increment send counter
                    isend = isend + 1
            
                end do !ielem_send
            end do !idom_send

        end do ! icomm



    end subroutine comm_send
    !*****************************************************************************************************












    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------------
    subroutine comm_recv(self)
        class(chidgVector_t),   intent(inout)   :: self

        integer(ik) :: icomm, idom_recv, ielem_recv, proc_recv, data_size, ierr
        integer(ik) :: idom_store, ielem_store


        !
        ! Receive data from each communicating processor
        !
        do icomm = 1,size(self%recv%comm)

            proc_recv = self%recv%comm(icomm)%proc
            
            do idom_recv = 1,size(self%recv%comm(icomm)%dom)

                ! Get storage index of incoming domain
                idom_store = self%recv%comm(icomm)%dom_store%at(idom_recv)


                do ielem_recv = 1,size(self%recv%comm(icomm)%dom(idom_recv)%vecs)

                    ! Get storage index of incoming element
                    ielem_store = self%recv%comm(icomm)%elem_store(idom_recv)%at(ielem_recv)

                    data_size = size(self%recv%comm(icomm)%dom(idom_store)%vecs(ielem_store)%vec)
                    call MPI_Recv(self%recv%comm(icomm)%dom(idom_store)%vecs(ielem_store)%vec, data_size, MPI_REAL8, proc_recv, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)

                    !data_size = size(self%recv%comm(icomm)%dom(idom_recv)%vecs(ielem_recv)%vec)
                    !call MPI_Recv(self%recv%comm(icomm)%dom(idom_recv)%vecs(ielem_recv)%vec, data_size, MPI_REAL8, proc_recv, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)

                end do ! ielem_recv
            end do ! idom_recv

        end do ! icomm



    end subroutine comm_recv
    !*****************************************************************************************************
















    !>  Wait for all non-blocking sends to complete.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine comm_wait(self)
        class(chidgVector_t),   intent(in)  :: self

        integer(ik) :: nwait, ierr


        nwait = size(self%send%isend_handles)
        call MPI_Waitall(nwait, self%send%isend_handles, MPI_STATUSES_IGNORE, ierr)

    end subroutine comm_wait
    !****************************************************************************************************













    !> Dump contents of the vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------------
    subroutine dump(self)
        class(chidgVector_t),   intent(in)   :: self

        integer(ik) :: idom

        ! Loop through domain vectors and compute contribution to vecotr L2-Norm
        do idom = 1,size(self%dom)
            call self%dom(idom)%dump()
        end do ! idom

    end subroutine dump
    !*****************************************************************************************************













    !---------------------------------------------------------------------------------------------
    !
    !
    !                              Operator Implementations
    !
    !---------------------------------------------------------------------------------------------



    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    function mult_real_chidgVector(left,right) result(res)
        real(rk),               intent(in)  :: left
        type(chidgVector_t),    intent(in)  :: right

        type(chidgVector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(right%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(right%dom)
            res%dom(idom) = left * right%dom(idom)
        end do

        res%send = right%send
        res%recv = right%recv

    end function mult_real_chidgVector
    !************************************************************************




    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------
    function mult_chidgVector_real(left,right) result(res)
        type(chidgVector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right

        type(chidgVector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(left%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(left%dom)
            res%dom(idom) = left%dom(idom) * right
        end do

        res%send = left%send
        res%recv = left%recv

    end function mult_chidgVector_real
    !************************************************************************





    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------
    function div_real_chidgVector(left,right) result(res)
        real(rk),               intent(in)  :: left
        type(chidgVector_t),    intent(in)  :: right

        type(chidgVector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(right%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(right%dom)
            res%dom(idom) = left / right%dom(idom)
        end do


        res%send = right%send
        res%recv = right%recv

    end function div_real_chidgVector
    !************************************************************************




    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------
    function div_chidgVector_real(left,right) result(res)
        type(chidgVector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right

        type(chidgVector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(left%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(left%dom)
            res%dom(idom) = left%dom(idom) / right
        end do

        res%send = left%send
        res%recv = left%recv

    end function div_chidgVector_real
    !*************************************************************************





    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------
    function add_chidgVector_chidgVector(left,right) result(res)
        type(chidgVector_t),    intent(in)  :: left
        type(chidgVector_t),    intent(in)  :: right

        type(chidgVector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(right%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(left%dom)
            res%dom(idom) = left%dom(idom) + right%dom(idom)
        end do


        res%send = right%send
        res%recv = right%recv

    end function add_chidgVector_chidgVector
    !*************************************************************************





    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------
    function sub_chidgVector_chidgVector(left,right) result(res)
        type(chidgVector_t),    intent(in)  :: left
        type(chidgVector_t),    intent(in)  :: right

        type(chidgVector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(right%dom)

        allocate(res%dom(ndom))


        do idom = 1,size(left%dom)
            res%dom(idom) = left%dom(idom) - right%dom(idom)
        end do

        res%send = right%send
        res%recv = right%recv

    end function sub_chidgVector_chidgVector
    !*************************************************************************
























end module type_chidgVector
