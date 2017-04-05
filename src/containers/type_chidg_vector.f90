module type_chidg_vector
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: ZERO, TWO
    use mod_chidg_mpi,              only: GROUP_MASTER, ChiDG_COMM, IRANK
    use type_mesh,                  only: mesh_t
    use type_mesh_new,              only: mesh_new_t
    use type_function,              only: function_t
    use type_chidg_vector_send,     only: chidg_vector_send_t
    use type_chidg_vector_recv,     only: chidg_vector_recv_t
    use type_blockvector
    use mpi_f08,                    only: MPI_AllReduce, MPI_Reduce, MPI_COMM, MPI_REAL8,    &
                                          MPI_SUM, MPI_STATUS_IGNORE, MPI_Recv, MPI_Request, &
                                          MPI_STATUSES_IGNORE, MPI_INTEGER4
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
    !------------------------------------------------------------------------------------------
    type, public :: chidg_vector_t

        type(blockvector_t),    allocatable :: dom(:)       !< Local block vector storage

        type(chidg_vector_send_t)           :: send         !< What to send to other processors
        type(chidg_vector_recv_t)           :: recv         !< Receive data from other processors

        integer(ik),            private     :: ntime_       !< No. of time instances stored

    contains

        generic,    public  :: init => initialize
        procedure,  private :: initialize

        procedure,  public  :: project                          ! Project function to basis
        procedure,  public  :: clear                            ! Zero the densevector data
        generic,    public  :: norm => norm_local, norm_comm    ! Compute L2 vector norm
        procedure,  public  :: norm_local                       ! proc-local L2 vector norm
        procedure,  public  :: norm_comm                        ! MPI group L2 vector norm
        generic,    public  :: norm_fields => norm_fields_comm  ! L2 norm of independent fields
        procedure,  public  :: norm_fields_comm                 ! MPI group L2 field norms
        procedure,  public  :: sumsqr                           ! Sum squared proc-local entries 
        procedure,  public  :: sumsqr_fields
        procedure,  public  :: dump

        procedure,  public  :: comm_send                        ! Nonblocking send to comm procs
        procedure,  public  :: comm_recv                        ! Blocking recv incomming data
        procedure,  public  :: comm_wait                        ! Wait to finish send data

        procedure,  public  :: release                          ! Release allocated resources
        procedure,  public  :: get_ntime                        ! Return ntime associated with
                                                                ! densevctors
        procedure,  public  :: set_ntime                        ! Set ntime in the associated
                                                                ! densevectors
!        generic :: assignment(=) => 

    end type chidg_vector_t
    !*****************************************************************************************










    !------------------------       OPERATORS       --------------------------------------

    public operator (*)
    interface operator(*)
        module procedure mult_real_chidg_vector          ! real * chidg_vector
        module procedure mult_chidg_vector_real          ! chidg_vector * real
    end interface


    public operator (/)
    interface operator (/)
        module procedure div_real_chidg_vector           ! real / chidg_vector
        module procedure div_chidg_vector_real           ! chidg_vector / real
    end interface


    public operator (-)
    interface operator (-)
        module procedure sub_chidg_vector_chidg_vector    ! chidg_vector - chidg_vector
    end interface

    public operator (+)
    interface operator (+)
        module procedure add_chidg_vector_chidg_vector    ! chidg_vector + chidg_vector
    end interface
















contains





    !>  Allocate and initialize chidg_vector_t storage and data.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    Array of mesh_t instances used to initialize each 
    !!                      blockvector_t subcomponent.
    !!
    !------------------------------------------------------------------------------------------
    subroutine initialize(self,mesh,ntime)
        class(chidg_vector_t),   intent(inout)   :: self
        !type(mesh_t),            intent(inout)   :: mesh(:)
        type(mesh_new_t),        intent(inout)   :: mesh
        integer(ik),             intent(in)      :: ntime

        integer(ik) :: ierr, ndomains, idom


        !
        ! Set ntime_ for the chidg_vector
        !
        self%ntime_ = ntime

        !
        ! Deallocate storage if necessary in case this is being called as a 
        ! reinitialization routine.
        !
        if (allocated(self%dom)) deallocate(self%dom)


        ! Allocate blockvector_t for each mesh
        ndomains = mesh%ndomains()
        allocate(self%dom(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError



        ! Call initialization procedure for each blockvector_t
        do idom = 1,ndomains
            call self%dom(idom)%init(mesh%domain(idom))
        end do


        ! Call initialization for determining what data to send and where
        call self%send%init(mesh)
        ! Call initialization for determining what data to receive and allocate storage for it
        call self%recv%init(mesh)


        ! Wait on outstanding mpi_reqests initiated during the send%init(mesh) call
        call self%send%init_wait()


    end subroutine initialize
    !******************************************************************************************







    !>  Project a function onto the global solution basis.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/25/2016
    !!
    !!  @author Mayank Sharma
    !!  @date   11/12/2016
    !!
    !!  TODO: Should itime be an input parameter here?
    !!
    !------------------------------------------------------------------------------------------
    subroutine project(self,mesh,fcn,ivar)
        class(chidg_vector_t),  intent(inout)   :: self
        !type(mesh_t),           intent(in)      :: mesh(:)
        type(mesh_new_t),       intent(in)      :: mesh
        class(function_t),      intent(inout)   :: fcn
        integer(ik),            intent(in)      :: ivar

        integer(ik)                 :: idom, ielem, ierr
        integer(ik)                 :: itime
        real(rk),       allocatable :: fmodes(:)
        character(:),   allocatable :: user_msg


        !
        ! Loop through elements in mesh and call function projection
        !
        do idom = 1,mesh%ndomains()

            ! Check that variable index 'ivar' is valid
            user_msg = 'project: variable index ivar exceeds the number of equations.'
            if (ivar > mesh%domain(idom)%neqns ) call chidg_signal(FATAL,user_msg)

            do ielem = 1,mesh%domain(idom)%nelem
                do itime = 1,mesh%domain(idom)%ntime
                    !
                    ! Call function projection
                    !
                    fmodes = mesh%domain(idom)%elems(ielem)%project(fcn)


                    !
                    ! Store the projected modes to the solution expansion
                    !
                    call self%dom(idom)%vecs(ielem)%setvar(ivar,itime,fmodes)

                end do ! itime
            end do ! ielem
        end do ! idomain


    end subroutine project
    !******************************************************************************************






    !>  Set all floating-point vector entries to zero.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(chidg_vector_t),   intent(inout)   :: self

        integer :: idom


        ! Call clear procedure for each blockvector_t
        do idom = 1,size(self%dom)
            call self%dom(idom)%clear()
        end do

        ! Call clear on recv storage
        call self%recv%clear()

    end subroutine clear
    !******************************************************************************************








    !>  Compute the process-local L2-Norm of the vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @return res     L2-norm of the vector
    !!
    !------------------------------------------------------------------------------------------
    function norm_local(self) result(res)
        class(chidg_vector_t),   intent(in)   :: self

        real(rk)    :: res
        integer(ik) :: idom, ielem


        ! Loop through domain vectors and compute contribution to vector sum of the squared elements
        res = self%sumsqr()

        ! Take the square root of the result
        res = sqrt(res)

    end function norm_local
    !******************************************************************************************











    !>  Compute the L2-Norm of the vector within the space of processors given by the MPI 
    !!  communicator.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   6/23/2016
    !!
    !!  @return res     L2-norm of the vector
    !!
    !------------------------------------------------------------------------------------------
    function norm_comm(self,comm) result(norm)
        class(chidg_vector_t),   intent(in)  :: self
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
    !******************************************************************************************













    !>  Compute the L2-Norm of the vector within the space of processors given by the MPI 
    !!  communicator.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   6/23/2016
    !!
    !!  @return res     L2-norm of the vector
    !!
    !------------------------------------------------------------------------------------------
    function norm_fields_comm(self,comm) result(norm)
        class(chidg_vector_t),   intent(in)  :: self
        type(mpi_comm),         intent(in)  :: comm

        real(rk), allocatable   :: sumsqr(:), norm(:)
        integer     :: ierr


        ! Compute sum of the squared elements of the processor-local vector
        sumsqr = self%sumsqr_fields()


        ! Alloate norm
        norm = sumsqr
        norm = ZERO


        ! Reduce sumsqr across all procs, distribute result back to all
        call MPI_AllReduce(sumsqr,norm,size(sumsqr),MPI_REAL8,MPI_SUM,comm,ierr)


        norm = sqrt(norm)

    end function norm_fields_comm
    !******************************************************************************************
















    !< Return the sum of the squared chidg_vector entries 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/23/2016
    !!
    !!  @return res     sum of the squared chidg_vector entries
    !!
    !------------------------------------------------------------------------------------------
    function sumsqr(self) result(res)
        class(chidg_vector_t),   intent(in)   :: self

        real(rk)    :: res
        integer(ik) :: idom, ielem


        res = ZERO

        ! Loop through domain vectors and compute contribution to vector sum of the squared elements
        do idom = 1,size(self%dom)
            res = res + self%dom(idom)%sumsqr()
        end do ! idom


    end function sumsqr
    !******************************************************************************************











    !< Return the sum of the squared chidg_vector entries for each field independently.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/17/2017
    !!
    !!  @return res     sum of the squared chidg_vector entries
    !!
    !------------------------------------------------------------------------------------------
    function sumsqr_fields(self) result(res)
        class(chidg_vector_t),   intent(in)   :: self

        real(rk),   allocatable :: res(:)
        integer(ik) :: idom, ielem


        ! Allocate size of res based on assumption of same equation set across domains.
        res = self%dom(1)%sumsqr_fields()
        res = ZERO


        ! Loop through domain vectors and compute contribution to vector sum of the squared elements
        do idom = 1,size(self%dom)
            res = res + self%dom(idom)%sumsqr_fields()
        end do ! idom


    end function sumsqr_fields
    !******************************************************************************************


















    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine comm_send(self)
        class(chidg_vector_t),   intent(inout)   :: self

        integer(ik)         :: icomm, iproc_send, idom_send, idom, ielem_send, &
                               ielem, ierr, data_size, isend
        type(mpi_request)   :: isend_handle


        !
        ! Loop through comms to send
        !
        isend = 1
        do icomm = 1,size(self%send%comm)

            ! Get processor rank we are sending to
            iproc_send = self%send%comm(icomm)%proc

            ! Loop through domains/elements to send
            do idom_send = 1,self%send%comm(icomm)%dom_send%size()
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
    !*****************************************************************************************












    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine comm_recv(self)
        class(chidg_vector_t),   intent(inout)   :: self

        integer(ik) :: icomm, idom_recv, ielem_recv, proc_recv, data_size, &
                       ierr, dparent_g, eparent_g

        real(rk), allocatable   :: test(:)

        !
        ! Receive data from each communicating processor
        !
        do icomm = 1,size(self%recv%comm)

            ! Get process we are receiving from
            proc_recv = self%recv%comm(icomm)%proc
            
            ! Recv each element chunk
            do idom_recv = 1,size(self%recv%comm(icomm)%dom)
                do ielem_recv = 1,size(self%recv%comm(icomm)%dom(idom_recv)%vecs)

                    data_size = size(self%recv%comm(icomm)%dom(idom_recv)%vecs(ielem_recv)%vec)
                    call MPI_Recv(self%recv%comm(icomm)%dom(idom_recv)%vecs(ielem_recv)%vec, data_size, MPI_REAL8, proc_recv, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)

                end do ! ielem_recv
            end do ! idom_recv


        end do ! icomm



    end subroutine comm_recv
    !*****************************************************************************************
















    !>  Wait for all non-blocking sends to complete.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine comm_wait(self)
        class(chidg_vector_t),   intent(in)  :: self

        integer(ik) :: nwait, ierr


        nwait = size(self%send%isend_handles)
        call MPI_Waitall(nwait, self%send%isend_handles, MPI_STATUSES_IGNORE, ierr)

    end subroutine comm_wait
    !*****************************************************************************************













    !> Dump contents of the vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine dump(self)
        class(chidg_vector_t),   intent(in)   :: self

        integer(ik) :: idom

        ! Loop through domain vectors and compute contribution to vecotr L2-Norm
        do idom = 1,size(self%dom)
            call self%dom(idom)%dump()
        end do ! idom

    end subroutine dump
    !****************************************************************************************







    !>  Release allocated resources.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/3/2017
    !!
    !----------------------------------------------------------------------------------------
    subroutine release(self)
        class(chidg_vector_t),  intent(inout)   :: self

        if (allocated(self%dom)) deallocate(self%dom)

    end subroutine release
    !****************************************************************************************







    !>  Return ntime
    !!
    !!  @author Mayank Sharma
    !!  @date   3/9/2017
    !!
    !----------------------------------------------------------------------------------------
    function get_ntime(self) result(ntime_out)
        class(chidg_vector_t),  intent(inout)   :: self

        integer(ik)     :: ntime_out

        !
        ! Get ntime 
        !
        ntime_out = self%ntime_


    end function get_ntime
    !****************************************************************************************







    !>  Set ntime in the densevectors
    !!
    !!  @author Mayank Sharma
    !!  @date   3/9/2017
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_ntime(self,ntime)
        class(chidg_vector_t),  intent(inout)   :: self
        integer(ik),            intent(in)      :: ntime

        integer(ik)     :: idom, ielem


        ! 
        ! Set ntime
        !
        do idom = 1,size(self%dom)
            do ielem = 1,size(self%dom(idom)%vecs)
                
                call self%dom(idom)%vecs(ielem)%set_ntime(ntime)

            end do
        end do

    end subroutine set_ntime
    !****************************************************************************************







    !-----------------------------------------------------------------------------------------
    !
    !
    !                              Operator Implementations
    !
    !-----------------------------------------------------------------------------------------



    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-----------------------------------------------------------------------------------------
    function mult_real_chidg_vector(left,right) result(res)
        real(rk),               intent(in)  :: left
        type(chidg_vector_t),    intent(in)  :: right

        type(chidg_vector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(right%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(right%dom)
            res%dom(idom) = left * right%dom(idom)
        end do

        res%send = right%send
        res%recv = right%recv

    end function mult_real_chidg_vector
    !*****************************************************************************************




    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function mult_chidg_vector_real(left,right) result(res)
        type(chidg_vector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right

        type(chidg_vector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(left%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(left%dom)
            res%dom(idom) = left%dom(idom) * right
        end do

        res%send = left%send
        res%recv = left%recv

    end function mult_chidg_vector_real
    !****************************************************************************************





    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function div_real_chidg_vector(left,right) result(res)
        real(rk),               intent(in)  :: left
        type(chidg_vector_t),    intent(in)  :: right

        type(chidg_vector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(right%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(right%dom)
            res%dom(idom) = left / right%dom(idom)
        end do


        res%send = right%send
        res%recv = right%recv

    end function div_real_chidg_vector
    !****************************************************************************************




    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function div_chidg_vector_real(left,right) result(res)
        type(chidg_vector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right

        type(chidg_vector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(left%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(left%dom)
            res%dom(idom) = left%dom(idom) / right
        end do

        res%send = left%send
        res%recv = left%recv

    end function div_chidg_vector_real
    !****************************************************************************************





    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function add_chidg_vector_chidg_vector(left,right) result(res)
        type(chidg_vector_t),    intent(in)  :: left
        type(chidg_vector_t),    intent(in)  :: right

        type(chidg_vector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(right%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(left%dom)
            res%dom(idom) = left%dom(idom) + right%dom(idom)
        end do


        res%send = right%send
        res%recv = right%recv

    end function add_chidg_vector_chidg_vector
    !****************************************************************************************





    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function sub_chidg_vector_chidg_vector(left,right) result(res)
        type(chidg_vector_t),    intent(in)  :: left
        type(chidg_vector_t),    intent(in)  :: right

        type(chidg_vector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(right%dom)

        allocate(res%dom(ndom))


        do idom = 1,size(left%dom)
            res%dom(idom) = left%dom(idom) - right%dom(idom)
        end do

        res%send = right%send
        res%recv = right%recv

    end function sub_chidg_vector_chidg_vector
    !****************************************************************************************












end module type_chidg_vector
